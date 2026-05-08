#!/usr/bin/env python3

# Houses core functions needed across pantree

import logging
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import gzip
import re
import ast

# Set up logger
logger = logging.getLogger(__name__)

# Function: readVCF
# Reads VCF file or gzipped file (must end in .gz), processes headers, then transforms data to pandas dataframe
# Input: vcf file path
# Returns pandas dataframe
def readVCF(vcfFile):
    comments = []
    if vcfFile.endswith('.gz'):
        with gzip.open(vcfFile, 'rt') as inFile:
            for line in inFile:
                comments.append(line.strip())
                if line.startswith('#CHROM'):
                    header = [x.strip() for x in line.split('\t')]
                    break
    else:
        with open(vcfFile, 'r') as inFile:
            for line in inFile:
                comments.append(line.strip())
                if line.startswith('#CHROM'):
                    header = [x.strip() for x in line.split('\t')]
                    break

    return pd.read_csv(vcfFile, dtype=str, sep='\t', header=None, names=header, comment="#"), comments

# Function: filterByMissing
# Removes unwanted rows from vcf dataframe with too many missing allele calls.
# Input: pandas dataframe
# Returns: pandas dataframe
def filterByMissing(tempVcfReport, missingAllowed=0, firstIndex=9):
    def is_missing(val):
        return val == "." or (isinstance(val, str) and "." in val.split("/")) or (isinstance(val, str) and "." in val.split("|"))

    missingCounts = tempVcfReport[tempVcfReport.columns[firstIndex:]].apply(lambda col: col.map(is_missing)).sum(axis=1)
    return tempVcfReport[missingCounts <= missingAllowed]

# Function: filterByMType
# Removes unwanted rows from vcf dataframe based on variant type. Defaults to biallelic snps
# Input: pandas dataframe
# Returns: pandas dataframe
# TODO: EXPERIMENTAL - work out in-program variant filtering. Preprocess is currently required.
def filterByType(tempVcfReport, variantTypes=['snp'], maxAlleles=2, missingAllowed=0):
    filteredVariantTypes = "_".join([a for a in ['snp', 'mnp', 'indel'] if a not in variantTypes])
    missingThreshold = len(tempVcfReport.columns[9:]) - missingAllowed
    typePattern = None
    goodRows = []
    match filteredVariantTypes:
        # only snp, then keep if ref and alt are not both 1 character long
        case "snp":
            pass
        # only mnp, then keep if alt and ref are not the same length longer than 1
        case "mnp":
            pass
        # only indel, then check if length of alt is same as ref
        case "indel":
            pass
        # snp and mnp, then keep if ref and alt are not same length
        case "snp_mnp":
            pass
        # snp and indel, then keep if alt and ref are the same length longer than 1
        case "snp_indel":
            pass
        # indel and mnp, then keep if ref and alt both 1 character long
        case "mnp_indel":
            pass
        case _:
            pass

    for index, row in tempVcfReport.iterrows():
        rowInfo = ast.literal_eval(
            re.sub('"', '', str([re.sub("^(.+)=(.+)$", r"'\1':'\2'", a) for a in row['INFO'].split(";")])).replace('[', '{').replace(
                ']', '}'))
        if int(rowInfo['NS']) >= missingThreshold:
            goodRows.append(index)

    tempVcfReport = tempVcfReport[tempVcfReport.index.isin(goodRows)]
    return tempVcfReport

# Function: windowScheme
# Plans the windowing scheme to be applied on a per-alignment basis
# Input: values based on grouped pandas dataframe and arguments
# Returns: List of tuples with 2 values
def windowScheme(length, size, spacing):
    currentPosition = 0
    outputIndices = []
    while currentPosition <= length:
        outputIndices.append((currentPosition, currentPosition + size if currentPosition + size < length else length))
        currentPosition += spacing + size
        if currentPosition >= length:
            break
    return outputIndices

# Function: initAlignment
# Sets up an alignment based on a given dataframe
# Input: pandas dataframe
# Returns: Biopython MultipleSeqAlignment
def initAlignment(df, ref=False):
    columnNames = df.columns[9:]
    alignment = MultipleSeqAlignment(records=[])
    if ref:
        currentRecord = SeqRecord(Seq(""), id="Reference", name="Reference")
        alignment.append(currentRecord)
    for column in columnNames:
        currentRecord = SeqRecord(Seq(""), id=column, name=column)
        alignment.append(currentRecord)
    return alignment

# Function: pdToAlignment
# Populates an alignment with data from a dataframe
# Input: pandas dataframe
# Returns: Biopython MultipleSeqAlignment
def pdToAlignment(df, ref=False):
    alignment = initAlignment(df, ref)
    columnNames = df.columns[9:]
    nucleotideDict = {a: "" for a in columnNames}

    def get_allele_index(val):
        if val == ".": return None
        if isinstance(val, str):
            if "/" in val:
                parts = val.split("/")
                return int(parts[0]) if parts[0] != "." else None
            if "|" in val:
                parts = val.split("|")
                return int(parts[0]) if parts[0] != "." else None
        try:
            return int(val)
        except (ValueError, TypeError):
            return None

    if ref:
        refSeq = ""
        for row, record in df.iterrows():
            alleles = [record['REF']] + record['ALT'].split(',')
            for column in columnNames:
                idx = get_allele_index(record[column])
                allele = alleles[idx] if idx is not None else '-'
                nucleotideDict[column] += allele
            refSeq += alleles[0]
        nucleotideDict["Reference"] = refSeq
        for a in range(len(alignment)):
            alignment[a].seq = Seq(nucleotideDict[alignment[a].name])
    else:
        for row, record in df.iterrows():
            alleles = [record['REF']] + record['ALT'].split(',')
            for column in columnNames:
                idx = get_allele_index(record[column])
                allele = alleles[idx] if idx is not None else '-'
                nucleotideDict[column] += allele

        for a in range(len(alignment)):
            alignment[a].seq = Seq(nucleotideDict[alignment[a].name])

    return alignment

# Function: buildTree
# Constructs a phylogenetic tree from alignment
# Input: Biopython MultipleSeqAlignment, method ('nj' or 'upgma')
# Returns: Biopython DistanceTree
def buildTree(alignment, method='nj'):
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator, method)
    return constructor.build_tree(alignment)

# Function: calculatePDistance
# Computes raw pairwise genetic distances (p-distance)
# Input: Biopython MultipleSeqAlignment
# Returns: DistanceMatrix
def calculatePDistance(alignment):
    calculator = DistanceCalculator('identity')
    return calculator.get_distance(alignment)

# Function: calculateJaccardDistance
# Computes Jaccard distance based on shared non-reference alleles
# Input: Biopython MultipleSeqAlignment, reference sequence name
# Returns: DistanceMatrix-like dictionary
def calculateJaccardDistance(alignment, refName="Reference"):
    names = [r.id for r in alignment]
    dm = {n1: {n2: 0.0 for n2 in names} for n1 in names}

    # Find reference sequence
    refSeq = None
    for record in alignment:
        if record.id == refName:
            refSeq = str(record.seq)
            break

    if refSeq is None:
        # If no reference, we can't easily define "non-reference"
        # Fallback to 1 - Jaccard similarity of all alleles?
        # Let's just use the first sequence as reference if not found
        refSeq = str(alignment[0].seq)

    num_sites = alignment.get_alignment_length()

    for i in range(len(names)):
        for j in range(i, len(names)):
            n1, n2 = names[i], names[j]
            s1 = str(alignment[i].seq)
            s2 = str(alignment[j].seq)

            intersection = 0
            union = 0
            for k in range(num_sites):
                a1, a2, r = s1[k], s2[k], refSeq[k]
                if a1 == '-' or a2 == '-' or r == '-':
                    continue

                if a1 != r or a2 != r:
                    union += 1
                    if a1 == a2:
                        intersection += 1

            dist = 1.0 - (intersection / union) if union > 0 else 0.0
            dm[n1][n2] = dm[n2][n1] = dist
    return dm

# Function: calculateDStatistic
# Implements the ABBA-BABA test (D-statistic) for a quadruple of lineages
# Input: Biopython MultipleSeqAlignment, names of (P1, P2, P3, O)
# Returns: D-statistic float
def calculateDStatistic(alignment, p1, p2, p3, o):
    seqs = {r.id: str(r.seq) for r in alignment}
    for name in [p1, p2, p3, o]:
        if name not in seqs:
            return None

    n_abba = 0
    n_baba = 0

    for k in range(alignment.get_alignment_length()):
        a1, a2, a3, ao = seqs[p1][k], seqs[p2][k], seqs[p3][k], seqs[o][k]

        if '-' in [a1, a2, a3, ao]:
            continue

        if a1 == ao and a2 == a3 and a1 != a2:
            n_baba += 1
        elif a2 == ao and a1 == a3 and a1 != a2:
            n_abba += 1

    if (n_abba + n_baba) == 0:
        return 0.0

    return (n_abba - n_baba) / (n_abba + n_baba)

# Function: beautifyTree
# Given a Biopython tree format, returns a human-readable Newick tree
# Input: Biopython DistanceTree
# Returns: String
def beautifyTree(tree):
    return re.sub(r'Inner\d+', '', str(tree.format(fmt='newick')).strip())
