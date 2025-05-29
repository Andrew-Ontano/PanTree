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
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

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
    #missingThreshold = len(tempVcfReport.columns[9:]) - missingAllowed
    missingCounts = (tempVcfReport[tempVcfReport.columns[firstIndex:]] == ".").sum(axis=1)
    return tempVcfReport[missingCounts <= missingAllowed]
    goodRows = []
    for index, row in tempVcfReport.iterrows():

        rowInfo = ast.literal_eval(
            re.sub('"', '', str([re.sub("^(.+)=(.+)$", r"'\1':'\2'", a) for a in row['INFO'].split(";")])).replace('[', '{').replace(
                ']', '}'))
        if int(rowInfo['NS']) >= missingThreshold:
            goodRows.append(index)
    tempVcfReport = tempVcfReport[tempVcfReport.index.isin(goodRows)]
    return tempVcfReport

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
    if ref:
        refSeq = ""
        for row, record in df.iterrows():
            alleles = [record['REF']] + record['ALT'].split(',')
            for column in columnNames:
                allele = alleles[int(record[column])] if record[column] != '.' else '-'
                nucleotideDict[column] += allele
            refSeq += alleles[0]
        nucleotideDict["Reference"] = refSeq
        for a in range(len(alignment)):
            alignment[a].seq = Seq(nucleotideDict[alignment[a].name])
    else:
        for row, record in df.iterrows():
            alleles = [record['REF']] + record['ALT'].split(',')
            for column in columnNames:
                allele = alleles[int(record[column])] if record[column] != '.' else '-'
                nucleotideDict[column] += allele

        for a in range(len(alignment)):
            alignment[a].seq = Seq(nucleotideDict[alignment[a].name])
    value = pd.Series()
    return alignment

# Function: buildTree
# Constructs a phylogenetic tree from alignment
# Input: Biopython MultipleSeqAlignment
# Returns: Biopython DistanceTree
# TODO: modify method to accept multiple tree construction methods
def buildTree(alignment):
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator, 'nj')
    return constructor.build_tree(alignment)

# Function: beautifyTree
# Given a Biopython tree format, returns a human-readable Newick tree
# Input: Biopython DistanceTree
# Returns: String
def beautifyTree(tree):
    return re.sub(r'Inner\d+', '', str(tree.format(fmt='newick')).strip())
