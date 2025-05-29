#!/usr/bin/env python3

import os
import logging
import pantreelib as pl
import re
import math
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    if not os.path.exists(args.input_vcf):
        logger.error(f"Couldn't find input VCF '{args.input_vcf}. Exiting.'")
        exit()

    # Load VCF
    vcfReport = pl.readVCF(args.input_vcf)
    if args.reference:
        sampleNames = list(vcfReport[0].columns[8:])
        sampleNames[0] = 'Reference'
    else:
        sampleNames = list(vcfReport[0].columns[9:])

    # Setup counters for total and baseline
    contigNames = vcfReport[0]['#CHROM'].unique()
    contigs = {a: [int((b.split("="))[-1][:-1]) for b in vcfReport[1] if b.startswith(f"##contig=<ID={a}")][0] for a in contigNames}
    outputDict = {a: {b: {c: 0 for c in sampleNames} for b in list(range(0, contigs[a], args.window))} for a in contigs}
    baselineDict = {a: {b: 0 for b in list(range(0, contigs[a], args.window))} for a in contigs}
    filteredVCF = pl.filterByMissing(vcfReport[0], args.missing)

    if args.groups is not None:
        groups = [group.split(',') for group in re.findall(r'\((.*?)\)', args.groups)]
    else:
        for name, group in filteredVCF.groupby('#CHROM'):
            group = group.reset_index().drop(columns='index')
            if args.reference:
                for rowNum, rowReport in group.iterrows():
                    window = int(math.floor(int(rowReport["POS"]) / args.window)) * args.window
                    baselineDict[name][window] += 1
                    if list(rowReport[9:]).count("0") == 0:
                        outputDict[name][window]["Reference"] += 1
                    else:
                        for sample in sampleNames[1:]:
                            if ([0]+list(rowReport[9:])).count(rowReport[sample]) == 1 and rowReport[sample] != ".":
                                outputDict[name][window][sample] += 1
            else:
                for rowNum, rowReport in group.iterrows():
                    window = int(math.floor(int(rowReport["POS"]) / args.window)) * args.window
                    baselineDict[name][window] += 1

                    for sample in sampleNames:
                        if list(rowReport[9:]).count(rowReport[sample]) == 1 and rowReport[sample] != ".":
                            outputDict[name][window][sample] += 1

    with open(f"{args.output_vcf}.tsv", 'w') as output:
        header = '\t'.join([a for a in sampleNames])
        output.write(f"Chrom\tWindow\t{header}\n")
        for chrom in outputDict.keys():
            for window in outputDict[chrom].keys():
                tallies = "\t".join([str(outputDict[chrom][window][a]) for a in outputDict[chrom][window].keys()])
                output.write(f"{chrom}\t{window}\t{tallies}\n")
    with open(f"{args.output_vcf}.index", "w") as output:
        output.write(f"Chrom\tWindow\tBaseline\n")
        for chrom in baselineDict.keys():
            for window in baselineDict[chrom].keys():
                output.write(f"{chrom}\t{window}\t{str(baselineDict[chrom][window])}\n")
            # window = int(math.floor(int(filteredVCF["POS"]) / args.window)) * args.window


#    outputDict = {a: {b: {c: 0 for c in sampleNames} for b in list(range(0, vcfReport.contigs[a].length, args.window))}
#                  for a in vcfReport.contigs}
#    baselineDict = {a: {b: 0 for b in list(range(0, vcfReport.contigs[a].length, args.window))} for a in
#                    vcfReport.contigs}
