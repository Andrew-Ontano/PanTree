#!/usr/bin/env python3

import os
import logging
import pandas as pd
import PanTree.pantreelib as pl
import re
import math

logger = logging.getLogger(__name__)

def execute(args):
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

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

    filteredVCF = pl.filterByMissing(vcfReport[0], args.missing)

    targetList = [a.strip() for a in args.targets.split(",")]
#    for record in filteredVCF:
#       pass
    if len([a for a in targetList if a in sampleNames]) != len(targetList):
        logger.error(f"Target names ({','.join(targetList)}) not found in samples ({','.join(sampleNames)}) . Exiting.")
        exit()

    with open(args.output_tsv, 'w') as outFile:
        header = "Chrom\tStart\tEnd\t" + "\t".join(targetList) + "\t" + "\t".join([s for s in sampleNames if s not in targetList])
        outFile.write(header + "\n")

        for name, group in filteredVCF.groupby('#CHROM'):
            schemes = pl.windowScheme(len(group), args.size, args.spacing)
            group = group.reset_index().drop(columns='index')

            for scheme in schemes:
                window = group[scheme[0]:scheme[1]]
                # Initialize tallies for this window
                # tallies[sample][target_allele_index]
                otherSamples = [s for s in sampleNames if s not in targetList]
                tallies = {s: [0] * len(targetList) for s in otherSamples}

                informativeSites = 0
                for _, row in window.iterrows():
                    # Get alleles for targets
                    targetAlleles = []
                    for t in targetList:
                        if t == "Reference":
                            targetAlleles.append("0")
                        else:
                            targetAlleles.append(row[t])

                    # Check if all target alleles are present and distinct
                    if "." in targetAlleles or len(set(targetAlleles)) != len(targetList):
                        continue

                    informativeSites += 1
                    for s in otherSamples:
                        allele = row[s]
                        if allele in targetAlleles:
                            tallies[s][targetAlleles.index(allele)] += 1

                if informativeSites > 0:
                    start_pos = group.POS[int(scheme[0])]
                    end_pos = group.POS[int(scheme[1]) - 1]
                    row_out = f"{name}\t{start_pos}\t{end_pos}"
                    # We can't really report target alleles per window if they change,
                    # but the requirement was "report how often other samples match each target allele"
                    # Maybe report percentages?

                    # Add dummy values for targets to match header
                    for t in targetList:
                        row_out += "\t-"

                    for s in otherSamples:
                        row_out += "\t" + ",".join([str(count) for count in tallies[s]])

                    outFile.write(row_out + "\n")