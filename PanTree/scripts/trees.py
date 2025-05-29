#!/usr/bin/env python3

import os
import logging
import pantreelib as pl
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    if not os.path.exists(args.input_vcf):
        logger.error(f"Couldn't find input VCF '{args.input_vcf}. Exiting.'")
        exit()

    # Process vcf then filter to the missing data threshold
    vcfReport = pl.readVCF(args.input_vcf)
    filteredVCF = pl.filterByMissing(vcfReport[0], args.missing)
    # Process vcf and window across chromosomes, finding alignments that match the thresholds
    if args.distance:
        if args.reference:
            fullNames = list(filteredVCF.columns[8:])
            fullNames[0] = 'Reference'
        else:
            fullNames = list(filteredVCF.columns[9:])
        sampleNames = fullNames if args.distance == "All" else [a.strip() for a in args.distance.split(',') if
                                                                a in fullNames]
        pairingList = [(a.strip(), b.strip()) for index, a in enumerate(sampleNames) for b in sampleNames[index + 1:]]
    else:
        pairingList = []
        sampleNames = []

    with open(args.output_tsv, 'w') as outFile:
        totalAlign = pl.initAlignment(filteredVCF, ref=args.reference)

        currentRow = f"Type\tChromosome\tStart\tEnd\tBP\tNewick"
        if args.distance:
            currentRow += "\t" + "\t".join([f"{a[0]}-{a[1]}" for a in pairingList])
        outFile.write(currentRow + "\n")

        for name, group in filteredVCF.groupby('#CHROM'):
            schemes = pl.windowScheme(len(group), args.size, args.spacing)
            group = group.reset_index().drop(columns='index')
            if args.compare:
                chromosomeAlign = pl.pdToAlignment(group, args.reference)
                totalAlign += chromosomeAlign
                chromosomeTree = pl.buildTree(chromosomeAlign)
                chromosomeBranchLength = chromosomeTree.total_branch_length()
                currentRow = f"Chromosome\t{name}\t1\t{len(chromosomeAlign[0])}\t{sum([len(a.seq.replace('-', '')) for a in chromosomeAlign])}\t{pl.beautifyTree(chromosomeTree)}"
                if args.distance:
                    currentRow += "\t" + "\t".join(
                        [str(chromosomeTree.distance(a[0], a[1]) / chromosomeBranchLength) for a in pairingList])
                outFile.write(currentRow + "\n")

            for scheme in schemes:
                schemeAlign = pl.pdToAlignment(group[scheme[0]:scheme[1]], args.reference)
                schemeTree = pl.buildTree(schemeAlign)
                schemeBranchLength = schemeTree.total_branch_length()

                if args.aberrant:
                    if len([a.branch_length / schemeBranchLength for a in schemeTree.depths() if
                            a.branch_length / schemeBranchLength >= args.aberrant]) > 0:
                        pass

                currentRow = f"Window\t{name}\t{group.POS[int(scheme[0])]}\t{group.POS[int(scheme[1]) - 1]}\t{sum([len(a.seq.replace('-', '')) for a in schemeAlign])}\t{pl.beautifyTree(schemeTree)}"
                if args.distance:
                    currentRow += "\t" + "\t".join([str(schemeTree.distance(a[0], a[1])/schemeBranchLength) for a in pairingList])
                outFile.write(currentRow+"\n")

        if args.compare:
            totalTree = pl.buildTree(totalAlign)
            totalBranchLength = totalTree.total_branch_length()
            currentRow = f"Full\tAll\t1\t{len(totalAlign[0])}\t{sum([len(a.seq.replace('-', '')) for a in totalAlign])}\t{pl.beautifyTree(totalTree)}"
            if args.distance:
                currentRow += "\t" + "\t".join(
                    [str(totalTree.distance(a[0], a[1]) / totalBranchLength) for a in pairingList])
            outFile.write(currentRow + "\n")
