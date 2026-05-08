#!/usr/bin/env python3

import os
import logging
import PanTree.pantreelib as pl

logger = logging.getLogger(__name__)

def execute(args):
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

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

    def get_metrics(alignment, pairingList, args):
        results = []
        if args.metric == 'tree':
            tree = pl.buildTree(alignment, method=args.method)
            totalBranchLength = tree.total_branch_length()
            newick = pl.beautifyTree(tree)
            distances = []
            if args.distance:
                for pair in pairingList:
                    try:
                        d = tree.distance(pair[0], pair[1]) / totalBranchLength if totalBranchLength > 0 else 0.0
                    except:
                        d = "-"
                    distances.append(str(d))
            return newick, distances
        elif args.metric == 'p-distance':
            dm = pl.calculatePDistance(alignment)
            distances = []
            if args.distance:
                for pair in pairingList:
                    distances.append(str(dm[pair[0], pair[1]]))
            return "-", distances
        elif args.metric == 'jaccard':
            dm = pl.calculateJaccardDistance(alignment)
            distances = []
            if args.distance:
                for pair in pairingList:
                    distances.append(str(dm[pair[0]][pair[1]]))
            return "-", distances
        elif args.metric == 'd-stat':
            distances = []
            if args.quadruple:
                p1, p2, p3, o = [x.strip() for x in args.quadruple.split(',')]
                d = pl.calculateDStatistic(alignment, p1, p2, p3, o)
                distances.append(str(d))
            return "-", distances
        return "-", []

    with open(args.output_tsv, 'w') as outFile:
        totalAlign = pl.initAlignment(filteredVCF, ref=args.reference)

        currentRow = f"Type\tChromosome\tStart\tEnd\tBP\tNewick"
        if args.metric == 'd-stat':
            currentRow += f"\tD-stat({args.quadruple})"
        elif args.distance:
            currentRow += "\t" + "\t".join([f"{a[0]}-{a[1]}" for a in pairingList])
        outFile.write(currentRow + "\n")

        for name, group in filteredVCF.groupby('#CHROM'):
            schemes = pl.windowScheme(len(group), args.size, args.spacing)
            group = group.reset_index().drop(columns='index')
            if args.compare:
                chromosomeAlign = pl.pdToAlignment(group, args.reference)
                totalAlign += chromosomeAlign
                newick, distances = get_metrics(chromosomeAlign, pairingList, args)
                currentRow = f"Chromosome\t{name}\t1\t{len(chromosomeAlign[0])}\t{sum([len(a.seq.replace('-', '')) for a in chromosomeAlign])}\t{newick}"
                if distances:
                    currentRow += "\t" + "\t".join(distances)
                outFile.write(currentRow + "\n")

            for scheme in schemes:
                schemeAlign = pl.pdToAlignment(group[scheme[0]:scheme[1]], args.reference)
                newick, distances = get_metrics(schemeAlign, pairingList, args)

                currentRow = f"Window\t{name}\t{group.POS[int(scheme[0])]}\t{group.POS[int(scheme[1]) - 1]}\t{sum([len(a.seq.replace('-', '')) for a in schemeAlign])}\t{newick}"
                if distances:
                    currentRow += "\t" + "\t".join(distances)
                outFile.write(currentRow+"\n")

        if args.compare:
            totalTree_newick, total_distances = get_metrics(totalAlign, pairingList, args)
            currentRow = f"Full\tAll\t1\t{len(totalAlign[0])}\t{sum([len(a.seq.replace('-', '')) for a in totalAlign])}\t{totalTree_newick}"
            if total_distances:
                currentRow += "\t" + "\t".join(total_distances)
            outFile.write(currentRow + "\n")
