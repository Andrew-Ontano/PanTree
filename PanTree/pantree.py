#!/usr/bin/env python3

# this is the script to house all components
# Core imports:
import argparse
import importlib

def moduleFromPath(path):
    return path[1:].replace(".py", "").replace("/", ".")


parser = argparse.ArgumentParser(description='PanTree tool for analyzing pangenomic VCFs.')
subparsers = parser.add_subparsers(title='pantree', dest='subparsers', required=True)

parserTrees = subparsers.add_parser("trees", description="Construct phylogenetic trees from traversals across variants")
parserTrees.add_argument('-i', '--input', dest='input_vcf', type=str, help='Input VCF file', required=True)
parserTrees.add_argument('-o', '--output', dest='output_tsv', type=str, help='Output alignment prefix', default='panTree-tree.tsv')
parserTrees.add_argument('-w', '--window-size', dest='size', type=int, help='Window size in # of VCF records to generate alignments from', default=10000)
parserTrees.add_argument('-s', '--window-overlap', dest='spacing', type=int, help='Window overlap in # of VCF records between each window. Can be negative for spaced windows', default=0)
parserTrees.add_argument('-m', '--missing-allowed', dest='missing', type=int, help="Maximum number of missing alleles per row", default=0)
parserTrees.add_argument('-b', '--blacklist', dest='blacklist', type=str, help='Remove species from analysis. Separated by commas', default=None)
parserTrees.add_argument('-c', '--compare', dest='compare', action='store_true', help='Perform comparison between species with windowed trees, then report distance from species tree', default=False)
parserTrees.add_argument('-r', '--reference', dest='reference', action='store_true', help='Incorporate reference lineage into alignments', default=False)
parserTrees.add_argument('-t', '--tree-file', dest='tree', type=str, help="Species tree in newick format. If absent, species tree is generated from VCF", default=None)
parserTrees.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
parserTrees.add_argument('-d', '--distance', dest='distance', type=str, help='Comma-separated list for genotypes for tip-distances. Calculate permutations between these samples. Use "All" to calculate all possible permutations.', default=None)
parserTrees.add_argument('-a', '--aberrant_threshold', dest='aberrant', type=float, help='Weighted branch length threshold to call a branch aberrant. If not in distance mode, prints windows with aberrant branches', default=None)

parserHeatmap = subparsers.add_parser("heatmap", description="Generate heatmap of variant identity")
parserHeatmap.add_argument('-i', '--input', dest='input_vcf', type=str, help='Input VCF file', required=True)
parserHeatmap.add_argument('-o', '--output', dest='output_vcf', type=str, help='Output VCF file', required=True)
parserHeatmap.add_argument('-w', '--window-size', dest='window', type=int, help='Window size in # of VCF records to generate alignments from', default=10000)
parserHeatmap.add_argument('-s', '--window-overlap', dest='spacing', type=int, help='Window overlap in # of VCF records between each window. Can be negative for spaced windows', default=0)
parserHeatmap.add_argument('-m', '--missing-allowed', dest='missing', type=int, help="Maximum number of missing alleles per row", default=0)
parserHeatmap.add_argument('-r', '--reference', dest='reference', action='store_true', help='Incorporate reference lineage into alignments', default=False)
parserHeatmap.add_argument('-g', '--groups', dest='groups', type=str, help="Sample groups to treat as same samples. Format: (A,B),(C,D),(E)", default=None)
#parser.add_argument('-t', '--targets', dest='target_genotypes', type=str, help='Genotype indices, separated by comma', nargs='*', required=True)
parserHeatmap.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)

parserNearest = subparsers.add_parser("nearest", description="Collate alleles that are polarizable between multiple samples")
parserNearest.add_argument('-i', '--input', dest='input_vcf', type=str, help='Input VCF file', required=True)
parserNearest.add_argument('-o', '--output', dest='output_tsv', type=str, help='Output alignment prefix', default='panTree-nearest.tsv')
parserNearest.add_argument('-t', '--targets', dest='targets', type=str, help='Target Samples/Reference to polarize against Format: Name1,Name2,...', required=True)
parserNearest.add_argument('-w', '--window-size', dest='size', type=int, help='Window size in # of VCF records to generate alignments from', default=10000)
parserNearest.add_argument('-s', '--window-overlap', dest='spacing', type=int, help='Window overlap in # of VCF records between each window. Can be negative for spaced windows', default=0)
parserNearest.add_argument('-m', '--missing-allowed', dest='missing', type=int, help="Maximum number of missing alleles per row", default=0)
parserNearest.add_argument('-r', '--reference', dest='reference', action='store_true', help='Incorporate reference lineage into alignments', default=False)
parserNearest.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)


args = parser.parse_args()

subparsersScript = f"/scripts/{args.subparsers}.py"
script = importlib.import_module(moduleFromPath(subparsersScript))
script.execute(args)
