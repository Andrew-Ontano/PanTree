#!/usr/bin/env python3
import argparse
import sys
import logging
from scripts import calculate, visualize


def main():
    parser = argparse.ArgumentParser(description="PanTree tool for analyzing pangenomic VCFs.")
    subparsers = parser.add_subparsers(title="subcommands", dest="subcommand", required=True)

    # Subcommand: calculate
    parser_calc = subparsers.add_parser(
        "calculate",
        description="Calculate sliding window stats, sequence alignments, and phylogenetic distances from VCF."
    )
    parser_calc.add_argument("-i", "--input", dest="input_vcf", required=True, help="Input VCF / VCF.gz (tabix-indexed)")
    parser_calc.add_argument("-o", "--output", dest="output_tsv", required=True, help="Output TSV file")
    parser_calc.add_argument("-s", "--samples", dest="samples", help="Sample list for comparisons (comma-separated)")
    parser_calc.add_argument("-w", "--window-size", dest="window_size", type=int, default=10000, help="Window bp or variant width")
    parser_calc.add_argument("-l", "--overlap", dest="overlap", type=int, default=0, help="Window bp or variant overlap")
    parser_calc.add_argument("-f", "--fixed-variants", dest="fixed_variants", action="store_true", help="Fixed-variant width mode instead of genomic position mode")
    parser_calc.add_argument("-r", "--include-ref", dest="include_ref", action="store_true", help="Should reference sequence be incorporated in unique counts/alignments")
    parser_calc.add_argument("-t", "--threads", dest="threads", type=int, default=1, help="Computing threads")
    parser_calc.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="Verbose step reporting mode")

    # Subcommand: visualize
    parser_vis = subparsers.add_parser(
        "visualize",
        description="Visualize window statistics and qualities."
    )
    parser_vis.add_argument("-i", "--input", dest="input_tsv", required=True, help="Input TSV file from calculate subroutine")
    parser_vis.add_argument("-o", "--output-dir", dest="outdir", default="plots", help="Output directory for plots")
    parser_vis.add_argument("-s", "--samples", dest="samples", help="Comma-separated pool of samples. Only pairs between these samples are plotted (both must match).")
    parser_vis.add_argument("-S", "--focus-samples", dest="focus", help="Comma-separated focus samples. Only pairs containing at least one are plotted.")
    parser_vis.add_argument("-f", "--frac", dest="frac", type=float, default=0.1, help="LOESS smoothing fraction (0.0 to 1.0, default: 0.1)")
    parser_vis.add_argument("--format", dest="format", default="png", choices=["png", "pdf", "svg"], help="Output image format (png, pdf, svg)")
    parser_vis.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="Verbose step reporting mode")

    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format="%(asctime)s [%(levelname)s] %(message)s")

    if args.subcommand == "calculate":
        calculate.run(args)
    elif args.subcommand == "visualize":
        visualize.run(args)


if __name__ == "__main__":
    main()
