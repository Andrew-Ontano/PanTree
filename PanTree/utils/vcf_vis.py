#!/usr/bin/env python3
import argparse
import os
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm


def parse_valid_pairs(columns, target_samples=None):
    """
    Extracts valid pair strings from DataFrame columns, matching exact sample names.
    """
    all_pair_strs = [c[:-5] for c in columns if c.endswith('_tree')]
    valid_pairs = []

    if target_samples:
        target_set = set(target_samples)
        for pair_str in all_pair_strs:
            parts = pair_str.split('-')
            matched = False
            for i in range(1, len(parts)):
                s1 = "-".join(parts[:i])
                s2 = "-".join(parts[i:])

                if s1 in target_set and s2 in target_set:
                    matched = True
                    break

            if matched:
                valid_pairs.append(pair_str)
    else:
        valid_pairs = all_pair_strs

    return valid_pairs


def main():
    parser = argparse.ArgumentParser(description="Visualize sliding window distances with Scatter + LOESS.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file from distance script.")
    parser.add_argument("-s", "--samples", type=str, help="Comma-separated list of samples to plot.")
    parser.add_argument("-f", "--frac", type=float, default=0.1,
                        help="LOESS smoothing fraction (0.0 to 1.0). Default: 0.1")
    parser.add_argument("-o", "--outdir", default="plots", help="Output directory for the plots.")
    parser.add_argument("--format", default="png", choices=['png', 'pdf', 'svg'], help="Output image format.")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading data from {args.input}...")
    df = pd.read_csv(args.input, sep='\t', na_values=["NA", "TREE_ERROR"])
    df['Position_Mb'] = ((df['Window_start'] + df['Window_end']) / 2) / 1_000_000

    target_samples = args.samples.split(',') if args.samples else None
    valid_pairs = parse_valid_pairs(df.columns, target_samples)

    if not valid_pairs:
        print("Error: No valid sample pairs found to plot. Check your --samples names.")
        return

    print(f"Found {len(valid_pairs)} sample pair(s) to plot.")

    chromosomes = df['Chromosome'].unique()
    measures = ['tree', 'jaccard', 'd']
    measure_labels = {
        'tree': 'Patristic Distance (NJ Tree)',
        'jaccard': 'Jaccard Distance',
        'd': 'Pairwise Divergence (d_xy)'
    }

    # Extract the lowess smoothing function
    lowess = sm.nonparametric.lowess

    for chrom in chromosomes:
        df_chrom = df[df['Chromosome'] == chrom]
        print(f"Plotting {chrom}...")

        for measure in measures:
            plt.figure(figsize=(14, 6))

            # Create a robust color iterator from the current style's color cycle
            color_cycle = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])

            for pair in valid_pairs:
                col_name = f"{pair}_{measure}"

                # Extract X and Y, and drop NaNs
                x = df_chrom['Position_Mb'].values
                y = df_chrom[col_name].values
                mask = ~np.isnan(x) & ~np.isnan(y)

                x_valid = x[mask]
                y_valid = y[mask]

                if len(x_valid) == 0:
                    continue

                # Get the next color robustly
                color = next(color_cycle)

                # 1. Plot the raw data as a faded scatter plot
                plt.scatter(x_valid, y_valid, color=color, alpha=0.2, s=15, edgecolors='none')

                # 2. Calculate and plot the LOESS trendline
                if len(x_valid) > 5:
                    smoothed = lowess(y_valid, x_valid, frac=args.frac)
                    plt.plot(smoothed[:, 0], smoothed[:, 1], color=color, linewidth=2.5, label=pair)
                else:
                    plt.plot(x_valid, y_valid, color=color, linewidth=2.5, label=pair)

            plt.title(f"Chromosome: {chrom} | {measure_labels[measure]}")
            plt.xlabel("Genomic Position (Mb)")
            plt.ylabel("Distance")

            plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0., fontsize='small')
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.tight_layout()

            out_file = os.path.join(args.outdir, f"{chrom}_{measure}.{args.format}")
            plt.savefig(out_file, dpi=300, bbox_inches="tight")
            plt.close()

    print(f"All plots saved to {args.outdir}/")


if __name__ == "__main__":
    main()