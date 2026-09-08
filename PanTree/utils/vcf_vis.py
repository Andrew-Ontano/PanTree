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
    parser = argparse.ArgumentParser(
        description="Visualize sliding window distances with Scatter + LOESS and Variant Density.")
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

    # Calculate position and sort to ensure fill_between draws correctly
    df['Position_Mb'] = ((df['Window_start'] + df['Window_end']) / 2) / 1_000_000
    df = df.sort_values(by=['Chromosome', 'Position_Mb'])

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

    lowess = sm.nonparametric.lowess

    for chrom in chromosomes:
        df_chrom = df[df['Chromosome'] == chrom]
        print(f"Plotting {chrom}...")

        # Extract density data for the background
        x_all = df_chrom['Position_Mb'].values
        # Fill missing variant counts with 0 just in case
        density = df_chrom['Variant_count'].fillna(0).values

        for measure in measures:
            fig, ax1 = plt.subplots(figsize=(14, 6))

            # --- 1. Background Variant Density (Right Y-Axis) ---
            ax2 = ax1.twinx()
            # step='mid' keeps the fill blocky like genomic windows, rather than a sloped line
            ax2.fill_between(x_all, 0, density, color='gray', alpha=0.2, step='mid', label='Variant Density')
            ax2.set_ylabel("Variant Count per Window", color='dimgray', fontsize=11, fontweight='bold')
            ax2.tick_params(axis='y', labelcolor='dimgray')
            ax2.set_ylim(bottom=0)

            # Fix Z-ordering so the density is drawn *behind* the lines and scatter points
            ax1.set_zorder(ax2.get_zorder() + 1)
            ax1.patch.set_visible(False)  # Makes ax1 transparent so ax2 is visible beneath it

            # --- 2. Foreground Distance Measures (Left Y-Axis) ---
            color_cycle = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])

            for pair in valid_pairs:
                col_name = f"{pair}_{measure}"

                x = df_chrom['Position_Mb'].values
                y = df_chrom[col_name].values
                mask = ~np.isnan(x) & ~np.isnan(y)

                x_valid = x[mask]
                y_valid = y[mask]

                if len(x_valid) == 0:
                    continue

                color = next(color_cycle)

                # Scatter
                ax1.scatter(x_valid, y_valid, color=color, alpha=0.2, s=15, edgecolors='none')

                # LOESS
                if len(x_valid) > 5:
                    smoothed = lowess(y_valid, x_valid, frac=args.frac)
                    ax1.plot(smoothed[:, 0], smoothed[:, 1], color=color, linewidth=2.5, label=pair)
                else:
                    ax1.plot(x_valid, y_valid, color=color, linewidth=2.5, label=pair)

            # Aesthetics for Primary Axis
            ax1.set_title(f"Chromosome: {chrom} | {measure_labels[measure]}")
            ax1.set_xlabel("Genomic Position (Mb)")
            ax1.set_ylabel("Distance", fontsize=11, fontweight='bold')
            ax1.grid(True, linestyle='--', alpha=0.4)

            # Combine legends from both axes
            lines1, labels1 = ax1.get_legend_handles_labels()
            # Add an empty proxy artist for spacing in the legend, then the density label
            from matplotlib.patches import Patch
            lines1.append(Patch(facecolor='gray', alpha=0.2))
            labels1.append('Variant Density')

            # Move the bounding box slightly further right (1.08) to clear the right Y-axis text
            ax1.legend(lines1, labels1, bbox_to_anchor=(1.08, 1), loc='upper left', borderaxespad=0., fontsize='small')

            plt.tight_layout()

            out_file = os.path.join(args.outdir, f"{chrom}_{measure}.{args.format}")
            plt.savefig(out_file, dpi=300, bbox_inches="tight")
            plt.close()

    print(f"All plots saved to {args.outdir}/")


if __name__ == "__main__":
    main()