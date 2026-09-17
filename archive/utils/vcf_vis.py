#!/usr/bin/env python3
import argparse
import os
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm


def parse_valid_pairs(columns, target_samples=None, focus_samples=None):
    """
    Extracts valid pair strings from DataFrame columns.
    target_samples (-s): BOTH samples must be in this list.
    focus_samples (-S): AT LEAST ONE sample must be in this list.
    """
    all_pair_strs = [c[:-5] for c in columns if c.endswith('_tree')]
    valid_pairs = []

    target_set = set(target_samples) if target_samples else None
    focus_set = set(focus_samples) if focus_samples else None

    # If neither filter is applied, return all available pairs
    if not target_set and not focus_set:
        return all_pair_strs

    for pair_str in all_pair_strs:
        parts = pair_str.split('-')
        is_valid_pair = False

        # Test all possible split combinations for hyphenated names
        for i in range(1, len(parts)):
            s1 = "-".join(parts[:i])
            s2 = "-".join(parts[i:])

            # Condition 1 (-s): Are both samples in the target pool?
            in_pool = True
            if target_set:
                if s1 not in target_set or s2 not in target_set:
                    in_pool = False

            # Condition 2 (-S): Is at least one sample in the focus set?
            in_focus = True
            if focus_set:
                if s1 not in focus_set and s2 not in focus_set:
                    in_focus = False

            # If a split satisfies all provided conditions, keep this pair
            if in_pool and in_focus:
                is_valid_pair = True
                break

        if is_valid_pair:
            valid_pairs.append(pair_str)

    return valid_pairs

def plot(chrom, measure, label, outformat, frac, outdir, lowess, valid_pairs, df_chrom, x_all, density):
    fig, ax1 = plt.subplots(figsize=(14, 6))

    # --- 1. Background Variant Density (Right Y-Axis) ---
    ax2 = ax1.twinx()
    ax2.fill_between(x_all, 0, density, color='gray', alpha=0.2, step='mid', label='Variant Density')
    ax2.set_ylabel("Variant Count per Window", color='dimgray', fontsize=11, fontweight='bold')
    ax2.tick_params(axis='y', labelcolor='dimgray')
    ax2.set_ylim(bottom=0)

    # Fix Z-ordering so the density is drawn behind the lines and scatter points
    ax1.set_zorder(ax2.get_zorder() + 1)
    ax1.patch.set_visible(False)

    # --- 2. Foreground Distance Measures (Left Y-Axis) ---
    color_cycle = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    for pair in valid_pairs:
        col_name = f"{pair}_tree" if measure == "tree_norm" else f"{pair}_{measure}"
        # Make sure the column actually exists in case the original TSV didn't include it
        if col_name not in df_chrom.columns:
            continue

        x = df_chrom['Position_Mb'].values
        if measure == "tree_norm":
            y = df_chrom[col_name].values / df_chrom['Total_Branch_Length']
        else:
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
            smoothed = lowess(y_valid, x_valid, frac=frac)
            ax1.plot(smoothed[:, 0], smoothed[:, 1], color=color, linewidth=2.5, label=pair)
        else:
            ax1.plot(x_valid, y_valid, color=color, linewidth=2.5, label=pair)

    # Aesthetics for Primary Axis
    ax1.set_title(f"Chromosome: {chrom} | {label}")
    ax1.set_xlabel("Genomic Position (Mb)")
    ax1.set_ylabel("Distance", fontsize=11, fontweight='bold')
    ax1.grid(True, linestyle='--', alpha=0.4)

    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    from matplotlib.patches import Patch
    lines1.append(Patch(facecolor='gray', alpha=0.2))
    labels1.append('Variant Density')

    ax1.legend(lines1, labels1, bbox_to_anchor=(1.08, 1), loc='upper left', borderaxespad=0., fontsize='small')

    plt.tight_layout()

    out_file = os.path.join(outdir, f"{chrom}_{measure}.{outformat}")
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description="Visualize sliding window distances with Scatter, LOESS, and Variant Density.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file from distance script.")
    parser.add_argument("-s", "--samples", type=str,
                        help="Comma-separated pool of samples. Only pairs between these samples are plotted (Both must match).")
    parser.add_argument("-S", "--focus-samples", dest="focus", type=str,
                        help="Comma-separated list of focus samples. Only pairs containing at least one of these are plotted.")
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

    # Process filtering arguments
    target_samples = args.samples.split(',') if args.samples else None
    focus_samples = args.focus.split(',') if args.focus else None

    valid_pairs = parse_valid_pairs(df.columns, target_samples, focus_samples)

    if not valid_pairs:
        print("Error: No valid sample pairs found to plot based on provided -s and -S filters.")
        return

    print(f"Found {len(valid_pairs)} sample pair(s) to plot after filtering.")

    chromosomes = df['Chromosome'].unique()
    measures = ['tree', 'tree_norm', 'jaccard', 'd']
    measure_labels = {
        'tree': 'Patristic Distance (NJ Tree)',
        'tree_norm': 'Patristic Distance (NJ Tree) normalized by total branch length',
        'jaccard': 'Jaccard Distance',
        'd': 'Pairwise Divergence (d_xy)'
    }

    lowess = sm.nonparametric.lowess

    for chrom in chromosomes:
        df_chrom = df[df['Chromosome'] == chrom]
        print(f"Plotting {chrom}...")

        # Extract density data for the background
        x_all = df_chrom['Position_Mb'].values
        density = df_chrom['Variant_count'].fillna(0).values

        for measure in measures:
            plot(chrom, measure, measure_labels[measure], args.format, args.frac, args.outdir, lowess, valid_pairs, df_chrom, x_all, density)

    print(f"All plots saved to {args.outdir}/")


if __name__ == "__main__":
    main()