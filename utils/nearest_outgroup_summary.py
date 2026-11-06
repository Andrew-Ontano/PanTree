#!/usr/bin/env python3
"""
Summarize how frequently selected samples are closest to non-selected samples
(outgroups) across genomic windows using distance metrics from `pantree.py calculate`.
Generates a tabulated TSV summary (per-sequence and pangenome-wide) and a grouped
stacked barplot for a specified distance metric.
"""

import argparse
import logging
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_samples_from_columns(columns):
    """
    Extract sample names from TSV columns matching `{sample}_called_sites`.
    """
    samples = []
    suffix = "_called_sites"
    for col in columns:
        if col.endswith(suffix):
            sample_name = col[:-len(suffix)]
            samples.append(sample_name)
    return samples


def get_pair_column(s1, s2, measure, columns):
    """
    Find column name for distance between s1 and s2 under measure.
    Checks both `{s1}-{s2}_{measure}` and `{s2}-{s1}_{measure}`.
    """
    cand1 = f"{s1}-{s2}_{measure}"
    if cand1 in columns:
        return cand1
    cand2 = f"{s2}-{s1}_{measure}"
    if cand2 in columns:
        return cand2
    return None


def calculate_nearest_outgroup_summary(df, selected_samples, non_selected_samples, measures):
    """
    Computes nearest outgroup metrics across windows per chromosome and across all chromosomes (pangenome).
    Returns a pandas DataFrame summary.
    """
    records = []

    # Regions to process: each unique chromosome + 'PANGENOME' (all chromosomes combined)
    chromosomes = list(df['Chromosome'].unique())
    scopes = [(c, df[df['Chromosome'] == c]) for c in chromosomes]
    scopes.append(('PANGENOME', df))

    for scope_name, sub_df in scopes:
        total_scope_windows = len(sub_df)
        total_scope_variants = pd.to_numeric(sub_df['Total_variants'], errors='coerce').fillna(0).sum()

        for measure in measures:
            for s in selected_samples:
                # Find available outgroups that have pairwise distance columns with s
                valid_outgroups = []
                outgroup_cols = {}
                for outgroup in non_selected_samples:
                    col = get_pair_column(s, outgroup, measure, df.columns)
                    if col:
                        valid_outgroups.append(outgroup)
                        outgroup_cols[outgroup] = col

                if not valid_outgroups:
                    continue

                # Tally windows and variants where outgroup is closest to s
                outgroup_win_counts = {og: 0.0 for og in valid_outgroups}
                outgroup_var_counts = {og: 0.0 for og in valid_outgroups}

                for _, row in sub_df.iterrows():
                    var_count = float(row.get('Total_variants', 0)) if pd.notna(row.get('Total_variants')) else 0.0

                    # Extract distances to all valid outgroups
                    dist_map = {}
                    for og in valid_outgroups:
                        val = row.get(outgroup_cols[og])
                        if pd.notna(val) and val != 'NA' and val != 'TREE_ERROR':
                            try:
                                dist_map[og] = float(val)
                            except (ValueError, TypeError):
                                pass

                    if not dist_map:
                        continue

                    min_dist = min(dist_map.values())
                    # All outgroups sharing min distance (handles ties evenly)
                    closest_ogs = [og for og, d in dist_map.items() if d == min_dist]
                    weight = 1.0 / len(closest_ogs)

                    for og in closest_ogs:
                        outgroup_win_counts[og] += weight
                        outgroup_var_counts[og] += weight * var_count

                for og in valid_outgroups:
                    win_count = outgroup_win_counts[og]
                    var_count = outgroup_var_counts[og]
                    prop_wins = (win_count / total_scope_windows) if total_scope_windows > 0 else 0.0
                    prop_vars = (var_count / total_scope_variants) if total_scope_variants > 0 else 0.0

                    records.append({
                        "Scope": scope_name,
                        "Sample": s,
                        "Outgroup": og,
                        "Measure": measure,
                        "Total_windows_closest": round(win_count, 4),
                        "Total_variants_closest": round(var_count, 4),
                        "Proportion_windows_closest": round(prop_wins, 6),
                        "Proportion_variants_closest": round(prop_vars, 6),
                    })

    summary_df = pd.DataFrame(records)
    return summary_df


def plot_grouped_stacked_bar(summary_df, measure, selected_samples, non_selected_samples, output_plot, value_col="Total_windows_closest"):
    """
    Creates a grouped barplot of stacked barplots.
    - Each group on the x-axis represents a chromosome (excluding 'PANGENOME').
    - Each bar within a group represents a selected sample.
    - Stacked segments represent totals where the selected sample was closest to each non-selected sample.
    """
    chrom_df = summary_df[(summary_df['Scope'] != 'PANGENOME') & (summary_df['Measure'] == measure)]
    if chrom_df.empty:
        logging.warning("No chromosome-level data to plot.")
        return

    chromosomes = list(chrom_df['Scope'].unique())
    num_chroms = len(chromosomes)
    num_samples = len(selected_samples)

    fig, ax = plt.subplots(figsize=(max(10, num_chroms * max(2.5, num_samples * 0.8)), 6))

    # Color mapping for non-selected samples (outgroups)
    cmap = plt.get_cmap('tab10' if len(non_selected_samples) <= 10 else 'tab20')
    outgroup_colors = {og: cmap(i % 20) for i, og in enumerate(non_selected_samples)}

    # Bar positions
    group_width = 0.8
    bar_width = group_width / max(num_samples, 1)
    indices = np.arange(num_chroms)

    legend_handles = {}

    for sample_idx, sample in enumerate(selected_samples):
        # Calculate x offsets for each selected sample within each chromosome group
        x_offsets = indices - (group_width / 2) + (sample_idx + 0.5) * bar_width

        # Track bottom offset for stacked bars
        bottoms = np.zeros(num_chroms)

        for og in non_selected_samples:
            vals = []
            for chrom in chromosomes:
                sub = chrom_df[(chrom_df['Scope'] == chrom) &
                               (chrom_df['Sample'] == sample) &
                               (chrom_df['Outgroup'] == og)]
                val = sub[value_col].sum() if not sub.empty else 0.0
                vals.append(val)
            vals = np.array(vals)

            color = outgroup_colors[og]
            bars = ax.bar(x_offsets, vals, bar_width, bottom=bottoms,
                          color=color, edgecolor='white', linewidth=0.5,
                          label=og if og not in legend_handles else "")
            if og not in legend_handles:
                legend_handles[og] = bars[0]
            bottoms += vals

    # Format Axes
    ax.set_xticks(indices)
    ax.set_xticklabels(chromosomes, rotation=45, ha='right', fontsize=10, fontweight='bold')
    ax.set_xlabel("Chromosome", fontsize=11, fontweight='bold')

    val_label = "Total Windows" if "windows" in value_col.lower() else "Total Variants"
    ax.set_ylabel(f"{val_label} Closest to Outgroup", fontsize=11, fontweight='bold')
    ax.set_title(f"Nearest Outgroup Breakdown by Chromosome ({measure.capitalize()} Distance)", fontsize=13, fontweight='bold')

    # Add sample labels below or above bars if multiple samples
    if num_samples > 1:
        # Custom secondary labels under x ticks or annotations
        sample_patches = [plt.Line2D([0], [0], color='none') for _ in selected_samples]
        # Or add a legend with outgroup colors
        pass

    # Unique legend for outgroups
    ax.legend(legend_handles.values(), legend_handles.keys(), title="Outgroup",
              bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True)

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_plot)), exist_ok=True)
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    plt.close()
    logging.info(f"Summary plot saved to {output_plot}")


def main():
    parser = argparse.ArgumentParser(
        description="Summarize nearest outgroup relationship frequency from PanTree calculate results."
    )
    parser.add_argument("-i", "--input", dest="input_tsv", required=True, help="Input TSV file from pantree.py calculate")
    parser.add_argument("-s", "--samples", dest="samples", required=True, help="Comma-separated list of selected samples to assess")
    parser.add_argument("-o", "--output-tsv", dest="output_tsv", required=True, help="Output TSV file for tabulated summary")
    parser.add_argument("-p", "--output-plot", dest="output_plot", help="Output image file for grouped stacked barplot (e.g. plot.png)")
    parser.add_argument(
        "-m", "--measure", dest="measure", default="manhattan",
        choices=["manhattan", "patristic", "jaccard"],
        help="Distance metric to use for the plot (default: manhattan)"
    )
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="Verbose logging")

    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format="%(asctime)s [%(levelname)s] %(message)s")

    logging.info(f"Reading calculate TSV from {args.input_tsv}...")
    df = pd.read_csv(args.input_tsv, sep='\t', na_values=["NA", "TREE_ERROR"])

    all_samples = parse_samples_from_columns(df.columns)
    if not all_samples:
        # Fallback: check pairwise columns
        for c in df.columns:
            if c.endswith("_manhattan") or c.endswith("_patristic") or c.endswith("_jaccard"):
                base = c.rsplit('_', 1)[0]
                parts = base.split('-')
                for p in parts:
                    if p not in all_samples and p:
                        all_samples.append(p)

    selected_samples = [s.strip() for s in args.samples.split(',') if s.strip()]
    non_selected_samples = [s for s in all_samples if s not in selected_samples]

    logging.info(f"Selected samples ({len(selected_samples)}): {selected_samples}")
    logging.info(f"Non-selected samples / outgroups ({len(non_selected_samples)}): {non_selected_samples}")

    if not non_selected_samples:
        logging.error("No non-selected samples found in TSV to compare against. Exiting.")
        sys.exit(1)

    measures = ["manhattan", "patristic", "jaccard"]
    # Filter measures that exist in df columns
    available_measures = [m for m in measures if any(c.endswith(f"_{m}") for c in df.columns)]
    if not available_measures:
        logging.error("No distance measures found in TSV columns. Exiting.")
        sys.exit(1)

    summary_df = calculate_nearest_outgroup_summary(df, selected_samples, non_selected_samples, available_measures)

    # Write summary TSV
    os.makedirs(os.path.dirname(os.path.abspath(args.output_tsv)), exist_ok=True)
    summary_df.to_csv(args.output_tsv, sep='\t', index=False)
    logging.info(f"Tabulated summary written to {args.output_tsv}")

    # Generate plot if requested
    if args.output_plot:
        plot_measure = args.measure
        if plot_measure not in available_measures:
            logging.warning(f"Requested plot measure '{plot_measure}' not found in TSV. Using '{available_measures[0]}'.")
            plot_measure = available_measures[0]
        plot_grouped_stacked_bar(summary_df, plot_measure, selected_samples, non_selected_samples, args.output_plot)


if __name__ == "__main__":
    main()
