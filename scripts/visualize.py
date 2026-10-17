import os
import logging
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm


def parse_valid_pairs(columns, measure_suffix, target_samples=None, focus_samples=None):
    """
    Extracts valid pair strings from DataFrame columns ending with `_{measure_suffix}`.
    target_samples (-s): BOTH samples must be in this list.
    focus_samples (-S): AT LEAST ONE sample must be in this list.
    """
    suffix = f"_{measure_suffix}"
    all_pair_strs = [c[:-len(suffix)] for c in columns if c.endswith(suffix)]
    valid_pairs = []

    target_set = set(target_samples) if target_samples else None
    focus_set = set(focus_samples) if focus_samples else None

    if not target_set and not focus_set:
        return all_pair_strs

    for pair_str in all_pair_strs:
        parts = pair_str.split('-')
        is_valid_pair = False

        for i in range(1, len(parts)):
            s1 = "-".join(parts[:i])
            s2 = "-".join(parts[i:])

            in_pool = True
            if target_set:
                if s1 not in target_set or s2 not in target_set:
                    in_pool = False

            in_focus = True
            if focus_set:
                if s1 not in focus_set and s2 not in focus_set:
                    in_focus = False

            if in_pool and in_focus:
                is_valid_pair = True
                break

        if is_valid_pair:
            valid_pairs.append(pair_str)

    return valid_pairs


def parse_samples_from_columns(columns):
    """
    Find sample names in TSV columns matching `{sample}_called_sites`.
    """
    samples = []
    suffix = "_called_sites"
    for col in columns:
        if col.endswith(suffix):
            sample_name = col[:-len(suffix)]
            samples.append(sample_name)
    return samples


def filter_samples(all_samples, target_samples=None, focus_samples=None):
    """
    Filter individual samples based on target (-s) or focus (-S) sample lists.
    """
    target_set = set(target_samples) if target_samples else None
    focus_set = set(focus_samples) if focus_samples else None

    if not target_set and not focus_set:
        return all_samples

    filtered = []
    for s in all_samples:
        if target_set and s not in target_set:
            continue
        if focus_set and s not in focus_set:
            continue
        filtered.append(s)
    return filtered


def build_color_map(keys):
    """
    Build a consistent color mapping for a list of keys (e.g. sample pairs or sample names).
    """
    prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    color_cycle = itertools.cycle(prop_cycle)
    color_map = {}
    for key in sorted(keys):
        color_map[key] = next(color_cycle)
    return color_map


def plot_distance(chrom, measure, measure_label, density_type, density_label, outformat, frac, outdir, valid_pairs, color_map, df_chrom, x_all, density):
    fig, ax1 = plt.subplots(figsize=(14, 6))

    # --- 1. Background Density (Right Y-Axis) ---
    ax2 = ax1.twinx()
    ax2.fill_between(x_all, 0, density, color='gray', alpha=0.2, step='mid', label=density_label)
    ax2.set_ylabel(density_label, color='dimgray', fontsize=11, fontweight='bold')
    ax2.tick_params(axis='y', labelcolor='dimgray')
    ax2.set_ylim(bottom=0)

    ax1.set_zorder(ax2.get_zorder() + 1)
    ax1.patch.set_visible(False)

    # --- 2. Pairwise Distance Measures (Left Y-Axis) ---
    lowess = sm.nonparametric.lowess

    for pair in valid_pairs:
        col_name = f"{pair}_{measure}"
        if col_name not in df_chrom.columns:
            continue

        x = df_chrom['Position_Mb'].values
        y = pd.to_numeric(df_chrom[col_name], errors='coerce').values

        mask = ~np.isnan(x) & ~np.isnan(y)
        x_valid = x[mask]
        y_valid = y[mask]

        if len(x_valid) == 0:
            continue

        color = color_map.get(pair, 'blue')

        # Scatter
        ax1.scatter(x_valid, y_valid, color=color, alpha=0.2, s=15, edgecolors='none')

        # LOESS curve
        if len(x_valid) > 5:
            try:
                smoothed = lowess(y_valid, x_valid, frac=frac)
                ax1.plot(smoothed[:, 0], smoothed[:, 1], color=color, linewidth=2.5, label=pair)
            except Exception:
                ax1.plot(x_valid, y_valid, color=color, linewidth=2.5, label=pair)
        else:
            ax1.plot(x_valid, y_valid, color=color, linewidth=2.5, label=pair)

    ax1.set_title(f"Chromosome: {chrom} | {measure_label} ({density_label})")
    ax1.set_xlabel("Genomic Position (Mb)")
    ax1.set_ylabel("Distance", fontsize=11, fontweight='bold')
    ax1.grid(True, linestyle='--', alpha=0.4)

    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    from matplotlib.patches import Patch
    lines1.append(Patch(facecolor='gray', alpha=0.2))
    labels1.append(density_label)

    ax1.legend(lines1, labels1, bbox_to_anchor=(1.08, 1), loc='upper left', borderaxespad=0., fontsize='small')

    plt.tight_layout()

    out_filename = f"{chrom}_{measure}_density_{density_type}.{outformat}"
    out_file = os.path.join(outdir, out_filename)
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()


def plot_proportion(chrom, prop_type, prop_label, y_col_suffix, outformat, frac, outdir, samples, color_map, df_chrom, x_all, density):
    fig, ax1 = plt.subplots(figsize=(14, 6))

    # --- 1. Background Density (Total_variants) ---
    ax2 = ax1.twinx()
    ax2.fill_between(x_all, 0, density, color='gray', alpha=0.2, step='mid', label='Total Variants')
    ax2.set_ylabel("Total Variants per Window", color='dimgray', fontsize=11, fontweight='bold')
    ax2.tick_params(axis='y', labelcolor='dimgray')
    ax2.set_ylim(bottom=0)

    ax1.set_zorder(ax2.get_zorder() + 1)
    ax1.patch.set_visible(False)

    # --- 2. Sample Proportions (Left Y-Axis) ---
    lowess = sm.nonparametric.lowess
    total_variants = pd.to_numeric(df_chrom['Total_variants'], errors='coerce').replace(0, np.nan)

    for sample in samples:
        col_name = f"{sample}_{y_col_suffix}"
        if col_name not in df_chrom.columns:
            continue

        raw_vals = pd.to_numeric(df_chrom[col_name], errors='coerce')
        y = (raw_vals / total_variants).values
        x = df_chrom['Position_Mb'].values

        mask = ~np.isnan(x) & ~np.isnan(y)
        x_valid = x[mask]
        y_valid = y[mask]

        if len(x_valid) == 0:
            continue

        color = color_map.get(sample, 'blue')

        ax1.scatter(x_valid, y_valid, color=color, alpha=0.2, s=15, edgecolors='none')

        if len(x_valid) > 5:
            try:
                smoothed = lowess(y_valid, x_valid, frac=frac)
                ax1.plot(smoothed[:, 0], smoothed[:, 1], color=color, linewidth=2.5, label=sample)
            except Exception:
                ax1.plot(x_valid, y_valid, color=color, linewidth=2.5, label=sample)
        else:
            ax1.plot(x_valid, y_valid, color=color, linewidth=2.5, label=sample)

    ax1.set_title(f"Chromosome: {chrom} | {prop_label}")
    ax1.set_xlabel("Genomic Position (Mb)")
    ax1.set_ylabel("Proportion", fontsize=11, fontweight='bold')
    ax1.set_ylim(0, 1.05)
    ax1.grid(True, linestyle='--', alpha=0.4)

    lines1, labels1 = ax1.get_legend_handles_labels()
    from matplotlib.patches import Patch
    lines1.append(Patch(facecolor='gray', alpha=0.2))
    labels1.append('Total Variants')

    ax1.legend(lines1, labels1, bbox_to_anchor=(1.08, 1), loc='upper left', borderaxespad=0., fontsize='small')

    plt.tight_layout()

    out_filename = f"{chrom}_{prop_type}.{outformat}"
    out_file = os.path.join(outdir, out_filename)
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()


def run(args):
    logging.info(f"Loading data from {args.input_tsv}...")
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.input_tsv, sep='\t', na_values=["NA", "TREE_ERROR"])

    # Position midpoint in Mb
    df['Position_Mb'] = ((df['Window_start'] + df['Window_end']) / 2) / 1_000_000
    df = df.sort_values(by=['Chromosome', 'Position_Mb'])

    target_samples = args.samples.split(',') if args.samples else None
    focus_samples = args.focus.split(',') if args.focus else None

    # Parse sample pairs
    all_pairs = parse_valid_pairs(df.columns, 'patristic', target_samples, focus_samples)
    pair_color_map = build_color_map(all_pairs)

    # Parse individual samples
    all_samples = parse_samples_from_columns(df.columns)
    active_samples = filter_samples(all_samples, target_samples, focus_samples)
    sample_color_map = build_color_map(active_samples)

    chromosomes = df['Chromosome'].unique()

    measures = ['patristic', 'jaccard', 'manhattan']
    measure_labels = {
        'patristic': 'Patristic Distance (NJ Tree)',
        'jaccard': 'Jaccard Distance',
        'manhattan': 'Manhattan Distance'
    }

    for chrom in chromosomes:
        df_chrom = df[df['Chromosome'] == chrom]
        logging.info(f"Generating plots for chromosome: {chrom}...")

        x_all = df_chrom['Position_Mb'].values
        var_density = pd.to_numeric(df_chrom['Total_variants'], errors='coerce').fillna(0).values
        bp_density = pd.to_numeric(df_chrom['Total_bp'], errors='coerce').fillna(0).values

        # 1. Distance Plots (in duplicate: variants background & bp background)
        for measure in measures:
            # Duplicate #1: Total_variants background
            plot_distance(
                chrom, measure, measure_labels[measure],
                'variants', 'Total Variants per Window',
                args.format, args.frac, args.outdir,
                all_pairs, pair_color_map, df_chrom, x_all, var_density
            )
            # Duplicate #2: Total_bp background
            plot_distance(
                chrom, measure, measure_labels[measure],
                'bp', 'Total Basepairs per Window',
                args.format, args.frac, args.outdir,
                all_pairs, pair_color_map, df_chrom, x_all, bp_density
            )

        # 2. Variant Proportion Plots
        # Variant occupancy: called_sites / Total_variants
        plot_proportion(
            chrom, 'variant_occupancy', 'Variant Occupancy per Sample',
            'called_sites', args.format, args.frac, args.outdir,
            active_samples, sample_color_map, df_chrom, x_all, var_density
        )
        # Unique variant proportion: unique_alleles / Total_variants
        plot_proportion(
            chrom, 'unique_variants_proportion', 'Unique Variants Proportion per Sample',
            'unique_alleles', args.format, args.frac, args.outdir,
            active_samples, sample_color_map, df_chrom, x_all, var_density
        )

    logging.info(f"All plots saved to {args.outdir}/")
