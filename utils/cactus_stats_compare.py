import sys
import io
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def load_graph_stats(filepath, label):
    """
    Loads a graph-stats.tsv file, strips comments/headers,
    and adds a graph label identifier.
    """
    clean_lines = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#ref_chrom') or line.startswith('# ref_chrom'):
                clean_lines.append(line.lstrip('#').strip())
            elif not line.startswith('#'):
                clean_lines.append(line.strip())

    df = pd.read_csv(io.StringIO('\n'.join(clean_lines)), sep=r'\s+')
    df['graph'] = label
    # Convert sequence length to Megabases for readability
    df['length_mb'] = df['length'] / 1e6
    # Compute graph complexity ratio (Edges / Nodes)
    df['edge_to_node_ratio'] = df['edges'] / df['nodes']
    return df


def plot_graph_comparison(file1, label1, file2, label2, output_prefix="pangenome_graph_comparison"):
    # Load both datasets
    df1 = load_graph_stats(file1, label1)
    df2 = load_graph_stats(file2, label2)

    # Combined long-format dataframe for grouped seaborn plotting
    df_combined = pd.concat([df1, df2], ignore_index=True)

    # Merge for wide-format comparative delta calculations
    merged = pd.merge(df1, df2, on='ref_chrom', suffixes=(f'_{label1}', f'_{label2}'))

    # Calculate percentage differences relative to Graph 2 (or base reference)
    merged['pct_diff_nodes'] = ((merged[f'nodes_{label1}'] - merged[f'nodes_{label2}']) / merged[
        f'nodes_{label2}']) * 100
    merged['pct_diff_edges'] = ((merged[f'edges_{label1}'] - merged[f'edges_{label2}']) / merged[
        f'edges_{label2}']) * 100
    merged['pct_diff_length'] = ((merged[f'length_{label1}'] - merged[f'length_{label2}']) / merged[
        f'length_{label2}']) * 100

    # Set up global plot aesthetics
    sns.set_theme(style="whitegrid", font_scale=1.0)
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1])

    # ---------------------------------------------------------
    # Panel A: Side-by-Side Node & Edge Counts
    # ---------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, :])

    # Melt combined dataframe for nodes and edges
    melted_counts = pd.melt(
        df_combined,
        id_vars=['ref_chrom', 'graph'],
        value_vars=['nodes', 'edges'],
        var_name='metric', value_name='count'
    )
    melted_counts['count_millions'] = melted_counts['count'] / 1e6

    sns.barplot(
        data=melted_counts,
        x='ref_chrom', y='count_millions',
        hue='graph', hue_order=[label1, label2],
        ax=ax_a, palette='Set2', edgecolor='black', linewidth=0.6
    )

    ax_a.set_title("A. Direct Comparison of Node and Edge Counts per Chromosome", fontsize=14, fontweight='bold',
                   pad=10)
    ax_a.set_xlabel("Reference Chromosome", fontsize=11)
    ax_a.set_ylabel("Count (Millions)", fontsize=11)
    ax_a.tick_params(axis='x', rotation=45)
    ax_a.legend(title="Pangenome Graph", loc="upper right")

    # ---------------------------------------------------------
    # Panel B: Percentage Change / Delta Metrics (% Difference)
    # ---------------------------------------------------------
    ax_b = fig.add_subplot(gs[1, 0])

    delta_melted = pd.melt(
        merged,
        id_vars=['ref_chrom'],
        value_vars=['pct_diff_nodes', 'pct_diff_edges', 'pct_diff_length'],
        var_name='Metric', value_name='pct_change'
    )

    # Format labels for legend
    metric_names = {
        'pct_diff_nodes': 'Nodes (%)',
        'pct_diff_edges': 'Edges (%)',
        'pct_diff_length': 'Total Length (%)'
    }
    delta_melted['Metric'] = delta_melted['Metric'].map(metric_names)

    sns.barplot(
        data=delta_melted,
        x='ref_chrom', y='pct_change', hue='Metric',
        ax=ax_b, palette='Blues_r', edgecolor='black', linewidth=0.5
    )

    ax_b.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax_b.set_title(f"B. Structural Delta ({label1} vs. {label2})", fontsize=13, fontweight='bold', pad=10)
    ax_b.set_xlabel("Reference Chromosome", fontsize=11)
    ax_b.set_ylabel(f"Percentage Gain / Loss (%)", fontsize=11)
    ax_b.tick_params(axis='x', rotation=45)
    ax_b.legend(title="Metric", loc="upper right")

    # ---------------------------------------------------------
    # Panel C: Graph Complexity Ratio (Edges / Nodes)
    # ---------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 1])

    sns.pointplot(
        data=df_combined,
        x='ref_chrom', y='edge_to_node_ratio', hue='graph',
        hue_order=[label1, label2], ax=ax_c, palette='Set2',
        markers=['o', 's'], linestyles=['-', '--'], scale=0.9
    )

    ax_c.set_title("C. Graph Reticulation Complexity (Edge-to-Node Ratio)", fontsize=13, fontweight='bold', pad=10)
    ax_c.set_xlabel("Reference Chromosome", fontsize=11)
    ax_c.set_ylabel("Complexity Ratio (E / V)", fontsize=11)
    ax_c.tick_params(axis='x', rotation=45)
    ax_c.legend(title="Pangenome Graph", loc="upper right")

    plt.tight_layout()
    plt.savefig(f"{output_prefix}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}.pdf", bbox_inches='tight')
    plt.close()

    print(f"Comparison plot successfully saved as {output_prefix}.png and {output_prefix}.pdf")


if __name__ == "__main__":
    # Example command-line execution:
    # python compare_graphs.py all_genomes_stats.tsv hq_subset_stats.tsv
    if len(sys.argv) >= 3:
        f1 = sys.argv[1]
        f2 = sys.argv[2]
        l1 = sys.argv[3] if len(sys.argv) > 3 else "All_Genomes"
        l2 = sys.argv[4] if len(sys.argv) > 4 else "HQ_Subset"
    else:
        # Default filenames
        f1 = "graph1_stats.tsv"
        f2 = "graph2_stats.tsv"
        l1 = "All_Genomes"
        l2 = "HQ_Subset"

    plot_graph_comparison(f1, l1, f2, l2)