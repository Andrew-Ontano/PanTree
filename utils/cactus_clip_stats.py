import sys
import io
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def load_cactus_stats(filepath):
    """
    Loads cactus-pangenome stats TSV, filtering metadata comments
    while preserving the #genome header line.
    """
    clean_lines = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#genome') or line.startswith('# genome'):
                clean_lines.append(line.lstrip('#').strip())
            elif not line.startswith('#'):
                clean_lines.append(line.strip())

    df = pd.read_csv(io.StringIO('\n'.join(clean_lines)), sep=r'\s+')
    return df


def plot_cactus_clipping_stats(df, output_prefix="pangenome_clipping"):
    # Set global aesthetic style
    sns.set_theme(style="whitegrid", font_scale=1.0)

    # ---------------------------------------------------------
    # 1. Stacked Bar Chart: Sequence Loss by Stage (%)
    # ---------------------------------------------------------
    df_stages = df[df['reason'] != 'TOTAL'].copy()

    # Pivot dataframe for stacked plotting
    stage_pivot = df_stages.pivot(index='genome', columns='reason', values='pct_of_input').fillna(0)

    # Order stages logically (chromosome assignment -> graph alignment -> clipping -> filtering)
    stage_order = ['ambiguous', 'unassigned', 'no_chromosome_graph', 'unaligned', 'clip', 'filter']
    present_stages = [s for s in stage_order if s in stage_pivot.columns]
    stage_pivot = stage_pivot[present_stages]

    # Calculate total loss for sorting genomes in plot
    stage_pivot['total_pct'] = stage_pivot.sum(axis=1)
    stage_pivot = stage_pivot.sort_values(by='total_pct', ascending=True)
    stage_pivot = stage_pivot.drop(columns=['total_pct'])

    # Plot horizontal stacked bar chart
    fig, ax = plt.subplots(figsize=(10, max(5, len(stage_pivot) * 0.4)))
    colors = sns.color_palette("Set2", len(present_stages))

    stage_pivot.plot(kind='barh', stacked=True, ax=ax, color=colors, edgecolor='black', linewidth=0.5)

    ax.set_title("Sequence Absent from Pangenome Graph by Pipeline Stage", fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel("Percentage of Input Sequence (%)", fontsize=12)
    ax.set_ylabel("Genome / Subspecies", fontsize=12)
    ax.legend(title="Loss Reason / Stage", bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_by_stage.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_by_stage.pdf", bbox_inches='tight')
    plt.close()

    # ---------------------------------------------------------
    # 2. Total Sequence Loss Ranking (% and Mb)
    # ---------------------------------------------------------
    df_total = df[df['reason'] == 'TOTAL'].copy()
    df_total['mb_lost'] = df_total['bp'] / 1e6
    df_total = df_total.sort_values(by='pct_of_input', ascending=False)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, max(5, len(df_total) * 0.4)), sharey=True)

    # Plot % of input sequence lost
    sns.barplot(data=df_total, y='genome', x='pct_of_input', ax=ax1, palette="Reds_r", edgecolor='black')
    ax1.set_title("Total Input Sequence Excluded (%)", fontsize=12, fontweight='bold')
    ax1.set_xlabel("Percentage of Input Genome (%)", fontsize=11)
    ax1.set_ylabel("Genome / Subspecies", fontsize=11)

    # Annotate bar percentages
    for p in ax1.patches:
        width = p.get_width()
        ax1.annotate(f"{width:.2f}%", (width + 0.3, p.get_y() + p.get_height() / 2.),
                     ha='left', va='center', fontsize=9)

    # Plot total Megabases lost
    sns.barplot(data=df_total, y='genome', x='mb_lost', ax=ax2, palette="Blues_r", edgecolor='black')
    ax2.set_title("Total Base Pairs Excluded (Mb)", fontsize=12, fontweight='bold')
    ax2.set_xlabel("Excluded Base Pairs (Mb)", fontsize=11)
    ax2.set_ylabel("")

    # Annotate bar Megabases
    for p in ax2.patches:
        width = p.get_width()
        ax2.annotate(f"{width:.1f} Mb", (width + (max(df_total['mb_lost']) * 0.01), p.get_y() + p.get_height() / 2.),
                     ha='left', va='center', fontsize=9)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_total_summary.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_total_summary.pdf", bbox_inches='tight')
    plt.close()

    print(f"Plots saved successfully as {output_prefix}_by_stage.png and {output_prefix}_total_summary.png")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        tsv_file = sys.argv[1]
    else:
        tsv_file = "stats.tsv"  # Replace with path to your cactus output TSV

    df_stats = load_cactus_stats(tsv_file)
    plot_cactus_clipping_stats(df_stats)