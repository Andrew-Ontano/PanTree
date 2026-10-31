import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_binned_coverage(bed_file, window_size=100, output_file="coverage_plot.png", primary_prefix=""):
    # 1. Load the mosdepth regions file
    # Mosdepth outputs bed.gz without headers, so we assign them
    df = pd.read_csv(bed_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'cov'])
    # Optional: Filter for primary reference chromosomes to exclude unplaced scaffolds
    # Adjust this list based on your specific Apis mellifera reference accessions
    primary_chroms = [f"NC_037{i}.1_RagTag" for i in range(638, 654)]
    df = df[df['chrom'].isin(primary_chroms)]

    # 2. Convert start positions to megabases for cleaner X-axis labels
    df['Position_MB'] = df['start'] / 1_000_000

    # 3. Apply a rolling mean to smooth the 500bp noise
    # window_size=100 on 500bp bins = 50kb smoothing window
    df['Smoothed_Cov'] = df.groupby('chrom')['cov'].transform(
        lambda x: x.rolling(window=window_size, min_periods=1, center=True).mean()
    )

    # 4. Plot using a Seaborn FacetGrid
    sns.set_theme(style="whitegrid")

    g = sns.relplot(
        data=df,
        x='Position_MB',
        y='Smoothed_Cov',
        col='chrom',
        col_wrap=4,
        kind='line',
        height=3,
        aspect=1.5,
        color='teal',
        linewidth=1
    )

    # 5. Formatting
    g.set_axis_labels("Position (Mb)", "Depth of Coverage (X)")
    g.set_titles(col_template="{col_name}")
    g.fig.suptitle("Smoothed Read Depth Coverage per Chromosome", y=1.02, fontsize=16)

    # Add a horizontal line representing the expected haploid coverage
    median_cov = df['cov'].median()
    for ax in g.axes.flat:
        ax.axhline(median_cov, color='red', linestyle='--', alpha=0.7, label=f"Median ({median_cov:.1f}X)")
        ax.set_ylim(0, df['Smoothed_Cov'].quantile(0.99) * 1.2)  # Cap Y-axis to ignore extreme peaks

    plt.savefig(output_file, dpi=300, bbox_inches='tight')

# Execution:
plot_binned_coverage(sys.argv[1], 100, output_file=sys.argv[3])
