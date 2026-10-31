import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def plot_retention_by_chromosome(tsv_path, output_filename="retention_plot.png"):
    # 1. Load and clean the dataset
    clean_lines = []
    with open(tsv_path, 'r') as f:
        for line in f:
            # Keep the header row but strip the leading '#'
            if line.startswith('#ref_chrom'):
                clean_lines.append(line.lstrip('#').strip())
            # Skip all other comment lines
            elif not line.startswith('#'):
                clean_lines.append(line.strip())

    import io
    df = pd.read_csv(io.StringIO('\n'.join(clean_lines)), sep='\t')

    # 2. Filter out unassigned sequences (e.g., _NONE_ or _AMBIGUOUS_)
    df = df[~df['ref_chrom'].isin(['_NONE_', '_AMBIGUOUS_'])]

    # 3. Calculate retention percentages relative to the input base pairs
    # input_bp is 100%, and each subsequent stage is a fraction of that input
    df['Input'] = 100.0
    df['Full Graph'] = (df['full_bp'] / df['input_bp']) * 100
    df['Clipped'] = (df['clip_bp'] / df['input_bp']) * 100
    df['Filtered'] = (df['filter_bp'] / df['input_bp']) * 100

    # 4. Melt the dataframe to long format for Seaborn plotting
    melted_df = df.melt(
        id_vars=['ref_chrom', 'genome'],
        value_vars=['Input', 'Full Graph', 'Clipped', 'Filtered'],
        var_name='Pipeline_Stage',
        value_name='Retained_Percent'
    )

    # Ensure categorical order of stages is preserved in the plot
    melted_df['Pipeline_Stage'] = pd.Categorical(
        melted_df['Pipeline_Stage'],
        categories=['Input', 'Full Graph', 'Clipped', 'Filtered'],
        ordered=True
    )

    # 5. Create a FacetGrid to plot each chromosome in a separate panel
    sns.set_theme(style="whitegrid", font_scale=1.1)

    # Using catplot to create a point plot (which beautifully shows the drop-off curve)
    g = sns.catplot(
        data=melted_df,
        x='Pipeline_Stage',
        y='Retained_Percent',
        hue='genome',
        col='ref_chrom',
        col_wrap=4,  # Adjust depending on how many chromosomes you have
        kind='point',  # 'point' draws lines between stages; use 'bar' for bar plots
        markers='o',
        scale=0.8,
        height=4,
        aspect=1.2,
        palette='tab20'  # Good for many distinct genomes
    )

    # 6. Formatting and labels
    g.set_axis_labels("Graph Pruning Stage", "Sequence Retained (%)")
    g.set_titles(col_template="{col_name}")
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle("Sequence Retention Across Pangenome Pruning Stages by Chromosome", fontsize=16, fontweight='bold')

    # Rotate x-axis labels for better readability
    for ax in g.axes.flat:
        ax.tick_params(axis='x', rotation=45)
        # Optional: Set y-limit if you only want to focus on the top percentages (e.g., 80-100%)
        # ax.set_ylim(80, 105)

    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.show()


# Run the function on your file
plot_retention_by_chromosome(sys.argv[1])