import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import traceback

"""
Script: Total Gene Coverage in Synteny Blocks

Purpose: Calculates and visualizes what percentage of total annotated genes are found in synteny blocks 
for each species pair.

Input:
- TSV files with genes in synteny blocks (final_output/*_synteny_tables_genes_fix.tsv)
- BED6 files with genome annotations (*_genes_only_stableid_fix.bed6)

Output:
- Pie charts showing gene coverage per species per pair (*_total_number_genes_piechart_*_fix.png)
- Summary TSV with all percentage values (summary_points_fix.tsv)
- Summary boxplot across all pairs (summary_boxplot_fix.png)

Prints:
- Error messages for failed pairs
- Number of collected data points
- Summary statistics for boxplot
"""

# Plot piechart of ratio of total number of genes and genes located in synteny blocks
# and summarize all per-direction percentages in a single boxplot.

# Species pairs for input
pairs = [
    ("can", "mac"),
    ("mac", "can"),
    ("can", "rat"),
    ("rat", "can"),
    ("can", "mus"),
    ("mus", "can"),
    ("mac", "rat"),
    ("rat", "mac"),
    ("mac", "mus"),
    ("mus", "mac"),
    ("rat", "mus"),
    ("mus", "rat"),
]

# Names of species in annotation files
annotation_names = [
    "Canis",
    "Macaca",
    "Rat",
    "Mus"
]

# Dictionary of annotation name to species abbreviation 
annotation_map = {
    "can": "Canis",
    "mac": "Macaca",
    "rat": "Rat",
    "mus": "Mus"
}

# ---- summary boxplot config ----
OUT_DIR = "results/data_check/total_number_genes"
os.makedirs(OUT_DIR, exist_ok=True)
METRIC_NAME = "genes_in_blocks_pct"
all_metrics = []   # rows: metric, pair_dir, species, value (percentage)


def count_genes_int(gene_table_path, source, target):
    """
    Count unique genes located in synteny blocks for both species.
    
    Input:
    - gene_table_path: str, path to TSV file with columns Region_ID, Species, Gene-ID
    - source: str, source species abbreviation (e.g., "can")
    - target: str, target species abbreviation (e.g., "mac")
    
    Output:
    - Tuple of (n_source: int, n_target: int), number of unique genes in blocks for each species
    
    Prints: Nothing
    """
    df = pd.read_csv(gene_table_path, sep="\t", dtype=str)
    df = df.dropna(subset=["Region_ID", "Species", "Gene-ID"])
    df = df.drop_duplicates(subset=["Region_ID", "Species", "Gene-ID"])

    src_name = annotation_map[source]
    tgt_name = annotation_map[target]

    # partial match (case-insensitive)
    n_source = df[df["Species"].str.contains(src_name, case=False, na=False)]["Gene-ID"].nunique()
    n_target = df[df["Species"].str.contains(tgt_name, case=False, na=False)]["Gene-ID"].nunique()

    return n_source, n_target


def count_genes_annotation(gene_annotation_path):
    """
    Count total number of unique genes in genome annotation.
    
    Input:
    - gene_annotation_path: str, path to BED6 file (gene-only annotation)
    
    Output:
    - int, number of unique genes
    
    Prints: Nothing
    """
    colnames = ["chrom", "start", "end", "GeneID", "score", "strand"]
    gene_annotation = pd.read_csv(gene_annotation_path, sep="\t", names=colnames, dtype=str)
    return gene_annotation["GeneID"].nunique()


def plot_total_number_genes(number_genes_ann, genes_in_blocks, species_role, source, target, output_plot):
    """
    Create pie chart showing proportion of genes in synteny blocks vs not in blocks.
    
    Input:
    - number_genes_ann: int, total number of genes from annotation
    - genes_in_blocks: int, number of genes found in synteny blocks
    - species_role: str, either "query" or "target"
    - source: str, source species abbreviation
    - target: str, target species abbreviation
    - output_plot: str, output path for pie chart (PNG)
    
    Output:
    - PNG pie chart
    
    Prints: Nothing
    """
    genes_not_in_blocks = number_genes_ann - genes_in_blocks
    labels = ['In synteny blocks', 'Not in synteny blocks']
    sizes = [genes_in_blocks, genes_not_in_blocks]
    colors = ['#66b3ff', '#ff9999']
    ratio = genes_in_blocks / number_genes_ann * 100 if number_genes_ann > 0 else 0.0

    plt.figure(figsize=(5,5), dpi=150)
    plt.pie(
        sizes,
        labels=labels,
        autopct='%1.1f%%',
        startangle=90,
        colors=colors,
        wedgeprops={'edgecolor': 'black'}
    )
    plt.title(f"{source}->{target} ({species_role}): {ratio:.1f}% of genes in synteny blocks")
    plt.savefig(output_plot)
    plt.close()


def finalize_boxplot(all_metrics, metric_name, out_csv, out_png):
    """
    Create summary boxplot showing distribution of gene coverage across all pairs.
    
    Input:
    - all_metrics: list of dict, each with keys: metric, pair_dir, species, value
    - metric_name: str, name of the metric to plot
    - out_csv: str, output path for TSV with all data points
    - out_png: str, output path for boxplot (PNG)
    
    Output:
    - TSV file with all percentage values
    - PNG boxplot with jittered points colored by direction
    
    Prints:
    - Warning if no data for metric
    - Number of data points and unique pairs used in boxplot
    """
    df = pd.DataFrame(all_metrics)
    df_m = df[df["metric"] == metric_name].copy()
    if df_m.empty:
        print(f"No rows for metric={metric_name}. Nothing to plot.")
        return

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df_m.to_csv(out_csv, sep="\t", index=False)

    n = len(df_m)
    n_pairs = df_m["pair_dir"].nunique()
    print(f"Boxplot will use {n} bars from {n_pairs} query->target pair-dirs.")

    plt.figure(figsize=(9,7), dpi=150)
    plt.boxplot(df_m["value"], vert=True, showmeans=True, widths=0.35, tick_labels=[metric_name])

    # jittered points with direction-based colors (A->B vs B->A by lexicographic order)
    rng = np.random.default_rng(42)
    x = rng.normal(1, 0.025, size=n)
    colors = df_m["pair_dir"].apply(
        lambda p: "tab:blue" if "->" in p and p.split("->")[0] < p.split("->")[1] else "tab:orange"
    )
    plt.scatter(x, df_m["value"], c=colors, s=35, alpha=0.9,
                edgecolors="black", linewidths=0.4, zorder=3)
    plt.scatter([], [], c="tab:blue", label="A->B")
    plt.scatter([], [], c="tab:orange", label="B->A")
    plt.legend(frameon=False, loc="lower left", fontsize=9)

    plt.ylim(0, 100)
    plt.ylabel("% of genes in synteny blocks", fontsize=12)
    plt.title(f"{metric_name}: summary across all query->target bars (n={n})", fontsize=13)
    plt.grid(axis="y", linestyle="--", alpha=0.4, zorder=0)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def process_pairs(source, target):
    """
    Process one species pair: count genes, create pie charts, collect summary data.
    
    Input:
    - source: str, source species abbreviation (e.g., "can")
    - target: str, target species abbreviation (e.g., "mac")
    
    Output:
    - Two PNG pie charts (one per species)
    - Appends percentage values to global all_metrics list
    
    Prints: Nothing (errors handled by main())
    """

    # --- inputs ---
    genes_tsv     = f"results/genes_in_blocks/mammalia/region_gene_tables/{source}_{target}_synteny_tables_genes_fix.tsv"
    annotation_src = f"results/ann_genes/{annotation_map[source]}_genes_only_stableid_fix.bed6"
    annotation_tgt = f"results/ann_genes/{annotation_map[target]}_genes_only_stableid_fix.bed6"

    # --- outputs ---
    piechart_src = f"{OUT_DIR}/{source}_{target}_total_number_genes_piechart_{source}_fix.png"
    piechart_tgt = f"{OUT_DIR}/{source}_{target}_total_number_genes_piechart_{target}_fix.png"

    # 1) Count total genes per annotation
    total_src = count_genes_annotation(annotation_src)
    total_tgt = count_genes_annotation(annotation_tgt)

    # 2) Count genes in synteny blocks
    in_blocks_src, in_blocks_tgt = count_genes_int(genes_tsv, source, target)

    # 3) Plot both species' coverage
    plot_total_number_genes(total_src, in_blocks_src, "query",  source, target, piechart_src)
    plot_total_number_genes(total_tgt, in_blocks_tgt, "target", source, target, piechart_tgt)

    # 4) Collect for summary boxplot (two points per direction)
    pct_src = (in_blocks_src / total_src * 100.0) if total_src else 0.0
    pct_tgt = (in_blocks_tgt / total_tgt * 100.0) if total_tgt else 0.0
    all_metrics.append({"metric": METRIC_NAME, "pair_dir": f"{source}->{target}", "species": source, "value": pct_src})
    all_metrics.append({"metric": METRIC_NAME, "pair_dir": f"{source}->{target}", "species": target, "value": pct_tgt})


def main():
    """
    Main function: process all species pairs and create summary visualizations.
    
    Prints:
    - Error messages for any failed pairs
    - Total number of collected data points
    - Summary statistics from finalize_boxplot()
    """
    for src, tgt in pairs:
        try:
            process_pairs(src, tgt)
        except Exception as e:
            print(f"[ERROR] {src}-{tgt}: {e}")
            traceback.print_exc()

    # Summary boxplot + CSV of points
    print(f"\nCollected {len(all_metrics)} values (expected {2*len(pairs)}).")
    finalize_boxplot(
        all_metrics,
        metric_name=METRIC_NAME,
        out_csv=os.path.join(OUT_DIR, "summary_points_fix.tsv"),
        out_png=os.path.join(OUT_DIR, "summary_boxplot_fix.png"),
    )

if __name__ == "__main__":
    main()