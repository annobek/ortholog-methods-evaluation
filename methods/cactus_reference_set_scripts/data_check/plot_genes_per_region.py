import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

"""
Script: Gene Distribution Analysis in Synteny Blocks

Purpose: Analyzes gene distribution across synteny blocks for mammalian species pairs.

Input:
- TSV files with gene assignments to synteny blocks (region_gene_tables/*_synteny_tables_genes_fix.tsv)
- TSV files with synteny block coordinates (region_gene_tables/*_synteny_tables_regions_fix.tsv)

Output:
- TSV files with gene counts per block (*_number_genes_in_blocks.tsv)
- PNG plots: length vs gene number (query and target), gene distribution histograms (overview and 0-100 detail)

Prints:
- Distribution statistics (mean, median, std dev, quartiles, min/max genes per block)
- Percentage of blocks with 0 genes
- Median gene density per Mb
"""

# Analyses number of genes per block:
# - count_genes() - Counts number of genes per synteny block, saves the table as TSV
# - check_results() - Calculates statistical data
# - Plots:
#   - plot_gene_distribution() - histogram of number of genes per block


# input: abbreviations for query, target and the directory where they are located
pairs = [
    ("can", "mac", "mammalia"),
    ("mac", "can", "mammalia"),
    ("can", "rat", "mammalia"),
    ("rat", "can", "mammalia"),
    ("can", "mus", "mammalia"),
    ("mus", "can", "mammalia"),
    ("mac", "rat", "mammalia"),
    ("rat", "mac", "mammalia"),
    ("mac", "mus", "mammalia"),
    ("mus", "mac", "mammalia"),
    ("rat", "mus", "mammalia"),
    ("mus", "rat", "mammalia"),
]

# Full names of the species

species_full = {
    "can": "Canis lupus familiaris",
    "mac": "Macaca fascicularis",
    "rat": "Rattus norvegicus",
    "mus": "Mus musculus",
}


def count_genes(gene_table, query, target):
    """
    Extract number of genes per synteny block as DataFrame.
    
    Input: 
    - gene_table: str, path to TSV file with columns Region_ID, Species, Gene-ID
    - query: str, full species name for query species
    - target: str, full species name for target species
    
    Output:
    - DataFrame with columns: Region_ID (str), n_genes_query (int), n_genes_target (int)
    
    Prints: Nothing
    """
    # Read and de-dup on (Region_ID, Species, Gene-ID)
    df = pd.read_csv(gene_table, sep="\t", dtype=str)
    df = df.dropna(subset=["Region_ID", "Species", "Gene-ID"])
    df = df.drop_duplicates(subset=["Region_ID", "Species", "Gene-ID"])

    # Count unique genes per (Region_ID, Species)
    counts = (df.groupby(["Region_ID", "Species"])["Gene-ID"]
                .nunique()
                .unstack(fill_value=0))
    
    # remove the columns-axis name ("Species") so it doesn't print
    counts.columns.name = None

    # Ensure both species are present as columns
    for sp in (query, target):
        if sp not in counts.columns:
            counts[sp] = 0

    # Build final table
    out = counts[[query, target]].astype(int).rename(
        columns={query: "n_genes_query", target: "n_genes_target"}
    ).reset_index()

    out = (out
             .sort_values('Region_ID',
                          key=lambda s: pd.to_numeric(s, errors='coerce'))
             .reset_index(drop=True))

    return out


def check_results(df_counts, query, target, input_regions_table, output_plot_length_vs_genes_query, output_plot_length_vs_genes_target):
    """
    Calculates and prints distribution statistics and creates scatter plots.
    
    Input:
    - df_counts: DataFrame from count_genes() with columns Region_ID, n_genes_query, n_genes_target
    - query: str, full species name for query
    - target: str, full species name for target
    - input_regions_table: str, path to TSV with Region_ID, Query_Start, Query_End, Target_Start, Target_End
    - output_plot_length_vs_genes_query: str, output path for query scatter plot (PNG)
    - output_plot_length_vs_genes_target: str, output path for target scatter plot (PNG)
    
    Output:
    - Two PNG scatter plots (block length vs number of genes)
    
    Prints:
    - Total blocks count, mean/median/std dev/min/max genes per block, quartiles
    - Number and percentage of blocks with 0 genes
    - Median gene density per Mb
    """

    print(f"\nChecking {query}-{target}:\n")

    # 1. Basic distribution statistics
    print("DISTRIBUTION STATISTICS:")
    
    for species, col in [(query, 'n_genes_query'), (target, 'n_genes_target')]:
        data = df_counts[col]
        print(f"\n{species}:")
        print(f"Total blocks: {len(data)}")
        print(f"Mean genes/block: {data.mean():.2f}")
        print(f"Median genes/block: {data.median():.0f}")
        print(f"Std dev: {data.std():.2f}")
        print(f"Min: {data.min()}, Max: {data.max()}")
        
        # Quartiles
        q25, q50, q75 = data.quantile([0.25, 0.5, 0.75])
        print(f"Quartiles: Q1={q25:.0f}, Q2={q50:.0f}, Q3={q75:.0f}")
    
    # 2. Check Blocks with 0 genes
    print("\nBLOCKS WITH 0 GENES CHECK:\n")

    for col in ("n_genes_query", "n_genes_target"):
        df_counts[col] = pd.to_numeric(df_counts[col], errors="coerce").fillna(0).astype(int)

    def pct(n, d):
        n = float(n)
        d = float(d)
        return 0.0 if d == 0 else 100.0 * n / d

    # denominators
    total_blocks = int(df_counts.shape[0])

    # numerators
    zero_q = int((df_counts["n_genes_query"] == 0).sum())
    zero_t = int((df_counts["n_genes_target"] == 0).sum())

    print(f"Blocks with 0 genes in {query}:  {zero_q} ({pct(zero_q, total_blocks):.1f}%)")
    print(f"Blocks with 0 genes in {target}: {zero_t} ({pct(zero_t, total_blocks):.1f}%)")

    # 3. Check block length vs number of genes in block

    # df_counts as we have + a sizes table:
    sizes = pd.read_csv(input_regions_table, sep="\t")  # Region_ID, q_len, t_len

    # df_counts has Region_ID + gene counts
    df_counts["Region_ID"] = df_counts["Region_ID"].astype(str)
    sizes["Region_ID"] = sizes["Region_ID"].astype(str)

    d = df_counts.merge(sizes[["Region_ID", "Query_Start", "Query_End", "Target_Start", "Target_End"]], on="Region_ID", how="left")

    for c in ["Query_Start","Query_End","Target_Start","Target_End"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")

    q_len = abs(d["Query_End"] - d["Query_Start"])
    t_len = abs(d["Target_End"] - d["Target_Start"])

    # Plot block length vs number of genes per block
    plt.subplots(figsize=(8,5), dpi=150)
    plt.scatter(q_len/1e3, d["n_genes_query"], s=8, alpha=.4)
    plt.xlabel("Query block length (kb)"); plt.ylabel("# genes in block (query)")
    plt.title("Genes vs block length (query)", wrap=True)
    plt.savefig(output_plot_length_vs_genes_query)
    plt.show()
    plt.close()
    plt.subplots(figsize=(8,5), dpi=150)
    plt.scatter(t_len/1e3, d["n_genes_target"], s=8, alpha=.4)
    plt.xlabel("Target block length (kb)"); plt.ylabel("# genes in block (target)")
    plt.title("Genes vs block length (target)", wrap=True)
    plt.savefig(output_plot_length_vs_genes_target)
    plt.show()
    plt.close()

    # 4. Gene density per Mb
    print("\nGENE DENSITY PER Mb:\n")
    d["q_gpm"] = d["n_genes_query"] / (q_len.replace(0, np.nan) / 1e6)
    d["t_gpm"] = d["n_genes_target"] / (t_len.replace(0, np.nan) / 1e6)
    print("Median genes/Mb (query, target):", d["q_gpm"].median(skipna=True), d["t_gpm"].median(skipna=True))


def plot_gene_distribution(df_counts, query, target, plot_output_overall, plot_output_closer):
    """
    Create histogram plots of gene distribution across synteny blocks.
    
    Input:
    - df_counts: DataFrame with columns Region_ID, n_genes_query, n_genes_target
    - query: str, full species name for query
    - target: str, full species name for target
    - plot_output_overall: str, output path for overall histogram (PNG)
    - plot_output_closer: str, output path for detailed 0-100 histogram (PNG)
    
    Output:
    - Two PNG histogram plots (overall distribution and 0-100 detail)
    
    Prints: Nothing
    """
    
    # Histograms

    # Plot the overall gene distribution
    plt.subplots(figsize=(8,5), dpi=150) 
    plt.hist(df_counts["n_genes_query"], bins=40, log=False, alpha=0.5, label='Query', histtype="stepfilled",color="blue", edgecolor='black', linewidth=1.2)
    plt.hist(df_counts["n_genes_target"], bins=40, log=False, alpha=0.5, label='Target', histtype="stepfilled",color="orange", edgecolor='black', linewidth=1.2)
    plt.legend(loc='upper right')
    plt.xlabel("Number of genes per block")
    plt.ylabel("Count")
    plt.title(f"Gene count per synteny block, {query} (Query) vs {target} (Target), overall", wrap=True)
    plt.savefig(plot_output_overall)
    plt.show()
    plt.close()


    # Plot the closer plot to investigate the interval 0-100
    q = pd.to_numeric(df_counts["n_genes_query"], errors="coerce").fillna(0).astype(int).to_numpy()
    t = pd.to_numeric(df_counts["n_genes_target"], errors="coerce").fillna(0).astype(int).to_numpy()

    kmax = int(max(q.max() if q.size else 0, t.max() if t.size else 0))

    bins = np.concatenate([
        np.arange(-0.5, 100.5, 1),           # 0..100 (unit bins)
        np.arange(100.5, kmax + 10.5, 10)    # 100+ (10-wide bins)
    ])

    plt.subplots(figsize=(10,7), dpi=150)
    plt.hist(q, bins=bins, log=False, alpha=0.5, label='Query', histtype="stepfilled",color="blue", edgecolor='black', linewidth=1.2)
    plt.hist(t, bins=bins, log=False, alpha=0.5, label='Target', histtype="stepfilled",color="orange", edgecolor='black', linewidth=1.2)
    #plt.yscale("log")
    plt.xlim(-0.5, 120.5)                   # show detailed zone clearly
    plt.xticks(np.arange(0, 121, 5))
    plt.legend(loc='upper right')
    plt.xlabel("Number of genes per block")
    plt.ylabel("Count")
    plt.title(f"Gene count per synteny block, {query} (Query) vs {target} (Target), detailed look at the interval 0-100", wrap=True)
    plt.savefig(plot_output_closer)
    plt.show()
    plt.close()


def plot_total_number_genes(number_genes_ann, genes_in_blocks, source, target, output_plot):
    """
    Create pie chart showing proportion of genes in synteny blocks.
    
    Input:
    - number_genes_ann: int, total number of annotated genes
    - genes_in_blocks: int, number of genes found in synteny blocks
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
    ratio = genes_in_blocks / number_genes_ann * 100

    plt.figure(figsize=(5,5), dpi=150)
    plt.pie(
        sizes,
        labels=labels,
        autopct='%1.1f%%',
        startangle=90,
        colors=colors,
        wedgeprops={'edgecolor': 'black'}
    )
    plt.title(f"{source}->{target}: {ratio:.1f}% of genes in synteny blocks")
    plt.savefig(output_plot)
    plt.close()


def process_pairs(source, target, type_species):
    """
    Runs the complete analysis pipeline for one species pair.
    
    Input:
    - source: str, source species abbreviation (e.g., "can")
    - target: str, target species abbreviation (e.g., "mac")
    - type_species: str, taxonomic group (e.g., "mammalia")
    
    Output:
    - TSV file with gene counts per block
    - 4 PNG plots (length vs genes for query/target, distribution overview, distribution 0-100)
    
    Prints:
    - Distribution statistics via check_results()
    """

    
    # --- inputs ---
    genes_tsv   = f"results/genes_in_blocks/{type_species}/region_gene_tables/{source}_{target}_synteny_tables_genes_fix.tsv"
    regions_tsv = f"results/genes_in_blocks/{type_species}/region_gene_tables/{source}_{target}_synteny_tables_regions_fix.tsv"
    

    # --- outputs ---
    output_dir = f"results/data_check/{type_species}/gene_distribution_in_blocks"
    os.makedirs(output_dir, exist_ok=True)
    
    counts_tsv      = f"{output_dir}/{source}_{target}_number_genes_in_blocks.tsv"
    len_vs_genes_q  = f"{output_dir}/{source}_{target}_length_vs_gene_number_query.png"
    len_vs_genes_t  = f"{output_dir}/{source}_{target}_length_vs_gene_number_target.png"
    dist_overview   = f"{output_dir}/{source}_{target}_gene_distribution_overview.png"
    dist_0_100      = f"{output_dir}/{source}_{target}_gene_distribution_0_100.png"
    
    full_query = species_full[source]
    full_target = species_full[target]


    # 1) count & save
    df_counts = count_genes(str(genes_tsv), full_query, full_target)
    df_counts.to_csv(counts_tsv, sep="\t", index=False)

    # 2) data check plots
    check_results(
        df_counts,
        full_query, full_target,
        str(regions_tsv),
        str(len_vs_genes_q),
        str(len_vs_genes_t),
    )

    # 3) distribution plots
    plot_gene_distribution(
        df_counts,
        full_query, full_target,
        str(dist_overview),
        str(dist_0_100),
    )
    

def main():
    for src, tgt, type_species in pairs:
        try:
            process_pairs(src, tgt, type_species)
        except Exception as e:
            print(f"[ERROR] {src}-{tgt}: {e}")

if __name__ == "__main__":
    main()