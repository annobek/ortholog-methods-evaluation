import pandas as pd
import matplotlib.pyplot as plt
import pyranges as pr
import os
import numpy as np
import traceback

"""
Script: Genome Coverage by Synteny Blocks

Purpose: Calculates what percentage of each genome is covered by synteny blocks for all species pairs.

Input:
- TSV files with synteny block coordinates (*_synteny_tables_regions_fix.tsv)
- GTF files with genome annotations (gtf_*.gtf)

Output:
- Bar plots showing genome coverage per pair (*_whole_genome_coverage_bar_fix.png)
- Summary TSV with all coverage percentages (summary_points_fix.tsv)
- Summary boxplot across all pairs (summary_boxplot_fix.png)

Prints:
- Processing status for each pair
- Block lengths (forward and reverse)
- Total genome bp for each species
- Coverage percentages for each species and direction
- Number of collected data points
- Summary statistics for boxplot
"""

# ============================================================
# CONFIG
# ============================================================

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

annotation_names = ["Canis", "Macaca", "Rat", "Mus"]

annotation_map = {
    "can": "Canis",
    "mac": "Macaca",
    "rat": "Rat",
    "mus": "Mus",
}

OUT_DIR = "results/data_check/coverage/whole_genome"
os.makedirs(OUT_DIR, exist_ok=True)
METRIC_NAME = "coverage_genome"
all_metrics = []  # to collect per-pair values for summary boxplot


# ============================================================
# HELPERS
# ============================================================

def unique_bp(pr_obj):
    """
    Calculate non-overlapping bp covered by genomic intervals.
    
    Input:
    - pr_obj: PyRanges object with genomic intervals
    
    Output:
    - int, total bp covered (after merging overlaps)
    
    Prints: Nothing
    """
    return 0 if pr_obj.empty else int(pr_obj.merge().lengths().sum())


def intersection_bp(pr_a, pr_b):
    """
    Calculate bp in common between two sets of genomic intervals.
    
    Input:
    - pr_a: PyRanges object with first set of intervals
    - pr_b: PyRanges object with second set of intervals
    
    Output:
    - int, total bp of intersection
    
    Prints: Nothing
    """
    j = pr_a.join(pr_b)
    if j.empty:
        return 0
    df = j.df
    inter = (df[["End", "End_b"]].min(axis=1) - df[["Start", "Start_b"]].max(axis=1)).clip(lower=0)
    return int(inter.sum())


def get_synteny_intervals(region_tables):
    """
    Read synteny region tables and convert to PyRanges interval objects.
    
    Input:
    - region_tables: list of str, paths to TSV files with columns Region_ID, Query_Scaffold, Query_Start, Query_End, Target_Scaffold, Target_Start, Target_End
    
    Output:
    - Tuple of (synteny_pr_query_array: list of PyRanges, synteny_pr_target_array: list of PyRanges)
    
    Prints:
    - Processing status for each region table
    """
    synteny_pr_query_array = []
    synteny_pr_target_array = []

    for region_table in region_tables:
        print(f"Processing synteny: {region_table}\n")

        region = pd.read_csv(region_table, sep="\t")

        synt_interval_query = region[["Region_ID", "Query_Scaffold", "Query_Start", "Query_End"]]
        synt_interval_target = region[["Region_ID", "Target_Scaffold", "Target_Start", "Target_End"]]

        synt_interval_query = synt_interval_query.rename(columns={
            "Query_Scaffold": "Chromosome", "Query_Start": "Start", "Query_End": "End"})
        synt_interval_target = synt_interval_target.rename(columns={
            "Target_Scaffold": "Chromosome", "Target_Start": "Start", "Target_End": "End"})

        for df in (synt_interval_query, synt_interval_target):
            df["Chromosome"] = df["Chromosome"].astype(str)
            df["Start"] = df["Start"].astype(int)
            df["End"] = df["End"].astype(int)

        synteny_pr_query_array.append(pr.PyRanges(synt_interval_query))
        synteny_pr_target_array.append(pr.PyRanges(synt_interval_target))

    return synteny_pr_query_array, synteny_pr_target_array


def get_block_length(synt_interval_query, synt_interval_target):
    """
    Compute total bp of all synteny blocks (including overlaps).
    
    Input:
    - synt_interval_query: list of PyRanges, query species intervals
    - synt_interval_target: list of PyRanges, target species intervals
    
    Output:
    - Tuple of (blq: list of int, blt: list of int), total bp for query and target
    
    Prints: Nothing
    """
    blq, blt = [], []
    for q, t in zip(synt_interval_query, synt_interval_target):
        blq.append((q.df["End"] - q.df["Start"]).sum())
        blt.append((t.df["End"] - t.df["Start"]).sum())
    return blq, blt


def count_genome_length(annotation_path, extra_pr_list=None):
    """
    Estimate total genome size (bp) using annotation + synteny block maxima.
    
    Input:
    - annotation_path: str, path to genome annotation file (GTF, GFF, BED, or TSV format, optionally gzipped)
    - extra_pr_list: list of PyRanges or None, additional intervals to consider for chromosome max positions
    
    Output:
    - int, estimated total genome length in bp
    
    Prints: Nothing
    """
    path = str(annotation_path)
    lower = path.lower()

    # --- choose parser by extension (handles .gz) ---
    if lower.endswith((".gtf", ".gtf.gz")):
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["Chromosome", "Source", "Feature", "Start", "End",
                   "Score", "Strand", "Frame", "Attributes"],
            comment="#",
            engine="python",
        )
        df = df[["Chromosome", "End"]].dropna()
    elif lower.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")):
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["Chromosome", "Source", "Feature", "Start", "End",
                   "Score", "Strand", "Phase", "Attributes"],
            comment="#",
            engine="python",
        )
        df = df[["Chromosome", "End"]].dropna()
    elif lower.endswith((".bed", ".bed.gz")):
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            usecols=[0, 2],
            names=["Chromosome", "End"],
            comment="#",
            engine="python",
        ).dropna()
    else:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            usecols=[0, 2],
            names=["Chromosome", "End"],
            comment="#",
            engine="python",
        ).dropna()

    df["Chromosome"] = df["Chromosome"].astype(str)
    df["End"] = pd.to_numeric(df["End"], errors="coerce").dropna().astype(int)
    chrom_max = df.groupby("Chromosome")["End"].max().astype(int).to_dict()

    if extra_pr_list:
        for pr_obj in extra_pr_list:
            if pr_obj is None or pr_obj.empty:
                continue
            dpr = pr_obj.df.groupby("Chromosome")["End"].max().astype(int).to_dict()
            for chrom, m in dpr.items():
                chrom_max[chrom] = max(chrom_max.get(chrom, 0), int(m))

    return int(sum(chrom_max.values()))


def plot_bp_coverage(direction_label, species_labels, bp_values, total_bp_dict, out_png):
    """
    Create bar plot showing genome coverage percentage for a species pair.
    
    Input:
    - direction_label: str, label for the pair direction (e.g., "can->mac (genome)")
    - species_labels: list of str, species abbreviations for x-axis
    - bp_values: list of int, bp covered for each species
    - total_bp_dict: dict, mapping species abbreviation to total genome bp
    - out_png: str, output path for bar plot (PNG)
    
    Output:
    - PNG bar plot
    
    Prints: Nothing
    """
    pct = [(bp / total_bp_dict[s]) * 100 for s, bp in zip(species_labels, bp_values)]

    plt.figure(figsize=(6, 4), dpi=150)
    plt.bar(species_labels, pct, color=["#66b3ff", "#ff9999"], edgecolor="black")
    plt.ylabel("% genome in synteny blocks")
    plt.title(f"Synteny block coverage: {direction_label}")
    plt.ylim(0, 100)
    for i, p in enumerate(pct):
        plt.text(i, p + 1, f"{p:.1f}%", ha="center", fontsize=10)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


# ============================================================
# BOX PLOT SUMMARY (same as conservativeness)
# ============================================================

def finalize_boxplot(all_metrics, metric_name, out_csv, out_png):
    """
    Create summary boxplot showing distribution of genome coverage across all pairs.
    
    Input:
    - all_metrics: list of dict, each with keys: metric, pair_dir, species, value
    - metric_name: str, name of the metric to plot
    - out_csv: str, output path for TSV with all data points
    - out_png: str, output path for boxplot (PNG)
    
    Output:
    - TSV file with all coverage percentages
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

    plt.figure(figsize=(9, 7), dpi=150)
    plt.boxplot(df_m["value"], vert=True, showmeans=True, widths=0.35, labels=[metric_name])

    # jittered points with direction-based colors
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

    '''
    # Label lowest few points (bottom 3)
    label_n = 3
    low_idx = np.argsort(df_m["value"])[:label_n]
    for i in low_idx:
        xi = x[i]
        val = df_m["value"].iloc[i]
        pair = df_m["pair_dir"].iloc[i]
        plt.annotate(pair, xy=(xi, val), xytext=(xi + 0.15, val - 0.5),
                     fontsize=9, ha="left", va="center",
                     arrowprops=dict(arrowstyle="->", color="gray", lw=0.8, alpha=0.6),
                     path_effects=[pe.withStroke(linewidth=3, foreground="white")])
    '''

    plt.ylim(0, 100)
    plt.ylabel("% genome in synteny blocks", fontsize=12)
    plt.title(f"{metric_name}: summary across all query→target bars (n={n})", fontsize=13)
    plt.grid(axis="y", linestyle="--", alpha=0.4, zorder=0)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


# ============================================================
# CORE
# ============================================================

def process_pairs(query, target):
    """
    Compute genome coverage for one species pair in both directions.
    
    Input:
    - query: str, query species abbreviation (e.g., "can")
    - target: str, target species abbreviation (e.g., "mac")
    
    Output:
    - PNG bar plot showing coverage for this pair
    - Appends coverage percentages to global all_metrics list
    
    Prints:
    - Processing status
    - Block lengths for forward and reverse directions
    - Total genome bp for each species
    - Coverage percentages for each species and direction
    - Output file path
    """
    fwd_regions = f"results/genes_in_blocks/mammalia/region_gene_tables/{query}_{target}_synteny_tables_regions_fix.tsv"
    rev_regions = f"results/genes_in_blocks/mammalia/region_gene_tables/{target}_{query}_synteny_tables_regions_fix.tsv"
    annotation_src_whole = f"material/sex_experiment/gtf_{annotation_map[query]}.gtf"
    annotation_tgt_whole = f"material/sex_experiment/gtf_{annotation_map[target]}.gtf"

    barplot_whole_genome = f"{OUT_DIR}/{query}_{target}_whole_genome_coverage_bar_fix.png"

    print(f"\n{'='*60}")
    print(f"Processing pair {query} -> {target}")
    print(f"{'='*60}\n")

    # Extract intervals
    interval_query_fwd, interval_target_fwd = get_synteny_intervals([fwd_regions])
    interval_query_rev, interval_target_rev = get_synteny_intervals([rev_regions])

    blq_fwd, blt_fwd = get_block_length(interval_query_fwd, interval_target_fwd)
    blq_rev, blt_rev = get_block_length(interval_query_rev, interval_target_rev)
    print(f"Block lengths:")
    print(f"Forward - Query: {blq_fwd[0]}, Target: {blt_fwd[0]}")
    print(f"Reverse - Query: {blq_rev[0]}, Target: {blt_rev[0]}\n")

    q_fwd = interval_query_fwd[0].merge()
    t_fwd = interval_target_fwd[0].merge()
    q_rev = interval_query_rev[0].merge()
    t_rev = interval_target_rev[0].merge()

    total_genome_bp_src = count_genome_length(annotation_src_whole, extra_pr_list=[q_fwd, t_rev])
    total_genome_bp_tgt = count_genome_length(annotation_tgt_whole, extra_pr_list=[t_fwd, q_rev])
    total_genome_bp_dict = {query: total_genome_bp_src, target: total_genome_bp_tgt}

    print(f"Denominators - genome bp: {query}={total_genome_bp_src:,}, {target}={total_genome_bp_tgt:,}\n")

    # bp coverage
    sp1_bp_fwd = unique_bp(q_fwd)
    sp2_bp_fwd = unique_bp(t_fwd)
    sp2_bp_rev = unique_bp(q_rev)
    sp1_bp_rev = unique_bp(t_rev)

    # Percent coverage
    pct_sp1_fwd_genome = sp1_bp_fwd / total_genome_bp_src * 100
    pct_sp2_fwd_genome = sp2_bp_fwd / total_genome_bp_tgt * 100
    pct_sp1_rev_genome = sp1_bp_rev / total_genome_bp_src * 100
    pct_sp2_rev_genome = sp2_bp_rev / total_genome_bp_tgt * 100

    print("Coverage (% of genome bp):")
    print(f"{query}->{target}: {query} {pct_sp1_fwd_genome:.2f}%, {target} {pct_sp2_fwd_genome:.2f}%")
    print(f"{target}->{query}: {target} {pct_sp2_rev_genome:.2f}%, {query} {pct_sp1_rev_genome:.2f}%\n")

    # Plot per-pair bar chart
    plot_bp_coverage(f"{query}->{target} (genome)",
                     [query, target],
                     [sp1_bp_fwd, sp2_bp_fwd],
                     total_genome_bp_dict,
                     barplot_whole_genome)

    print(f"Plot saved to: {barplot_whole_genome}\n")

    # Add metrics for boxplot
    all_metrics.append({"metric": METRIC_NAME, "pair_dir": f"{query}->{target}", "species": query, "value": pct_sp1_fwd_genome})
    all_metrics.append({"metric": METRIC_NAME, "pair_dir": f"{query}->{target}", "species": target, "value": pct_sp2_fwd_genome})


# ============================================================
# MAIN
# ============================================================

def main():
    """
    Main function: process all species pairs and create summary visualizations.
    
    Prints:
    - Error messages for any failed pairs
    - Total number of collected data points vs expected
    - Summary statistics from finalize_boxplot()
    """
    for src, tgt in pairs:
        try:
            process_pairs(src, tgt)
        except Exception as e:
            print(f"[ERROR] {src}-{tgt}: {e}")
            traceback.print_exc()

    expected = 2 * len(pairs)
    print(f"\nCollected {len(all_metrics)} bar values; expected {expected}.\n")

    finalize_boxplot(
        all_metrics,
        metric_name=METRIC_NAME,
        out_csv=os.path.join(OUT_DIR, "summary_points_fix.tsv"),
        out_png=os.path.join(OUT_DIR, "summary_boxplot_fix.png"),
    )


if __name__ == "__main__":
    main()