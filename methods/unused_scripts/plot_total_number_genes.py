import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import traceback
import glob

# ==============================================
# Configuration
# ==============================================

'''
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
    ("ar", "card", "plants"),
    ("card", "ar", "plants"),
    ("dmel", "dsec", "drosophila"),
    ("dsec", "dmel", "drosophila"),
    ("dmel", "dsim", "drosophila"),
    ("dsim", "dmel", "drosophila"),
    ("dsec", "dsim", "drosophila"),
    ("dsim", "dsec", "drosophila")
]
'''
pairs = [
    ("ar", "card", "plants"),
    ("card", "ar", "plants"),
    ("dmel", "dsec", "drosophila"),
    ("dsec", "dmel", "drosophila"),
    ("dmel", "dsim", "drosophila"),
    ("dsim", "dmel", "drosophila"),
    ("dsec", "dsim", "drosophila"),
    ("dsim", "dsec", "drosophila")
]



species_full = {
    "can": "Canis lupus familiaris",
    "mac": "Macaca fascicularis",
    "rat": "Rattus norvegicus",
    "mus": "Mus musculus",
    "ar": "Arabidopsis thaliana",
    "card": "Cardamine hirsuta",
    "dmel": "Drosophila melanogaster",
    "dsec": "Drosophila sechellia",
    "dsim": "Drosophila simulans"
}

annotation_map = {
    "can": "Canis",
    "mac": "Macaca",
    "rat": "Rat",
    "mus": "Mus",
    "ar": "Arabidopsis",
    "card": "Cardamine",
    "dmel": "dmel",
    "dsec": "dsec",
    "dsim": "dsim"
}

METRIC_NAME = "genes_in_blocks_pct"
all_metrics = []


# ==============================================
# Helper functions
# ==============================================

def smart_find_annotation(base_dir, species_code):
    """Try to find an annotation file for a species within a group directory."""
    name = annotation_map[species_code]
    patterns = [
        f"{name}*genes_only_stableid_fix.bed6",  # mammals
        f"{name}.gff3", f"{name}.gff",           # plants
        f"{name}.gtf"                            # flies
    ]
    for pat in patterns:
        matches = glob.glob(os.path.join(base_dir, pat))
        if matches:
            return matches[0]
    raise FileNotFoundError(f"No annotation file found for {species_code} in {base_dir}")


def count_genes_int(gene_table_path, source, target):
    """Count genes in synteny blocks from synteny table."""
    df = pd.read_csv(gene_table_path, sep="\t", dtype=str)
    df = df.dropna(subset=["Region_ID", "Species", "Gene-ID"]).drop_duplicates(subset=["Region_ID", "Species", "Gene-ID"])
    n_source = df[df["Species"].str.contains(annotation_map[source], case=False, na=False)]["Gene-ID"].nunique()
    n_target = df[df["Species"].str.contains(annotation_map[target], case=False, na=False)]["Gene-ID"].nunique()
    return n_source, n_target


def count_genes_annotation(path):
    """Count total genes in an annotation, automatically detecting format."""
    ext = os.path.splitext(path)[1].lower()

    if ext == ".bed6" or ext == ".bed":
        cols = ["chrom", "start", "end", "GeneID", "score", "strand"]
        df = pd.read_csv(path, sep="\t", names=cols, dtype=str)
        return df["GeneID"].nunique()

    elif ext in [".gff", ".gff3"]:
        df = pd.read_csv(path, sep="\t", comment="#", header=None, dtype=str)
        if df.shape[1] < 9:
            return 0
        df.columns = ["seqid","source","type","start","end","score","strand","phase","attributes"]
        genes = df[df["type"].str.lower()=="gene"].copy()
        genes["GeneID"] = genes["attributes"].str.extract(r'ID=([^;]+)')
        return genes["GeneID"].nunique()

    elif ext == ".gtf":
        df = pd.read_csv(path, sep="\t", comment="#", header=None, dtype=str)
        if df.shape[1] < 9:
            return 0
        df.columns = ["seqname","source","feature","start","end","score","strand","frame","attributes"]
        genes = df[df["feature"].str.lower()=="gene"].copy()
        genes["GeneID"] = genes["attributes"].str.extract(r'gene_id "([^"]+)"')
        return genes["GeneID"].nunique()

    else:
        print(f"[WARN] Unrecognized file format: {path}")
        return 0


def plot_total_number_genes(total_genes, in_blocks, role, source, target, out_path):
    """Piechart for a single species."""
    missing = total_genes - in_blocks
    sizes = [in_blocks, missing]
    labels = ["In synteny blocks", "Not in synteny blocks"]
    ratio = in_blocks / total_genes * 100 if total_genes else 0

    plt.figure(figsize=(5,5), dpi=150)
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90,
            colors=["#66b3ff","#ff9999"], wedgeprops={'edgecolor':'black'})
    plt.title(f"{source}->{target} ({role}): {ratio:.1f}% in synteny blocks")
    plt.savefig(out_path)
    plt.close()


def finalize_boxplot(df, metric_name, out_csv, out_png, title_suffix=""):
    """Boxplot summarizing all species results."""
    if df.empty:
        return
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df.to_csv(out_csv, sep="\t", index=False)

    plt.figure(figsize=(9,7), dpi=150)
    plt.boxplot(df["value"], vert=True, showmeans=True, widths=0.35, tick_labels=[metric_name])

    rng = np.random.default_rng(42)
    x = rng.normal(1, 0.025, size=len(df))
    colors = df["pair_dir"].apply(lambda p: "tab:blue" if p.split("->")[0] < p.split("->")[1] else "tab:orange")
    plt.scatter(x, df["value"], c=colors, s=35, alpha=0.9, edgecolors="black", linewidths=0.4, zorder=3)
    plt.scatter([], [], c="tab:blue", label="A->B")
    plt.scatter([], [], c="tab:orange", label="B->A")
    plt.legend(frameon=False, loc="lower left", fontsize=9)

    plt.ylim(0, 100)
    plt.ylabel("% of genes in synteny blocks", fontsize=12)
    plt.title(f"{metric_name} summary {title_suffix}", fontsize=13)
    plt.grid(axis="y", linestyle="--", alpha=0.4, zorder=0)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def process_pair(src, tgt, group):
    print(f"\n=== Processing {src}->{tgt} ({group}) ===")

    # Input files
    gene_table = f"results/genes_in_blocks/{group}/region_gene_tables/{src}_{tgt}_synteny_tables_genes.tsv"
    ann_dir = f"results/ann_genes/{group}"
    ann_src = smart_find_annotation(ann_dir, src)
    ann_tgt = smart_find_annotation(ann_dir, tgt)

    # Output directory
    out_dir = f"results/data_check/{group}/total_number_genes"
    os.makedirs(out_dir, exist_ok=True)

    pie_src = f"{out_dir}/{src}_{tgt}_pie_{src}.png"
    pie_tgt = f"{out_dir}/{src}_{tgt}_pie_{tgt}.png"

    # Counts
    total_src = count_genes_annotation(ann_src)
    total_tgt = count_genes_annotation(ann_tgt)
    in_blocks_src, in_blocks_tgt = count_genes_int(gene_table, src, tgt)

    # Plots
    plot_total_number_genes(total_src, in_blocks_src, "query", src, tgt, pie_src)
    plot_total_number_genes(total_tgt, in_blocks_tgt, "target", src, tgt, pie_tgt)

    # Save results
    pct_src = (in_blocks_src / total_src * 100) if total_src else 0
    pct_tgt = (in_blocks_tgt / total_tgt * 100) if total_tgt else 0
    all_metrics.append({"group": group, "metric": METRIC_NAME, "pair_dir": f"{src}->{tgt}", "species": src, "value": pct_src})
    all_metrics.append({"group": group, "metric": METRIC_NAME, "pair_dir": f"{src}->{tgt}", "species": tgt, "value": pct_tgt})


# ==============================================
# Main
# ==============================================

def main():
    for src, tgt, group in pairs:
        try:
            process_pair(src, tgt, group)
        except Exception as e:
            print(f"[ERROR] {src}-{tgt} ({group}): {e}")
            traceback.print_exc()

    df_all = pd.DataFrame(all_metrics)
    print(f"\nCollected {len(df_all)} results from {len(pairs)} pairs.")

    # Global boxplot
    finalize_boxplot(
        df_all[df_all["metric"] == METRIC_NAME],
        METRIC_NAME,
        out_csv="results/data_check/summary_total_genes_all.tsv",
        out_png="results/data_check/summary_total_genes_all.png",
        title_suffix="(all groups)"
    )

    # Per-group boxplots
    for grp in df_all["group"].unique():
        sub = df_all[df_all["group"] == grp]
        finalize_boxplot(
            sub,
            METRIC_NAME,
            out_csv=f"results/data_check/{grp}/summary_total_genes_{grp}.tsv",
            out_png=f"results/data_check/{grp}/summary_total_genes_{grp}.png",
            title_suffix=f"({grp})"
        )


if __name__ == "__main__":
    main()
