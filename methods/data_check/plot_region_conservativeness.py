import pandas as pd
import matplotlib.pyplot as plt
import pyranges as pr
import traceback
import os
import numpy as np

# =========================
# CONFIG
# =========================
THRESH = 0.90
METRIC_NAME = f"conservativeness_{int(THRESH*100)}"
OUT_DIR = f"results/data_check/conservativeness/overlap{int(THRESH*100)}"
os.makedirs(OUT_DIR, exist_ok=True)

all_metrics = []

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
#pairs = [("can", "mac"), ("mac", "can")]

# =========================
# HELPERS
# =========================
def _normalize_species_label(pr_obj, _label_unused):
    df = pr_obj.df.copy()
    df["Chromosome"] = df["Chromosome"].astype(str)
    df["Start"] = df["Start"].astype(int)
    df["End"]   = df["End"].astype(int)
    return pr.PyRanges(df)

def _pairwise_overlap(pr_left, pr_right):
    j = pr_left.join(pr_right)
    df = j.df.copy()
    df["int_len"] = (df[["End","End_b"]].min(axis=1) -
                     df[["Start","Start_b"]].max(axis=1)).clip(lower=0)
    return df

def get_synteny_intervals(region_tables):
    q_arr, t_arr = [], []
    for region_table in region_tables:
        print(f"Processing synteny: {region_table}\n")
        region = pd.read_csv(region_table, sep="\t")

        q = region[["Region_ID","Query_Scaffold","Query_Start","Query_End"]].rename(
            columns={"Query_Scaffold":"Chromosome","Query_Start":"Start","Query_End":"End"}
        )
        t = region[["Region_ID","Target_Scaffold","Target_Start","Target_End"]].rename(
            columns={"Target_Scaffold":"Chromosome","Target_Start":"Start","Target_End":"End"}
        )
        for df in (q, t):
            df["Chromosome"] = df["Chromosome"].astype(str)
            df["Start"] = df["Start"].astype(int)
            df["End"] = df["End"].astype(int)

        q_arr.append(pr.PyRanges(q))
        t_arr.append(pr.PyRanges(t))
    return q_arr, t_arr

def get_block_length(q_arr, t_arr):
    blq, blt = [], []
    for q, t in zip(q_arr, t_arr):
        blq.append((q.df["End"] - q.df["Start"]).sum())
        blt.append((t.df["End"] - t.df["Start"]).sum())
    return blq, blt

def identify_reciprocal_pairs(q_fwd, t_fwd, q_rev, t_rev, thr=0.9):
    def lens(pr_obj):
        d = pr_obj.df.copy()
        d["len"] = d["End"] - d["Start"]
        return d.groupby("Region_ID")["len"].sum().to_dict()

    L_qf, L_tf, L_qr, L_tr = map(lens, (q_fwd, t_fwd, q_rev, t_rev))

    sp1_df = _pairwise_overlap(q_fwd, t_rev)   # fwd query vs rev target
    sp2_df = _pairwise_overlap(t_fwd, q_rev)   # fwd target vs rev query

    sp1 = (sp1_df.groupby(["Region_ID","Region_ID_b"])["int_len"]
                 .sum().rename("ovl_sp1").reset_index())
    sp2 = (sp2_df.groupby(["Region_ID","Region_ID_b"])["int_len"]
                 .sum().rename("ovl_sp2").reset_index())

    pairs = sp1.merge(sp2, on=["Region_ID","Region_ID_b"], how="outer").fillna(0)

    pairs["len_fwd_query"] = pairs["Region_ID"].map(L_qf)
    pairs["len_fwd_target"] = pairs["Region_ID"].map(L_tf)
    pairs["len_rev_query"]  = pairs["Region_ID_b"].map(L_qr)
    pairs["len_rev_target"] = pairs["Region_ID_b"].map(L_tr)

    pairs["fwd_query_frac"] = pairs["ovl_sp1"] / pairs["len_fwd_query"]
    pairs["fwd_target_frac"] = pairs["ovl_sp2"] / pairs["len_fwd_target"]
    pairs["rev_target_frac"] = pairs["ovl_sp1"] / pairs["len_rev_target"]
    pairs["rev_query_frac"]  = pairs["ovl_sp2"] / pairs["len_rev_query"]

    mask = (
        (pairs["fwd_query_frac"] >= thr) &
        (pairs["fwd_target_frac"] >= thr) &
        (pairs["rev_query_frac"]  >= thr) &
        (pairs["rev_target_frac"] >= thr)
    )
    reciprocal_pairs = pairs.loc[mask].copy()

    return reciprocal_pairs, set(reciprocal_pairs["Region_ID"]), set(reciprocal_pairs["Region_ID_b"])

def count_genes_in_blocks(gene_table_path, source, target):
    df = (pd.read_csv(gene_table_path, sep="\t", dtype=str)
            .dropna(subset=["Region_ID","Species","Gene-ID"])
            .drop_duplicates(subset=["Region_ID","Species","Gene-ID"]))
    n_source = df[df["Species"] == source]["Gene-ID"].nunique()
    n_target = df[df["Species"] == target]["Gene-ID"].nunique()
    return n_source, n_target

def count_genes_in_reciprocal_ids_single(df_genes, ids, species):
    df = df_genes.copy()
    df["Region_ID"] = df["Region_ID"].astype(str)
    keep = set(map(str, ids))
    df = df[df["Region_ID"].isin(keep)]
    return df[df["Species"] == species]["Gene-ID"].nunique()

def compute_ratios_per_species(n_query_recipr, n_target_recipr, n_query_all_fwd, n_target_all_fwd):
    r1 = n_query_recipr / n_query_all_fwd if n_query_all_fwd > 0 else 0.0
    r2 = n_target_recipr / n_target_all_fwd if n_target_all_fwd > 0 else 0.0
    return r1, r2

def plot_conservativeness_bar(r1, r2, query, target, output_plot, thr):
    labels = [query, target]
    values = [r1 * 100, r2 * 100]
    plt.figure(figsize=(11,9), dpi=150)
    plt.bar(labels, values, color=['steelblue', 'violet'], edgecolor='black')
    #plt.ylabel("% of genes in reciprocal synteny")
    plt.ylabel("% of synteny length reciprocally conserved")
    plt.title(f"Reciprocal coverage: {query} → {target}")
    plt.suptitle(f"overlap threshold = {int(thr*100)}%")
    plt.ylim(0, 100)
    for i, v in enumerate(values):
        plt.text(i, min(98, v + 1.5), f"{v:.1f}%", ha='center', fontsize=10)
    plt.tight_layout()
    plt.savefig(output_plot)
    plt.close()

def finalize_boxplot(all_metrics, metric_name, out_csv, out_png):
    #import matplotlib.patheffects as pe

    df = pd.DataFrame(all_metrics)
    df_m = df[df["metric"] == metric_name].copy()
    if df_m.empty:
        print(f"No rows for metric={metric_name}. Nothing to plot.")
        return

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df_m.to_csv(out_csv, sep="\t", index=False)

    n = len(df_m)
    n_pairs = df_m["pair_dir"].nunique()
    print(f"Boxplot will use {n} bars from {n_pairs} query→target pair-dirs.")

    plt.figure(figsize=(9, 7), dpi=150)
    plt.boxplot(df_m["value"], vert=True, showmeans=True, widths=0.35, labels=[metric_name])

    # jittered points with direction-based colors
    rng = np.random.default_rng(42)
    x = rng.normal(1, 0.025, size=n)

    # Assign color based on lexicographic order of species names
    colors = df_m["pair_dir"].apply(
        lambda p: "tab:blue" if "->" in p and p.split("->")[0] < p.split("->")[1] else "tab:orange"
    )

    plt.scatter(
        x, df_m["value"],
        c=colors, s=35, alpha=0.9, edgecolors="black", linewidths=0.4, zorder=3
    )

    # Add legend for direction meaning
    plt.scatter([], [], c="tab:blue", label="A->B")
    plt.scatter([], [], c="tab:orange", label="B->A")
    plt.legend(frameon=False, loc="lower left", fontsize=9)

    '''
    # label low outliers (<95%)
    # label low outliers (<95%) with staggered offsets
    y_min_label = 95
    label_positions = []
    for xi, val, pair in zip(x, df_m["value"], df_m["pair_dir"]):
        if val < y_min_label:
            # check for existing nearby labels and stagger if needed
            y_offset = 0
            for _, existing_y in label_positions:
                if abs(val - existing_y) < 0.4:
                    y_offset -= 0.5  # move a bit lower to avoid overlap

            plt.annotate(
                pair,
                xy=(xi, val),
                xytext=(xi + 0.15, val + y_offset),
                fontsize=9,
                ha="left",
                va="center",
                arrowprops=dict(arrowstyle="->", color="gray", lw=0.8, alpha=0.6),
                path_effects=[pe.withStroke(linewidth=3, foreground="white")],
            )
            label_positions.append((xi, val + y_offset))
    '''


    plt.ylim(80, 100)
    plt.ylabel("% of synteny length reciprocally conserved", fontsize=12)
    plt.title(f"{metric_name}: summary across all query->target bars (n={n})", fontsize=13)
    plt.grid(axis="y", linestyle="--", alpha=0.4, zorder=0)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()




# =========================
# CORE
# =========================
def process_pairs(query, target, thr=THRESH):
    fwd_regions = f"results/genes_in_blocks/mammalia/region_gene_tables/{query}_{target}_synteny_tables_regions_fix.tsv"
    rev_regions = f"results/genes_in_blocks/mammalia/region_gene_tables/{target}_{query}_synteny_tables_regions_fix.tsv"
  

    # ensure pair output dir exists
    os.makedirs(OUT_DIR, exist_ok=True)
    barplot = os.path.join(OUT_DIR, f"{query}_{target}_conservativeness_bar_overlap{int(thr*100)}_fix.png")

    print(f"\n{'='*60}")
    print(f"Processing pair {query} -> {target}")
    print(f"{'='*60}\n")

    interval_query_fwd, interval_target_fwd = get_synteny_intervals([fwd_regions])
    interval_query_rev, interval_target_rev = get_synteny_intervals([rev_regions])

    blq_fwd, blt_fwd = get_block_length(interval_query_fwd, interval_target_fwd)
    blq_rev, blt_rev = get_block_length(interval_query_rev, interval_target_rev)
    print("Block lengths:")
    print(f"  Forward - Query: {blq_fwd[0]}, Target: {blt_fwd[0]}")
    print(f"  Reverse - Query: {blq_rev[0]}, Target: {blt_rev[0]}\n")

    q_fwd_n = _normalize_species_label(interval_query_fwd[0], query)
    t_fwd_n = _normalize_species_label(interval_target_fwd[0], target)
    q_rev_n = _normalize_species_label(interval_query_rev[0], target)
    t_rev_n = _normalize_species_label(interval_target_rev[0], query)

    pairs_df, fwd_ids, rev_ids = identify_reciprocal_pairs(
        q_fwd_n, t_fwd_n, q_rev_n, t_rev_n, thr=thr
    )
    print("Overlap analysis:")
    print(f"  Reciprocal block pairs: {len(pairs_df)}")

    '''
    df_fwd = (pd.read_csv(fwd_genes, sep="\t", dtype=str)
            .dropna(subset=["Region_ID","Species","Gene-ID"])
            .drop_duplicates(subset=["Region_ID","Species","Gene-ID"]))
    df_rev = (pd.read_csv(rev_genes, sep="\t", dtype=str)
                .dropna(subset=["Region_ID","Species","Gene-ID"])
                .drop_duplicates(subset=["Region_ID","Species","Gene-ID"]))

    n_query_recipr  = count_genes_in_reciprocal_ids_single(df_fwd, fwd_ids, query)
    n_target_recipr = count_genes_in_reciprocal_ids_single(df_rev,  rev_ids,  target)

    n_query_all_fwd, n_target_all_fwd = count_genes_in_blocks(fwd_genes, query, target)
    n_query_all_rev, n_target_all_rev = count_genes_in_blocks(rev_genes, target, query)

    print("Gene counts in all blocks:")
    print(f"  Forward - {query}: {n_query_all_fwd}, {target}: {n_target_all_fwd}")
    print(f"  Reverse - {target}: {n_query_all_rev}, {query}: {n_target_all_rev}\n")

    print("Gene counts in reciprocal regions:")
    print(f"  {query}: {n_query_recipr}, {target}: {n_target_recipr}\n")

    ratio_sp1, ratio_sp2 = compute_ratios_per_species(
        n_query_recipr, n_target_recipr, n_query_all_fwd, n_query_all_rev
    )
    print("Ratios:")
    print(f"  {query}: {ratio_sp1:.3f} ({ratio_sp1*100:.1f}%)")
    print(f"  {target}: {ratio_sp2:.3f} ({ratio_sp2*100:.1f}%)\n")
    '''

    # --- Compute total and reciprocal block lengths (bp-based conservativeness) ---
    # Total synteny bp in forward direction
    total_bp_query  = (interval_query_fwd[0].df["End"] - interval_query_fwd[0].df["Start"]).sum()
    total_bp_target = (interval_target_fwd[0].df["End"] - interval_target_fwd[0].df["Start"]).sum()

    # Reciprocal block bp (only those Region_IDs that appear in both directions)
    recip_bp_query  = interval_query_fwd[0].df.query("Region_ID in @fwd_ids").eval("End - Start").sum()
    recip_bp_target = interval_target_fwd[0].df.query("Region_ID in @rev_ids").eval("End - Start").sum()

    # Compute ratios (reciprocal fraction of synteny length)
    ratio_sp1 = recip_bp_query / total_bp_query if total_bp_query > 0 else 0.0
    ratio_sp2 = recip_bp_target / total_bp_target if total_bp_target > 0 else 0.0

    print("Basepair coverage (synteny length):")
    print(f"  {query}: total={total_bp_query:,}, reciprocal={recip_bp_query:,}, "
          f"ratio={ratio_sp1*100:.2f}%")
    print(f"  {target}: total={total_bp_target:,}, reciprocal={recip_bp_target:,}, "
          f"ratio={ratio_sp2*100:.2f}%\n")


    plot_conservativeness_bar(ratio_sp1, ratio_sp2, query, target, barplot, thr)

    all_metrics.append({
        "metric": METRIC_NAME,
        "pair_dir": f"{query}->{target}",
        "species": query,
        "value": ratio_sp1 * 100.0,
    })
    all_metrics.append({
        "metric": METRIC_NAME,
        "pair_dir": f"{query}->{target}",
        "species": target,
        "value": ratio_sp2 * 100.0
    })
    print(f"Plot saved to: {barplot}\n")

# =========================
# MAIN
# =========================
if __name__ == "__main__":
    for src, tgt in pairs:
        try:
            process_pairs(src, tgt, thr=THRESH)
        except Exception as e:
            print(f"\n[ERROR] {src}-{tgt}: {type(e).__name__}")
            print(f"Message: {e}")
            print("Traceback:")
            traceback.print_exc()

    expected = 2 * len(pairs)
    print(f"\nCollected {len(all_metrics)} bar values; expected {expected}.\n")

    finalize_boxplot(
        all_metrics,
        metric_name=METRIC_NAME,
        out_csv=os.path.join(OUT_DIR, "summary_points_fix.tsv"),
        out_png=os.path.join(OUT_DIR, "summary_boxplot_fix.png"),
    )
