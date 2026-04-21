import pandas as pd
from pathlib import Path

# Evaluate methods vs CactBUSComp-set (union of cactus + BUSComp)
#
# Output mode (new logic):
# - Bidirectional comparison for metrics/FN/FP (A-B == B-A) using undirected keys.
# - Preserve direction for outputs:
#     * FN written in REF direction (CactBUSComp)
#     * FP written in TEST direction (method output)
# - Write ONLY combined outputs:
#       Protein1, Protein2, found_in_map
# - No species/gene/pair/split outputs (kept out to avoid mapping-driven NaNs).
#
# found_in_map=True iff BOTH Protein1 and Protein2 exist in MAP_PATH (ProteinID column).


# ----------------------------
# INPUT
# ----------------------------

# Reference set: CactBUSComp-set
REF_PATH = "results/evaluation/mammalia/vs_buscomp_set/CactBUSComp_set.tsv"

# Mapping file (only used for found_in_map flag)
MAP_PATH = "material/sex_experiment/gene_to_protein_map_FIXED.tsv"

# Test set:
# Uncomment corresponding path
# Prot-syn - n-score
#TEST_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_nscore.tsv"
# Prot-syn - sum
TEST_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_sum.tsv"
# Tree - majority
#TEST_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/majority/complete_tree_conserved_orthologs_2.tsv"
# Tree - max_score
#TEST_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/standard/complete_tree_conserved_orthologs_2.tsv"
# Tree - whitelist
#TEST_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/whitelist/complete_tree_conserved_orthologs_2.tsv"
# PTP
#TEST_PATH = "/storage/EasyVectorOmics/phylotreepruner/results/mammalia/ptp_one_to_one_pairs_new.tsv"
# OrthoFinder
#TEST_PATH = "results/evaluation/mammalia/vs_orthofinder/orthofinder_1to1.tsv"
# Output:

# Uncomment corresponding path:

# prot-syn - nscore
#OUT_METRICS = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/nscore_approach/test_no_transitivity/CactBUSComp_prot_syn_nscore_metrics.txt"
#OUT_METRICS_TSV = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/nscore_approach/test_no_transitivity/CactBUSComp_prot_syn_nscore_metrics.tsv"
#OUT_MISSED_BY_TEST = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/nscore_approach/test_no_transitivity/all_false_negatives_prot_syn_nscore.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/nscore_approach/test_no_transitivity/all_false_positives_prot_syn_nscore.tsv"
#OUT_SPLIT_DIR_FN = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/nscore_approach/test_no_transitivity"
#OUT_SPLIT_DIR_FP = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/nscore_approach/test_no_transitivity"
#OUT_UNMAPPED_FN = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/nscore_approach/test_no_transitivity/false_negatives_UNMAPPED.tsv"
#OUT_UNMAPPED_FP = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/nscore_approach/test_no_transitivity/false_positives_UNMAPPED.tsv"

# prot-syn - sum
OUT_METRICS = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/sum_approach/test_no_transitivity/CactBUSComp_prot_syn_sum_metrics.txt"
OUT_METRICS_TSV = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/sum_approach/test_no_transitivity/CactBUSComp_prot_syn_sum_metrics.tsv"
OUT_MISSED_BY_TEST = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/sum_approach/test_no_transitivity/all_false_negatives_prot_syn_sum.tsv"
OUT_MISSED_BY_REF = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/sum_approach/test_no_transitivity/all_false_positives_prot_syn_sum.tsv"
#OUT_SPLIT_DIR_FN = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/sum_approach/test_no_transitivity"
#OUT_SPLIT_DIR_FP = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/sum_approach/test_no_transitivity"
#OUT_UNMAPPED_FN = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/sum_approach/test_no_transitivity/false_negatives_UNMAPPED.tsv"
#OUT_UNMAPPED_FP = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/sum_approach/test_no_transitivity/false_positives_UNMAPPED.tsv"

# tree - majority
#OUT_METRICS = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/majority/CactBUSComp_tree_majority_metrics.txt"
#OUT_METRICS_TSV = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/majority/CactBUSComp_tree_majority_metrics.tsv"
#OUT_MISSED_BY_TEST = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/majority/all_false_negatives_tree_majority.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/majority/all_false_positives_tree_majority.tsv"
#OUT_SPLIT_DIR_FN = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/majority"
#OUT_SPLIT_DIR_FP = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/majority"
#OUT_UNMAPPED_FN = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/majority/false_negatives_UNMAPPED.tsv"
#OUT_UNMAPPED_FP = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/majority/false_positives_UNMAPPED.tsv"

# tree - max_score
#OUT_METRICS = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/max_score/CactBUSComp_tree_max_score_metrics.txt"
#OUT_METRICS_TSV = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/max_score/CactBUSComp_tree_max_score_metrics.tsv"
#OUT_MISSED_BY_TEST = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/max_score/all_false_negatives_tree_max_score.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/max_score/all_false_positives_tree_max_score.tsv"
#OUT_SPLIT_DIR_FN = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/max_score"
#OUT_SPLIT_DIR_FP = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/max_score"
#OUT_UNMAPPED_FN = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/max_score/false_negatives_UNMAPPED.tsv"
#OUT_UNMAPPED_FP = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/max_score/false_positives_UNMAPPED.tsv"

# tree - whitelist
#OUT_METRICS = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/whitelist/CactBUSComp_tree_whitelist_metrics.txt"
#OUT_METRICS_TSV = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/whitelist/CactBUSComp_tree_whitelist_metrics.tsv"
#OUT_MISSED_BY_TEST = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/whitelist/all_false_negatives_tree_whitelist.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/whitelist/all_false_positives_tree_whitelist.tsv"
#OUT_SPLIT_DIR_FN = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/whitelist"
#OUT_SPLIT_DIR_FP = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/whitelist"
#OUT_UNMAPPED_FN = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/whitelist/false_negatives_UNMAPPED.tsv"
#OUT_UNMAPPED_FP = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/whitelist/false_positives_UNMAPPED.tsv"

# ptp
#OUT_METRICS = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_ptp/CactBUSComp_ptp_metrics.txt"
#OUT_METRICS_TSV = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_ptp/CactBUSComp_ptp_metrics.tsv"
#OUT_MISSED_BY_TEST = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_ptp/all_false_negatives_ptp.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_ptp/all_false_positives_ptp.tsv"
#OUT_SPLIT_DIR_FN = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_ptp"
#OUT_SPLIT_DIR_FP = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_ptp"
#OUT_UNMAPPED_FN = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_ptp/false_negatives_UNMAPPED.tsv"
#OUT_UNMAPPED_FP = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_ptp/false_positives_UNMAPPED.tsv"

# OrthoFinder
#OUT_METRICS = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_orthofinder/CactBUSComp_orthofinder_metrics.txt"
#OUT_METRICS_TSV = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_orthofinder/CactBUSComp_orthofinder_metrics.tsv"
#OUT_MISSED_BY_TEST = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_orthofinder/all_false_negatives_orthofinder.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_orthofinder/all_false_positives_orthofinder.tsv"
#OUT_SPLIT_DIR_FN = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_orthofinder"
#OUT_SPLIT_DIR_FP = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_orthofinder"
#OUT_UNMAPPED_FN = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_orthofinder/false_negatives_UNMAPPED.tsv"
#OUT_UNMAPPED_FP = "results/evaluation/mammalia/all_vs_cactbuscomp/vs_orthofinder/false_positives_UNMAPPED.tsv"
# ----------------------------



# ----------------------------
# IO: reference loader
# ----------------------------

def load_cactbuscomp_pairs(file_path, sep="\t"):
    """
    Load CactBUSComp file which may have different column names.
    Accepts:
      - Query_Protein / Target_Protein
      - Protein1 / Protein2 (auto-renamed to Query_/Target_)
    """
    df = pd.read_csv(file_path, sep=sep, dtype=str)

    if "Query_Protein" in df.columns and "Target_Protein" in df.columns:
        return df
    if "Protein1" in df.columns and "Protein2" in df.columns:
        return df.rename(columns={"Protein1": "Query_Protein", "Protein2": "Target_Protein"})

    raise ValueError(
        "CactBUSComp file must contain either (Query_Protein, Target_Protein) "
        f"or (Protein1, Protein2). Found: {df.columns.tolist()}"
    )


# ----------------------------
# Robust pair handling (bidirectional matching)
# ----------------------------

def build_pair_index(df, qcol, tcol):
    """
    Build:
      - keys: set of undirected pair keys (tuple(sorted((q, t))))
      - directed: dict mapping undirected key -> (q, t) as it appears in df (first occurrence)

    Bidirectional matching: A-B == B-A for metrics/FN/FP.
    """
    keys = set()
    directed = {}
    flipped = 0

    for q, t in zip(df[qcol], df[tcol]):
        if pd.isna(q) or pd.isna(t) or q == t:
            continue

        key = tuple(sorted((q, t)))
        keys.add(key)

        if key in directed and directed[key] != (q, t):
            flipped += 1

        directed.setdefault(key, (q, t))

    if flipped:
        print(f"WARNING: {flipped} pair(s) appear with BOTH orientations in columns {qcol}/{tcol}. Keeping first occurrence.")
    return keys, directed


def compare_sets(ref_set, test_set):
    tp = len(ref_set & test_set)
    fn = len(ref_set - test_set)
    fp = len(test_set - ref_set)
    recall = tp / (tp + fn) if tp + fn else 0
    precision = tp / (tp + fp) if tp + fp else 0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0
    union_size = len(ref_set | test_set)
    jaccard = len(ref_set & test_set) / union_size if union_size else 0
    return dict(tp=tp, fn=fn, fp=fp, recall=recall, precision=precision, f1=f1, jaccard=jaccard)


def save_metrics_as_tsv(metrics_dict, ref_count, test_count, output_path):
    metrics_df = pd.DataFrame({
        "Num_Pairs_Reference": [ref_count],
        "Num_Pairs_Test": [test_count],
        "Num_TP": [metrics_dict["tp"]],
        "Num_FN": [metrics_dict["fn"]],
        "Num_FP": [metrics_dict["fp"]],
        "Precision": [metrics_dict["precision"]],
        "Recall": [metrics_dict["recall"]],
        "F1_score": [metrics_dict["f1"]],
        "Jaccard_Similarity": [metrics_dict["jaccard"]],
        "%_of_ref_recovered": [metrics_dict["recall"] * 100],
    })
    metrics_df.to_csv(output_path, sep="\t", index=False)
    print(f"Saved metrics TSV to: {output_path}")


# ----------------------------
# Minimal mapping flag helper
# ----------------------------

def add_found_in_map_flag(df, map_protein_ids):
    """
    Adds found_in_map=True iff BOTH proteins exist in map_protein_ids.
    Expects columns: Query_Protein, Target_Protein
    """
    df = df.copy()
    df["found_in_map"] = df["Query_Protein"].isin(map_protein_ids) & df["Target_Protein"].isin(map_protein_ids)
    return df


# ----------------------------
# MAIN
# ----------------------------

def main():
    # --- Load ref (CactBUSComp) ---
    print("Loading CactBUSComp reference set...")
    ref_df = load_cactbuscomp_pairs(REF_PATH)

    # --- Load test (methods) ---
    print("Loading test set...")
    test_df = pd.read_csv(TEST_PATH, sep="\t", dtype=str)

    # Choose Protein1/Protein2 if present, else Gene1/Gene2
    if "Protein1" in test_df.columns and "Protein2" in test_df.columns:
        test_qcol, test_tcol = "Protein1", "Protein2"
        print(f"Test file columns used for comparison: {test_qcol}, {test_tcol}")
    elif "Gene1" in test_df.columns and "Gene2" in test_df.columns:
        test_qcol, test_tcol = "Gene1", "Gene2"
        print(f"Test file columns used for comparison: {test_qcol}, {test_tcol}")
    else:
        raise ValueError(
            "Test file must contain either (Protein1, Protein2) or (Gene1, Gene2). "
            f"Found: {test_df.columns.tolist()}"
        )

    # --- Build undirected keys + keep directed orientation for outputs ---
    ref_set, ref_dir = build_pair_index(ref_df, "Query_Protein", "Target_Protein")
    test_set, test_dir = build_pair_index(test_df, test_qcol, test_tcol)

    print(f"\nPairs in ref_set (CactBUSComp): {len(ref_set)}")
    print(f"Pairs in test_set: {len(test_set)}")

    # --- FN / FP (undirected) ---
    missed_by_test = ref_set - test_set   # FN: ref has, test misses
    missed_by_ref = test_set - ref_set    # FP: test has, ref misses

    # --- Preserve direction for outputs ---
    fn_rows = [ref_dir[k] for k in missed_by_test]   # REF direction
    fp_rows = [test_dir[k] for k in missed_by_ref]   # TEST direction

    Path(OUT_MISSED_BY_TEST).parent.mkdir(parents=True, exist_ok=True)
    fn_df = pd.DataFrame(fn_rows, columns=["Query_Protein", "Target_Protein"])
    fp_df = pd.DataFrame(fp_rows, columns=["Query_Protein", "Target_Protein"])

    # --- Metrics ---
    m = compare_sets(ref_set, test_set)
    text = (
        f"Pairs in ref_set (CactBUSComp): {len(ref_set)}\n"
        f"Pairs in test_set: {len(test_set)}\n\n"
        f"TP: {m['tp']}, FN: {m['fn']}, FP: {m['fp']}\n"
        f"Recall: {m['recall']:.3f}\n"
        f"Precision: {m['precision']:.3f}\n"
        f"F1 Score: {m['f1']:.3f}\n"
        f"Jaccard Similarity: {m['jaccard']:.3f}\n"
        f"Recovered {m['recall']:.1%} of reference 1:1 orthologs "
        f"({m['tp']} out of {m['tp'] + m['fn']})\n"
    )

    Path(OUT_METRICS).parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_METRICS, "w") as f:
        f.write(text)
    print("\n" + text)

    save_metrics_as_tsv(m, len(ref_set), len(test_set), OUT_METRICS_TSV)

    # --- found_in_map flag ---
    map_df = pd.read_csv(MAP_PATH, sep="\t", dtype=str)
    if "ProteinID" not in map_df.columns:
        raise ValueError("Map file must contain column: ProteinID")

    map_protein_ids = set(map_df["ProteinID"].dropna().astype(str))

    fn_df = add_found_in_map_flag(fn_df, map_protein_ids)
    fp_df = add_found_in_map_flag(fp_df, map_protein_ids)

    # --- Final minimal outputs ---
    fn_out = fn_df.rename(columns={"Query_Protein": "Protein1", "Target_Protein": "Protein2"})[
        ["Protein1", "Protein2", "found_in_map"]
    ]
    fp_out = fp_df.rename(columns={"Query_Protein": "Protein1", "Target_Protein": "Protein2"})[
        ["Protein1", "Protein2", "found_in_map"]
    ]

    fn_out.to_csv(OUT_MISSED_BY_TEST, sep="\t", index=False)
    fp_out.to_csv(OUT_MISSED_BY_REF, sep="\t", index=False)

    print("Saved minimal FN table:", OUT_MISSED_BY_TEST)
    print("Saved minimal FP table:", OUT_MISSED_BY_REF)


if __name__ == "__main__":
    main()
