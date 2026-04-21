import pandas as pd
from pathlib import Path

# ----------------------------
# CONFIG (keep yours as-is)
# ----------------------------

REF_FILES = [
    "results/reciprocal_pairs/can_mac_one2one_fix.tsv",
    "results/reciprocal_pairs/can_rat_one2one_fix.tsv",
    "results/reciprocal_pairs/can_mus_one2one_fix.tsv",
    "results/reciprocal_pairs/mac_rat_one2one_fix.tsv",
    "results/reciprocal_pairs/mac_mus_one2one_fix.tsv",
    "results/reciprocal_pairs/rat_mus_one2one_fix.tsv",
]

MAP_PATH = "material/sex_experiment/gene_to_protein_map_FIXED.tsv"

TEST_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_nscore.tsv"

OUT_METRICS = "results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/simple_cactus_prot_syn_nscore_metrics.txt"
OUT_METRICS_TSV = "results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/simple_cactus_prot_syn_nscore_metrics.tsv"
OUT_MISSED_BY_TEST = "results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/simple_all_false_negatives_nscore.tsv"
OUT_MISSED_BY_REF = "results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/simple_all_false_positives_nscore.tsv"


# ----------------------------
# IO helpers
# ----------------------------

def load_ref_pairs(files, query_col="Query_Protein", target_col="Target_Protein", sep="\t"):
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep=sep, dtype=str)
        if query_col not in df.columns or target_col not in df.columns:
            raise ValueError(f"Missing required columns in {f}: expected '{query_col}' and '{target_col}'")
        dfs.append(df[[query_col, target_col]].rename(columns={query_col: "Protein1", target_col: "Protein2"}))
    return pd.concat(dfs, ignore_index=True)


def load_test_pairs(path, sep="\t"):
    # test uses Protein1/Protein2
    return pd.read_csv(path, sep=sep, usecols=["Protein1", "Protein2"], dtype=str)


# ----------------------------
# Core: undirected matching (bidirectional)
# ----------------------------

def build_undirected_index(df, col1="Protein1", col2="Protein2"):
    """
    Returns:
      keys: set of undirected keys (tuple(sorted((p1,p2))))
      first_dir: dict undirected key -> (p1,p2) as first seen in df
    """
    keys = set()
    first_dir = {}

    for p1, p2 in zip(df[col1], df[col2]):
        if pd.isna(p1) or pd.isna(p2) or p1 == p2:
            continue
        key = tuple(sorted((p1, p2)))
        keys.add(key)
        first_dir.setdefault(key, (p1, p2))

    return keys, first_dir


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
# Annotation: species (+ gene)
# ----------------------------

def annotate_pairs(df, prot2species, prot2gene=None):
    df = df.copy()
    df["Species1"] = df["Protein1"].map(prot2species)
    df["Species2"] = df["Protein2"].map(prot2species)

    if prot2gene is not None:
        df["Gene1"] = df["Protein1"].map(prot2gene)
        df["Gene2"] = df["Protein2"].map(prot2gene)

    # Put a nice, stable column order
    cols = ["Protein1", "Protein2", "Gene1", "Gene2", "Species1", "Species2"]
    cols = [c for c in cols if c in df.columns]
    return df[cols]


def main():
    # --- Load data ---
    ref_df = load_ref_pairs(REF_FILES)
    test_df = load_test_pairs(TEST_PATH)

    # --- Build undirected sets (bidirectional matching) ---
    ref_set, ref_first_dir = build_undirected_index(ref_df, "Protein1", "Protein2")
    test_set, test_first_dir = build_undirected_index(test_df, "Protein1", "Protein2")

    # --- FN/FP by undirected key ---
    missed_by_test = ref_set - test_set  # FN
    missed_by_ref = test_set - ref_set   # FP

    fn_rows = [ref_first_dir[k] for k in missed_by_test]   # keep REF orientation
    fp_rows = [test_first_dir[k] for k in missed_by_ref]   # keep TEST orientation

    fn_df = pd.DataFrame(fn_rows, columns=["Protein1", "Protein2"])
    fp_df = pd.DataFrame(fp_rows, columns=["Protein1", "Protein2"])

    # --- Metrics ---
    m = compare_sets(ref_set, test_set)
    text = (
        f"Pairs in ref_set: {len(ref_set)}\n"
        f"Pairs in test_set: {len(test_set)}\n\n"
        f"TP: {m['tp']}, FN: {m['fn']}, FP: {m['fp']}\n"
        f"Recall: {m['recall']:.3f}\n"
        f"Precision: {m['precision']:.3f}\n"
        f"F1 Score: {m['f1']:.3f}\n"
        f"Jaccard Similarity: {m['jaccard']:.3f}\n"
        f"Recovered {m['recall']:.1%} of ref 1:1 orthologs "
        f"({m['tp']} out of {m['tp'] + m['fn']})\n"
    )

    # --- Save outputs (mkdirs) ---
    Path(OUT_MISSED_BY_TEST).parent.mkdir(parents=True, exist_ok=True)
    Path(OUT_MISSED_BY_REF).parent.mkdir(parents=True, exist_ok=True)
    Path(OUT_METRICS).parent.mkdir(parents=True, exist_ok=True)
    Path(OUT_METRICS_TSV).parent.mkdir(parents=True, exist_ok=True)

    with open(OUT_METRICS, "w") as f:
        f.write(text)
    print(text)

    save_metrics_as_tsv(m, len(ref_set), len(test_set), OUT_METRICS_TSV)

    # --- Load mapping and annotate FN/FP with species (+ gene if available) ---
    map_df = pd.read_csv(MAP_PATH, sep="\t", dtype=str)
    if "ProteinID" not in map_df.columns or "Species" not in map_df.columns:
        raise ValueError("Map file must contain columns: ProteinID, Species")

    prot2species = (
        map_df.dropna(subset=["ProteinID"])
              .drop_duplicates(subset=["ProteinID"], keep="first")
              .set_index("ProteinID")["Species"]
              .to_dict()
    )

    prot2gene = None
    if "GeneID" in map_df.columns:
        prot2gene = (
            map_df.dropna(subset=["ProteinID"])
                  .drop_duplicates(subset=["ProteinID"], keep="first")
                  .set_index("ProteinID")["GeneID"]
                  .to_dict()
        )

    fn_df = annotate_pairs(fn_df, prot2species, prot2gene)
    fp_df = annotate_pairs(fp_df, prot2species, prot2gene)

    fn_df.to_csv(OUT_MISSED_BY_TEST, sep="\t", index=False)
    fp_df.to_csv(OUT_MISSED_BY_REF, sep="\t", index=False)

    print("Saved FN table:", OUT_MISSED_BY_TEST)
    print("Saved FP table:", OUT_MISSED_BY_REF)


if __name__ == "__main__":
    main()
