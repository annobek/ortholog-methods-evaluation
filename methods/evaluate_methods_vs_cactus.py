import pandas as pd
from pathlib import Path
import re

# Uncomment corresponding path:

#=========
# Mammals
#=========

# Input


CACTUS_FILES = [
    "results/reciprocal_pairs/can_mac_one2one_fix.tsv",
    "results/reciprocal_pairs/can_rat_one2one_fix.tsv",
    "results/reciprocal_pairs/can_mus_one2one_fix.tsv",
    "results/reciprocal_pairs/mac_rat_one2one_fix.tsv",
    "results/reciprocal_pairs/mac_mus_one2one_fix.tsv",
    "results/reciprocal_pairs/rat_mus_one2one_fix.tsv",
]

MAP_PATH = "material/sex_experiment/gene_to_protein_map_FIXED.tsv"

# Prot-synt - n-score
#SECOND_PATH = "material/sex_experiment/one_to_one_nscore.tsv"
# Prot-synt - sum
SECOND_PATH = "material/sex_experiment/one_to_one_sum.tsv"
# Neighborhood - gm
#SECOND_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_gm.tsv"
# Neighborhood - sum
#SECOND_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_sum.tsv"
# Tree - majority
#SECOND_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/majority/complete_tree_conserved_orthologs_2.tsv"
# Tree - standard
#SECOND_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/standard/complete_tree_conserved_orthologs_2.tsv"
# Tree - whitelist
#SECOND_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/whitelist/complete_tree_conserved_orthologs_2.tsv"
# PTP
#SECOND_PATH = "/storage/EasyVectorOmics/phylotreepruner/results/mammalia/ptp_one_to_one_pairs_new.tsv"
# OrthoFinder
#SECOND_PATH = "results/evaluation/vs_orthofinder/orthofinder_1to1.tsv"
#--------------------------------------------
# Output:
# Uncomment corresponding path:
#--------------------------------------------
# prot-synt - nscore
#OUT_METRICS = "results/evaluation/mammalia/vs_prot_syn/nscore_approach/cactus_prot_syn_nscore_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_prot_syn/nscore_approach/all_false_negatives_nscore.tsv"
#OUT_MISSED_BY_CACTUS = "results/evaluation/mammalia/vs_prot_syn/nscore_approach/all_false_positives_nscore.tsv"

#OUT_SPLIT_DIR_FN = "results/evaluation/mammalia/vs_prot_syn/nscore_approach"
#OUT_SPLIT_DIR_FP = "results/evaluation/mammalia/vs_prot_syn/nscore_approach"
#OUT_UNMAPPED_FN = "results/evaluation/mammalia/vs_prot_syn/nscore_approach/false_negatives_UNMAPPED.tsv"
#OUT_UNMAPPED_FP = "results/evaluation/mammalia/vs_prot_syn/nscore_approach/false_positives_UNMAPPED.tsv"

# prot-synt - sum
OUT_METRICS = "results/evaluation/mammalia/vs_prot_syn/sum_approach/cactus_prot_syn_sum_metrics.txt"
OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_prot_syn/sum_approach/all_false_negatives_sum.tsv"
OUT_MISSED_BY_CACTUS = "results/evaluation/mammalia/vs_prot_syn/sum_approach/all_false_positives_sum.tsv"

OUT_SPLIT_DIR_FN = "results/evaluation/mammalia/vs_prot_syn/sum_approach"
OUT_SPLIT_DIR_FP = "results/evaluation/mammalia/vs_prot_syn/sum_approach"
OUT_UNMAPPED_FN = "results/evaluation/mammalia/vs_prot_syn/sum_approach/false_negatives_UNMAPPED.tsv"
OUT_UNMAPPED_FP = "results/evaluation/mammalia/vs_prot_syn/sum_approach/false_positives_UNMAPPED.tsv"

# neighborhood - gm
#OUT_METRICS = "results/evaluation/mammalia/vs_neighborhood/cactus_nbh_gm_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_neighborhood/missed_by_nbh_gm.tsv"
#OUT_MISSED_BY_CACTUS = "results/evaluation/mammalia/vs_neighborhood/missed_by_cactus_vs_nbh_gm.tsv"
#neighborhood - sum
#OUT_METRICS = "results/evaluation/mammalia/vs_neighborhood/cactus_nbh_sum_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_neighborhood/missed_by_nbh_sum.tsv"
#OUT_MISSED_BY_CACTUS = "results/evaluation/mammalia/vs_neighborhood/missed_by_cactus_vs_nbh_sum.tsv"
# tree - majority
#OUT_METRICS = "results/evaluation/mammalia/vs_tree/majority/cactus_majority_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_tree/majority/missed_by_tree_majority_vs_cactus.tsv"
#OUT_MISSED_BY_CACTUS = "results/evaluation/mammalia/vs_tree/majority/missed_by_cactus_vs_tree_majority.tsv"
# tree - standard
#OUT_METRICS = "results/evaluation/mammalia/vs_tree/standard/cactus_standard_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_tree/standard/missed_by_tree_standard_vs_cactus.tsv"
#OUT_MISSED_BY_CACTUS = "results/evaluation/mammalia/vs_tree/standard/missed_by_cactus_vs_tree_standard.tsv"
# tree - whitelist
#OUT_METRICS = "results/evaluation/mammalia/vs_tree/whitelist/cactus_whitelist_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_tree/whitelist/missed_by_tree_whitelist_vs_cactus.tsv"
#OUT_MISSED_BY_CACTUS = "results/evaluation/mammalia/vs_tree/whitelist/missed_by_cactus_vs_tree_whitelist.tsv"
# ptp
#OUT_METRICS = "results/evaluation/mammalia/vs_ptp/cactus_ptp_metrics_new.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_ptp/missed_by_ptp_vs_cactus_new.tsv"
#OUT_MISSED_BY_CACTUS = "results/evaluation/mammalia/vs_ptp/missed_by_cactus_vs_ptp_new.tsv"
# OrthoFinder
#OUT_METRICS = "results/evaluation/mammalia/vs_orthofinder/cactus_orthofinder_metrics_new.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_orthofinder/missed_by_orthofinder_vs_cactus_new.tsv"
#OUT_MISSED_BY_CACTUS = "results/evaluation/mammalia/vs_orthofinder/missed_by_cactus_vs_orthofinder_new.tsv"
# -----------------------------------


def load_pairs(files, query_col="Query_Protein", target_col="Target_Protein", sep="\t"):
    """
    Load multiple Cactus CSVs and extract only the columns needed for the protein pairs.
    Automatically ignores other columns present in the file.
    """
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep=sep, dtype=str)
        if query_col not in df.columns or target_col not in df.columns:
            raise ValueError(f"Missing required columns in {f}: expected '{query_col}' and '{target_col}'")
        dfs.append(df[[query_col, target_col]])
    return pd.concat(dfs, ignore_index=True)


# ----------------------------
# Robust pair handling
# ----------------------------
# Goal:
# - Compare algorithms using UNORDERED keys (so A-B == B-A for metrics).
# - Preserve ORIGINAL DIRECTION for outputs:
#     * False Negatives (missed_by_second): output in CACTUS direction
#     * False Positives (missed_by_cactus): output in SECOND direction
#
# This prevents "Query/Target swapped" artifacts caused by sorting pairs
# and then labelling them as Query/Target.
# ----------------------------

def build_pair_index(df, qcol, tcol):
    """
    Build:
      - keys: set of undirected pair keys (tuple(sorted((q, t))))
      - directed: dict mapping undirected key -> (q, t) as it appears in df (first occurrence)

    This lets us:
      - compute metrics using keys (undirected)
      - write FN/FP using directed orientation from the appropriate source
    """
    keys = set()
    directed = {}
    flipped = 0

    for q, t in zip(df[qcol], df[tcol]):
        if pd.isna(q) or pd.isna(t) or q == t:
            continue

        key = tuple(sorted((q, t)))
        keys.add(key)

        # If the same undirected key appears with opposite orientation, warn.
        if key in directed and directed[key] != (q, t):
            flipped += 1

        # Keep first occurrence so output is stable.
        directed.setdefault(key, (q, t))

    if flipped:
        print(f"WARNING: {flipped} pair(s) appear with BOTH orientations in columns {qcol}/{tcol}. Keeping first occurrence.")
    return keys, directed


def compare_sets(cactus_set, second_set):
    """Return metrics dictionary comparing cactus vs nbh pairs."""
    tp = len(cactus_set & second_set)
    fn = len(cactus_set - second_set)
    fp = len(second_set - cactus_set)
    recall = tp / (tp + fn) if tp + fn else 0
    precision = tp / (tp + fp) if tp + fp else 0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0
    # FIX: Handle empty sets
    union_size = len(cactus_set | second_set)
    jaccard = len(cactus_set & second_set) / union_size if union_size else 0
    return dict(tp=tp, fn=fn, fp=fp, recall=recall, precision=precision, f1=f1, jaccard=jaccard)


# ----------------------------
# Species / gene mapping helpers
# ----------------------------

def species_code(species_name: str) -> str:
    """Map species name (robust by first two words) to a short code."""
    if not isinstance(species_name, str) or species_name.strip() == "":
        return ""
    key = " ".join(species_name.split()[:2]).lower()
    code_map = {
        "canis lupus": "can",
        "macaca fascicularis": "mac",
        "mus musculus": "mus",
        "rattus norvegicus": "rat",
    }
    return code_map.get(key, re.sub(r"[^a-z]+", "_", key)[:8])


def expected_pair_directions_from_cactus_files(cactus_files):
    """
    Build per-pair expected direction from filenames like:
      can_mac_one2one_fix.tsv -> expected ('can' -> 'mac')
    Returns dict: frozenset({a,b}) -> (a,b)  (expected query_code, target_code)
    """
    exp = {}
    for f in cactus_files:
        stem = Path(f).stem  # e.g. "rat_mus_one2one_fix"
        parts = stem.split("_")
        if len(parts) < 2:
            continue
        left, right = parts[0], parts[1]
        exp[frozenset((left, right))] = (left, right)
    return exp


def add_map_columns(df, prot2species, prot2gene=None):
    """Add Query/Target Species (+ code) and Gene columns to df."""
    df = df.copy()
    df["Query_Species"] = df["Query_Protein"].map(prot2species)
    df["Target_Species"] = df["Target_Protein"].map(prot2species)
    df["Q_code"] = df["Query_Species"].apply(species_code)
    df["T_code"] = df["Target_Species"].apply(species_code)

    if prot2gene is not None:
        df["Query_Gene"] = df["Query_Protein"].map(prot2gene)
        df["Target_Gene"] = df["Target_Protein"].map(prot2gene)

    return df


def canonicalize_direction_per_pair(df, expected_dir):
    """
    Canonicalize direction using expected cactus pair direction PER PAIR (not global).
    Keeps original columns for auditing:
      Orig_Query_Protein, Orig_Target_Protein, Orig_Query_Species, Orig_Target_Species, ...
    Then swaps Query/Target when row is reversed vs expected_dir for that species pair.

    Adds:
      Pair = "<left>_<right>" according to expected cactus direction
      Swapped = True/False
    """
    df = df.copy()

    # Preserve original for auditing
    df["Orig_Query_Protein"] = df["Query_Protein"]
    df["Orig_Target_Protein"] = df["Target_Protein"]
    df["Orig_Query_Species"] = df["Query_Species"]
    df["Orig_Target_Species"] = df["Target_Species"]
    df["Orig_Q_code"] = df["Q_code"]
    df["Orig_T_code"] = df["T_code"]
    if "Query_Gene" in df.columns and "Target_Gene" in df.columns:
        df["Orig_Query_Gene"] = df["Query_Gene"]
        df["Orig_Target_Gene"] = df["Target_Gene"]

    swapped = []

    for i, row in df.iterrows():
        qc, tc = row["Q_code"], row["T_code"]

        # If mapping missing, we can't canonicalize reliably.
        if not qc or not tc:
            df.at[i, "Pair"] = ""
            swapped.append(False)
            continue

        key = frozenset((qc, tc))
        if key not in expected_dir:
            # Not one of the expected species pairs
            df.at[i, "Pair"] = ""
            swapped.append(False)
            continue

        exp_q, exp_t = expected_dir[key]
        df.at[i, "Pair"] = f"{exp_q}_{exp_t}"

        # If current is reversed vs expected, swap.
        need_swap = (qc == exp_t and tc == exp_q)
        swapped.append(bool(need_swap))

        if need_swap:
            # swap proteins
            qprot, tprot = row["Query_Protein"], row["Target_Protein"]
            df.at[i, "Query_Protein"] = tprot
            df.at[i, "Target_Protein"] = qprot

            # swap species
            qsp, tsp = row["Query_Species"], row["Target_Species"]
            df.at[i, "Query_Species"] = tsp
            df.at[i, "Target_Species"] = qsp

            # swap codes
            df.at[i, "Q_code"] = exp_q
            df.at[i, "T_code"] = exp_t

            # swap genes (if present)
            if "Query_Gene" in df.columns and "Target_Gene" in df.columns:
                qg, tg = row["Query_Gene"], row["Target_Gene"]
                df.at[i, "Query_Gene"] = tg
                df.at[i, "Target_Gene"] = qg

    df["Swapped"] = swapped
    return df


def write_split_outputs_fullinfo(df, out_dir, unmapped_path, prefix):
    """
    Write:
      - per species-pair files (prefix_<pair>.tsv)
      - unmapped rows
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    unmapped = df[df["Pair"] == ""].copy()
    if len(unmapped) > 0:
        unmapped.to_csv(unmapped_path, sep="\t", index=False)
        print("Saved unmapped rows to:", unmapped_path)

    mappable = df[df["Pair"] != ""].copy()
    for pair_name, sub in mappable.groupby("Pair"):
        out_path = str(Path(out_dir) / f"{prefix}_{pair_name}.tsv")

        cols = [
            "Query_Protein", "Target_Protein",
            "Query_Species", "Target_Species",
        ]
        if "Query_Gene" in sub.columns and "Target_Gene" in sub.columns:
            cols += ["Query_Gene", "Target_Gene"]

        # include auditing columns
        cols += [
            "Orig_Query_Protein", "Orig_Target_Protein",
            "Orig_Query_Species", "Orig_Target_Species",
            "Swapped"
        ]
        if "Orig_Query_Gene" in sub.columns and "Orig_Target_Gene" in sub.columns:
            cols += ["Orig_Query_Gene", "Orig_Target_Gene"]

        sub[cols].to_csv(out_path, sep="\t", index=False)
        print("Saved:", out_path, "n=", len(sub))

def write_split_outputs(df, out_dir, unmapped_path, prefix):
    """
    Write:
      - per species-pair files (prefix_<pair>.tsv)
      - unmapped rows

    Output columns are kept minimal:
      Query_Protein, Target_Protein, Query_Species, Target_Species, Query_Gene, Target_Gene
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    unmapped = df[df["Pair"] == ""].copy()
    if len(unmapped) > 0:
        unmapped.to_csv(unmapped_path, sep="\t", index=False)
        print("Saved unmapped rows to:", unmapped_path)

    mappable = df[df["Pair"] != ""].copy()
    for pair_name, sub in mappable.groupby("Pair"):
        out_path = str(Path(out_dir) / f"{prefix}_{pair_name}.tsv")

        cols = [
            "Query_Protein", "Target_Protein",
            "Query_Species", "Target_Species",
            "Query_Gene", "Target_Gene"
        ]

        # If Gene columns don't exist (e.g. map lacks GeneID), fall back gracefully
        cols_existing = [c for c in cols if c in sub.columns]

        sub[cols_existing].to_csv(out_path, sep="\t", index=False)
        print("Saved:", out_path, "n=", len(sub))


def validate_canonical_direction(df, expected_dir, label):
    """
    Ensures every row with a known Pair is oriented exactly as expected.
    Raises AssertionError if any violations exist.
    """
    bad_rows = []

    for i, row in df.iterrows():
        if row.get("Pair", "") == "":
            continue  # unmapped / unknown pair
        qc, tc = row["Q_code"], row["T_code"]
        exp_q, exp_t = expected_dir[frozenset((qc, tc))]
        if (qc, tc) != (exp_q, exp_t):
            bad_rows.append((i, row["Query_Protein"], row["Target_Protein"], qc, tc, exp_q, exp_t))

    if bad_rows:
        print(f"\n[{label}] ERROR: {len(bad_rows)} row(s) not in canonical direction.")
        print("First 10:")
        for r in bad_rows[:10]:
            print("  idx=", r[0], "Q=", r[1], "T=", r[2], "got=", (r[3], r[4]), "expected=", (r[5], r[6]))
        raise AssertionError(f"{label}: canonical direction validation failed")

    print(f"[{label}] OK: all mappable rows are in canonical direction.")


def main():
    # --- Load data ---
    cactus_df = load_pairs(CACTUS_FILES)
    # Uncomment if second file has Protein columns (for neighborhood, PTP, and OrthoFinder):
    second_df = pd.read_csv(SECOND_PATH, sep="\t", usecols=["Protein1", "Protein2"], dtype=str)
    # Uncomment if second file has Gene columns (For tree and original plants orthlogs):
    # second_df = pd.read_csv(SECOND_PATH, sep="\t", usecols=["Gene1", "Gene2"], dtype=str)

    # --- Convert to undirected keys, and also keep directed orientation per source ---
    cactus_set, cactus_dir = build_pair_index(cactus_df, "Query_Protein", "Target_Protein")
    second_set, second_dir = build_pair_index(second_df, "Protein1", "Protein2")
    # For Gene columns: second_set, second_dir = build_pair_index(second_df, "Gene1", "Gene2")

    # --- Identify missed pairs (by undirected key) ---
    missed_by_second = cactus_set - second_set  # False Negatives (cactus pairs second missed)
    missed_by_cactus = second_set - cactus_set  # False Positives (second pairs not in cactus)

    # --- Save FN/FP with preserved direction ---
    # FN should be written in CACTUS direction (Query_Protein -> Target_Protein as cactus reports it)
    fn_rows = [cactus_dir[k] for k in missed_by_second]
    # FP should be written in SECOND direction (Protein1 -> Protein2 as second reports it)
    fp_rows = [second_dir[k] for k in missed_by_cactus]

    Path(OUT_MISSED_BY_SECOND).parent.mkdir(parents=True, exist_ok=True)
    fn_df = pd.DataFrame(fn_rows, columns=["Query_Protein", "Target_Protein"])
    fp_df = pd.DataFrame(fp_rows, columns=["Query_Protein", "Target_Protein"])

    fn_df.to_csv(OUT_MISSED_BY_SECOND, sep="\t", index=False)
    fp_df.to_csv(OUT_MISSED_BY_CACTUS, sep="\t", index=False)

    # --- Compute metrics ---
    m = compare_sets(cactus_set, second_set)
    text = (
        f"Pairs in cactus_set: {len(cactus_set)}\n"
        f"Pairs in second_set: {len(second_set)}\n\n"
        f"TP: {m['tp']}, FN: {m['fn']}, FP: {m['fp']}\n"
        f"Recall: {m['recall']:.3f}\n"
        f"Precision: {m['precision']:.3f}\n"
        f"F1 Score: {m['f1']:.3f}\n"
        f"Jaccard Similarity: {m['jaccard']:.3f}\n"
        f"Recovered {m['recall']:.1%} of cactus 1:1 orthologs "
        f"({m['tp']} out of {m['tp'] + m['fn']})\n"
    )

    # --- Save metrics ---
    with open(OUT_METRICS, "w") as f:
        f.write(text)
    print(text)

    # ------------------------------------------------------------
    # NEW: annotate FN/FP with species (+ gene) and split by pair
    # ------------------------------------------------------------

    map_df = pd.read_csv(MAP_PATH, sep="\t", dtype=str)
    if "ProteinID" not in map_df.columns or "Species" not in map_df.columns:
        raise ValueError("Map file must contain columns: ProteinID, Species")

    prot2species = dict(zip(map_df["ProteinID"].astype(str), map_df["Species"].astype(str)))

    prot2gene = None
    if "GeneID" in map_df.columns:
        prot2gene = (
            map_df.dropna(subset=["ProteinID"])
                  .drop_duplicates(subset=["ProteinID"], keep="first")
                  .set_index("ProteinID")["GeneID"]
                  .to_dict()
        )

    expected_dir = expected_pair_directions_from_cactus_files(CACTUS_FILES)
    # print expected directions once and verify
    print("Expected directions (from CACTUS filenames):")
    for k, (q, t) in sorted(expected_dir.items(), key=lambda x: (x[1][0], x[1][1])):
        print(f"  {q} -> {t}")


    # Add mapping columns
    fn_df = add_map_columns(fn_df, prot2species, prot2gene)
    fp_df = add_map_columns(fp_df, prot2species, prot2gene)

    # Canonicalize direction PER species pair using cactus file names (robust)
    fn_df = canonicalize_direction_per_pair(fn_df, expected_dir)
    fp_df = canonicalize_direction_per_pair(fp_df, expected_dir)

    # assert every mappable row matches the expected direction
    validate_canonical_direction(fn_df, expected_dir, "FN")
    validate_canonical_direction(fp_df, expected_dir, "FP")


    # Split and write per-pair files
    write_split_outputs(fn_df, OUT_SPLIT_DIR_FN, OUT_UNMAPPED_FN, prefix="false_negatives")
    write_split_outputs(fp_df, OUT_SPLIT_DIR_FP, OUT_UNMAPPED_FP, prefix="false_positives")


if __name__ == "__main__":
    main()