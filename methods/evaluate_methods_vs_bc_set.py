import pandas as pd
from pathlib import Path

# Evaluate methods vs BC-set

# Input:

# Reference set: BC-set
REF_PATH = "results/evaluation/mammalia/vs_b_set/bc_set.tsv"
# Second set:
# Uncomment corresponding path
# Neighborhood - gm
SECOND_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_gm.tsv"
# Neighborhood - sum
#SECOND_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_sum.tsv"
# Tree - majority
#SECOND_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/majority/complete_tree_conserved_orthologs_2.tsv"
# Tree - standard
#SECOND_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/standard/complete_tree_conserved_orthologs_2.tsv"
# Tree - whitelist
#SECOND_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/whitelist/complete_tree_conserved_orthologs_2.tsv"
# PTP
#SECOND_PATH = "/storage/EasyVectorOmics/phylotreepruner/results/ptp_one_to_one_pairs_new.tsv"
# OrthoFinder
#SECOND_PATH = "results/evaluation/vs_orthofinder/orthofinder_1to1.tsv"

# Output:

# Uncomment corresponding path:

# neighborhood - gm
OUT_METRICS = "results/evaluation/mammalia/all_vs_bc_set/vs_neighborhood/bc_nbh_gm_metrics.txt"
OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/all_vs_bc_set/vs_neighborhood/missed_by_nbh_gm_vs_bc.tsv"
OUT_MISSED_BY_REF = "results/evaluation/mammalia/all_vs_bc_set/vs_neighborhood/missed_by_bc_vs_nbh_gm.tsv"
#neighborhood - sum
#OUT_METRICS = "results/evaluation/mammalia/all_vs_bc_set/vs_neighborhood/bc_nbh_sum_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/all_vs_bc_set/vs_neighborhood/missed_by_nbh_sum_vs_bc.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/all_vs_bc_set/vs_neighborhood/missed_by_bc_vs_nbh_sum.tsv"
# tree - majority
#OUT_METRICS = "results/evaluation/mammalia/all_vs_bc_set/vs_tree/majority/bc_majority_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/all_vs_bc_set/vs_tree/majority/missed_by_tree_majority_vs_bc.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/all_vs_bc_set/vs_tree/majority/missed_by_bc_vs_tree_majority.tsv"
# tree - standard
#OUT_METRICS = "results/evaluation/mammalia/all_vs_bc_set/vs_tree/standard/bc_standard_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/all_vs_bc_set/vs_tree/standard/missed_by_tree_standard_vs_bc.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/all_vs_bc_set/vs_tree/standard/missed_by_bc_vs_tree_standard.tsv"
# tree - whitelist
#OUT_METRICS = "results/evaluation/mammalia/all_vs_bc_set/vs_tree/whitelist/bc_whitelist_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/all_vs_bc_set/vs_tree/whitelist/missed_by_tree_whitelist_vs_bc.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/all_vs_bc_set/vs_tree/whitelist/missed_by_bc_vs_tree_whitelist.tsv"
# ptp
#OUT_METRICS = "results/evaluation/mammalia/all_vs_bc_set/vs_ptp/bc_ptp_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/all_vs_bc_set/vs_ptp/missed_by_ptp_vs_bc.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/all_vs_bc_set/vs_ptp/missed_by_bc_vs_ptp.tsv"
# OrthoFinder
#OUT_METRICS = "results/evaluation/mammalia/all_vs_bc_set/vs_orthofinder/bc_orthofinder_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/all_vs_bc_set/vs_orthofinder/missed_by_orthofinder_vs_bc.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/all_vs_bc_set/vs_orthofinder/missed_by_bc_vs_orthofinder.tsv"
# ----------------------------


def to_pair_set(df, a, b):
    """Convert two columns into a set of unordered pairs."""
    pairs = {tuple(sorted((x, y))) for x, y in zip(df[a], df[b]) if x != y}
    return pairs


def compare_sets(ref_set, second_set):
    """Return metrics dictionary comparing reference vs second pairs."""
    tp = len(ref_set & second_set)
    fn = len(ref_set - second_set)
    fp = len(second_set - ref_set)
    recall = tp / (tp + fn) if tp + fn else 0
    precision = tp / (tp + fp) if tp + fp else 0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0
    jaccard = len(ref_set & second_set) / len(ref_set | second_set)
    return dict(tp=tp, fn=fn, fp=fp, recall=recall, precision=precision, f1=f1, jaccard=jaccard)

    

def main():
    # --- Load data ---
    ref_df = pd.read_csv(REF_PATH, sep="\t")
    second_df = pd.read_csv(SECOND_PATH, sep="\t")

    # --- Convert to sets ---
    ref_set = to_pair_set(ref_df, "Protein1", "Protein2")
    # Uncomment if second file has Protein columns (for neighborhood, PTP and OrthoFinder):
    second_set = to_pair_set(second_df, "Protein1", "Protein2")
    # Uncomment if second file has Gene columns (For tree):
    #second_set = to_pair_set(second_df, "Gene1", "Gene2")
  
    # --- Identify missed pairs ---
    missed_by_second = ref_set - second_set
    missed_by_ref = second_set - ref_set

    Path(OUT_MISSED_BY_SECOND).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(list(missed_by_second), columns=["Query_Protein", "Target_Protein"]).to_csv(OUT_MISSED_BY_SECOND, sep="\t", index=False)
    pd.DataFrame(list(missed_by_ref), columns=["Query_Protein", "Target_Protein"]).to_csv(OUT_MISSED_BY_REF, sep="\t", index=False)

    # --- Compute metrics ---
    m = compare_sets(ref_set, second_set)
    text = (
        f"Pairs in ref_set: {len(ref_set)}\n"
        f"Pairs in second_set: {len(second_set)}\n\n"
        f"TP: {m['tp']}, FN: {m['fn']}, FP: {m['fp']}\n"
        f"Recall: {m['recall']:.3f}\n"
        f"Precision: {m['precision']:.3f}\n"
        f"F1 Score: {m['f1']:.3f}\n"
        f"Jaccard Similarity: {m['jaccard']:.3f}\n"
        f"Recovered {m['recall']:.1%} of ref 1:1 orthologs "
        f"({m['tp']} out of {m['tp'] + m['fn']})\n"
    )

    # --- Save metrics ---
    with open(OUT_METRICS, "w") as f:
        f.write(text)
    print(text)



if __name__ == "__main__":
    main()
