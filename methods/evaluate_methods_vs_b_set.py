import pandas as pd
from pathlib import Path


# Uncomment corresponding path:

#=========
# Mammals
#=========

# Input

# reference set: b-set (BUSCO+compleasm)
REF_FILES = "/storage/EasyVectorOmics/ortholog_evaluation/compare_busco_compleasm/b_set.tsv"


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
SECOND_PATH = "results/evaluation/mammalia/vs_orthofinder/orthofinder_1to1.tsv"
#--------------------------------------------
# Output:
# Uncomment corresponding path:
#--------------------------------------------
# neighborhood - gm
#OUT_METRICS = "results/evaluation/mammalia/vs_b_set/vs_neighborhood/b_set_nbh_gm_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_b_set/vs_neighborhood/missed_by_nbh_gm_vs_b_set.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/vs_b_set/vs_neighborhood/missed_by_b_set_vs_nbh_gm.tsv"
#neighborhood - sum
#OUT_METRICS = "results/evaluation/mammalia/vs_b_set/vs_neighborhood/b_set_nbh_sum_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_b_set/vs_neighborhood/missed_by_nbh_sum_vs_b_set.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/vs_b_set/vs_neighborhood/missed_by_b_set_vs_nbh_sum.tsv"
# tree - majority
#OUT_METRICS = "results/evaluation/mammalia/vs_b_set/vs_tree/majority/ref_majority_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_b_set/vs_tree/majority/missed_by_tree_majority_vs_ref.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/vs_b_set/vs_tree/majority/missed_by_ref_vs_tree_majority.tsv"
# tree - standard
#OUT_METRICS = "results/evaluation/mammalia/vs_b_set/vs_tree/standard/ref_standard_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_b_set/vs_tree/standard/missed_by_tree_standard_vs_ref.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/vs_b_set/vs_tree/standard/missed_by_ref_vs_tree_standard.tsv"
# tree - whitelist
#OUT_METRICS = "results/evaluation/mammalia/vs_b_set/vs_tree/whitelist/ref_whitelist_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_b_set/vs_tree/whitelist/missed_by_tree_whitelist_vs_ref.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/vs_b_set/vs_tree/whitelist/missed_by_ref_vs_tree_whitelist.tsv"
# ptp
#OUT_METRICS = "results/evaluation/mammalia/vs_b_set/vs_ptp/b_set_ptp_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_b_set/vs_ptp/missed_by_ptp_vs_b_set.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/mammalia/vs_b_set/vs_ptp/missed_by_b_set_vs_ptp.tsv"
# OrthoFinder
OUT_METRICS = "results/evaluation/mammalia/vs_b_set/vs_orthofinder/b_set_orthofinder_metrics.txt"
OUT_MISSED_BY_SECOND = "results/evaluation/mammalia/vs_b_set/vs_orthofinder/missed_by_orthofinder_vs_b_set.tsv"
OUT_MISSED_BY_REF = "results/evaluation/mammalia/vs_b_set/vs_orthofinder/missed_by_b_set_vs_orthofinder.tsv"
# -----------------------------------

#=========
# Plants
#=========

# Input

'''
REF_FILES = [
    "results/reciprocal_pairs/ar_card_one2one.tsv"
]
'''

# Neighborhood - gm (global pairs)
#SECOND_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/cardamine/global_pairs_gm.tsv"
# Neighborhood - gm (one to one)
#SECOND_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/cardamine/one_to_one_gm.tsv"
# Neighborhood - sum (global pairs)
#SECOND_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/cardamine/global_pairs_sum.tsv"
# Neighborhood - sum (one to one)
#SECOND_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/cardamine/one_to_one_sum.tsv"
# Neighborhood - pairs_original
#SECOND_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/cardamine/pairs_original.tsv"
# Neighborhood - pairs_transitive_local.tsv
#SECOND_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/cardamine/pairs_transitive_local.tsv"
# Neighborhood - pairs_unique
#SECOND_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/cardamine/pairs_unique.tsv"
# Original orthologs
#SECOND_PATH = "material/cardamine/filtered_orthologs.tsv"
######################################
# Tree - majority
#SECOND_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/majority/complete_tree_conserved_orthologs_2.tsv"
# Tree - standard
#SECOND_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/standard/complete_tree_conserved_orthologs_2.tsv"
# Tree - whitelist
#SECOND_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/whitelist/complete_tree_conserved_orthologs_2.tsv"
#########################################################
# PTP
#SECOND_PATH = "/storage/EasyVectorOmics/phylotreepruner/results/plants/ptp_one_to_one_pairs.tsv"
# OrthoFinder
#SECOND_PATH = "results/evaluation/plants/vs_orthofinder/orthofinder_1to1.tsv"
#----------------------------------------------
# Output:
# Uncomment corresponding path:
#----------------------------------------------
# neighborhood - gm (global pairs)
#OUT_METRICS = "results/evaluation/plants/vs_neighborhood/ref_nbh_gm_global_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/plants/vs_neighborhood/missed_by_nbh_gm_global.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/plants/vs_neighborhood/missed_by_ref_vs_nbh_gm_global.tsv"
# neighborhood - gm (one to one)
#OUT_METRICS = "results/evaluation/plants/vs_neighborhood/ref_nbh_gm_one_to_one_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/plants/vs_neighborhood/missed_by_nbh_gm_one_to_one.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/plants/vs_neighborhood/missed_by_ref_vs_nbh_gm_one_to_one.tsv"
# neighborhood - sum (global pairs)
#OUT_METRICS = "results/evaluation/plants/vs_neighborhood/ref_nbh_sum_global_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/plants/vs_neighborhood/missed_by_nbh_sum_global.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/plants/vs_neighborhood/missed_by_ref_vs_nbh_sum_global.tsv"
# neighborhood - sum (one to one)
#OUT_METRICS = "results/evaluation/plants/vs_neighborhood/ref_nbh_sum_one_to_one_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/plants/vs_neighborhood/missed_by_nbh_sum_one_to_one.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/plants/vs_neighborhood/missed_by_ref_vs_nbh_sum_one_to_one.tsv"
# neighborhood - pairs_original
#OUT_METRICS = "results/evaluation/plants/vs_neighborhood/ref_nbh_pairs_original_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/plants/vs_neighborhood/missed_by_nbh_pairs_original.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/plants/vs_neighborhood/missed_by_ref_vs_nbh_pairs_original.tsv"
# neighborhood - pairs_transitive_local
#OUT_METRICS = "results/evaluation/plants/vs_neighborhood/ref_nbh_pairs_transitive_local_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/plants/vs_neighborhood/missed_by_nbh_pairs_transitive_local.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/plants/vs_neighborhood/missed_by_ref_vs_nbh_pairs_transitive_local.tsv"
# neighborhood - pairs_unique
#OUT_METRICS = "results/evaluation/plants/vs_neighborhood/ref_nbh_pairs_unique_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/plants/vs_neighborhood/missed_by_nbh_pairs_unique.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/plants/vs_neighborhood/missed_by_ref_vs_nbh_pairs_unique.tsv"
#----------------------------------------------
# original orthologs
#OUT_METRICS = "results/evaluation/plants/vs_original_orth/ref_vs_original_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/plants/vs_original_orth/missed_by_original.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/plants/vs_original_orth/missed_by_ref_vs_original.tsv"
# tree - majority
#OUT_METRICS = "results/evaluation/vs_tree/majority/ref_majority_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/vs_tree/majority/missed_by_tree_majority_vs_ref.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/vs_tree/majority/missed_by_ref_vs_tree_majority.tsv"
# tree - standard
#OUT_METRICS = "results/evaluation/vs_tree/standard/ref_standard_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/vs_tree/standard/missed_by_tree_standard_vs_ref.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/vs_tree/standard/missed_by_ref_vs_tree_standard.tsv"
# tree - whitelist
#OUT_METRICS = "results/evaluation/vs_tree/whitelist/ref_whitelist_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/vs_tree/whitelist/missed_by_tree_whitelist_vs_ref.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/vs_tree/whitelist/missed_by_ref_vs_tree_whitelist.tsv"
# ptp
#OUT_METRICS = "results/evaluation/plants/vs_ptp/ref_ptp_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/plants/vs_ptp/missed_by_ptp_vs_ref.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/plants/vs_ptp/missed_by_ref_vs_ptp.tsv"
# OrthoFinder
#OUT_METRICS = "results/evaluation/plants/vs_orthofinder/ref_orthofinder_metrics_new.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/plants/vs_orthofinder/missed_by_orthofinder_vs_ref_new.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/plants/vs_orthofinder/missed_by_ref_vs_orthofinder_new.tsv"
#-----------------------------

#=============
# Drosophila
#=============

# Input

# Drosophila

'''
REF_FILES = [
    "results/reciprocal_pairs/dmel_dsec_one2one.tsv",
    "results/reciprocal_pairs/dmel_dsim_one2one.tsv",
    "results/reciprocal_pairs/dsec_dsim_one2one.tsv"
]
'''
# Neighborhood - gm (global pairs)
#SECOND_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/drosophila/global_pairs_gm.tsv"
# Neighborhood - gm (one to one)
#SECOND_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/drosophila/one_to_one_gm.tsv"
######################################
# Tree - majority
#SECOND_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/majority/complete_tree_conserved_orthologs_2.tsv"
# Tree - standard
#SECOND_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/standard/complete_tree_conserved_orthologs_2.tsv"
# Tree - whitelist
#SECOND_PATH = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/whitelist/complete_tree_conserved_orthologs_2.tsv"
#########################################################
# PTP
#SECOND_PATH = "/storage/EasyVectorOmics/phylotreepruner/results/drosophila/ptp_one_to_one_pairs.tsv"
# OrthoFinder
#SECOND_PATH = "results/evaluation/drosophila/vs_orthofinder/orthofinder_1to1.tsv"
#----------------------------------
# Output:
# Uncomment corresponding path:
#----------------------------------
# neighborhood - gm (global pairs)
#OUT_METRICS = "results/evaluation/drosophila/vs_neighborhood/ref_nbh_gm_global_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/drosophila/vs_neighborhood/missed_by_nbh_gm_global.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/drosophila/vs_neighborhood/missed_by_ref_vs_nbh_gm_global.tsv"
# neighborhood - gm (one to one)
#OUT_METRICS = "results/evaluation/drosophila/vs_neighborhood/ref_nbh_gm_one_to_one_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/drosophila/vs_neighborhood/missed_by_nbh_gm_one_to_one.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/drosophila/vs_neighborhood/missed_by_ref_vs_nbh_gm_one_to_one.tsv"
# tree - majority
#OUT_METRICS = "results/evaluation/vs_tree/majority/ref_majority_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/vs_tree/majority/missed_by_tree_majority_vs_ref.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/vs_tree/majority/missed_by_ref_vs_tree_majority.tsv"
# tree - standard
#OUT_METRICS = "results/evaluation/vs_tree/standard/ref_standard_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/vs_tree/standard/missed_by_tree_standard_vs_ref.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/vs_tree/standard/missed_by_ref_vs_tree_standard.tsv"
# tree - whitelist
#OUT_METRICS = "results/evaluation/vs_tree/whitelist/ref_whitelist_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/vs_tree/whitelist/missed_by_tree_whitelist_vs_ref.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/vs_tree/whitelist/missed_by_ref_vs_tree_whitelist.tsv"
# ptp
#OUT_METRICS = "results/evaluation/drosophila/vs_ptp/ref_ptp_metrics.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/drosophila/vs_ptp/missed_by_ptp_vs_ref.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/drosophila/vs_ptp/missed_by_ref_vs_ptp.tsv"
# OrthoFinder
#OUT_METRICS = "results/evaluation/drosophila/vs_orthofinder/ref_orthofinder_metrics_new.txt"
#OUT_MISSED_BY_SECOND = "results/evaluation/drosophila/vs_orthofinder/missed_by_orthofinder_vs_ref_new.tsv"
#OUT_MISSED_BY_REF = "results/evaluation/drosophila/vs_orthofinder/missed_by_ref_vs_orthofinder_new.tsv"
#-----------------------------

'''
def load_pairs(files, query_col="Query_Protein", target_col="Target_Protein", sep="\t"):
    """
    Load multiple ref CSVs and extract only the columns needed for the protein pairs.
    Automatically ignores other columns present in the file.
    """
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep=sep)
        if query_col not in df.columns or target_col not in df.columns:
            raise ValueError(f"Missing required columns in {f}: expected '{query_col}' and '{target_col}'")
        dfs.append(df[[query_col, target_col]])
    return pd.concat(dfs, ignore_index=True)
'''



def to_pair_set(df, a, b):
    """Convert two columns into a set of unordered pairs."""
    pairs = {tuple(sorted((x, y))) for x, y in zip(df[a], df[b]) if x != y}
    return pairs


def compare_sets(ref_set, second_set):
    """Return metrics dictionary comparing ref vs nbh pairs."""
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
    ref_df = pd.read_csv(REF_FILES, sep="\t")
    # Uncomment if second file has Protein columns (for neighborhood, PTP, and OrthoFinder):
    second_df = pd.read_csv(SECOND_PATH, sep="\t", usecols=["Protein1", "Protein2"])
    # Uncomment if second file has Gene columns (For tree and original plants orthlogs):
    #second_df = pd.read_csv(SECOND_PATH, sep="\t", usecols=["Gene1", "Gene2"])

    # --- Convert to sets ---
    ref_set = to_pair_set(ref_df, "Protein1", "Protein2")
    # Uncomment if second file has Protein columns (for neighborhood, PTP and OrthoFinder):
    second_set = to_pair_set(second_df, "Protein1", "Protein2")
    # Uncomment if second file has Gene columns (For tree and original plants orthlogs):
    #second_set = to_pair_set(second_df, "Gene1", "Gene2")

    # --- Identify missed pairs ---
    missed_by_second = ref_set - second_set
    missed_by_ref = second_set - ref_set

    Path(OUT_MISSED_BY_SECOND).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(list(missed_by_second), columns=["Protein1", "Protein2"]).to_csv(OUT_MISSED_BY_SECOND, sep="\t", index=False)
    pd.DataFrame(list(missed_by_ref), columns=["Protein1", "Protein2"]).to_csv(OUT_MISSED_BY_REF, sep="\t", index=False)

    # --- Compute metrics ---
    m = compare_sets(ref_set, second_set)
    text = (
        f"Pairs in cactus_set: {len(ref_set)}\n"
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
