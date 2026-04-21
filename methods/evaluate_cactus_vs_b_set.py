import pandas as pd
from pathlib import Path
import pickle
#from collections import Counter 

# Evaluate cactus vs BUSCO+Compleasm (reference); Create BC-set 

# Input:

# reference set: b-set (BUSCO+compleasm)
REF_PATH = "/storage/EasyVectorOmics/ortholog_evaluation/compare_busco_compleasm/b_set.tsv"
#REF_PATH = "material/compleasm/b_set.tsv"

# second set: cactus files
SECOND_PATH = [
    "results/reciprocal_pairs/can_mac_one2one_fix.tsv",
    "results/reciprocal_pairs/can_rat_one2one_fix.tsv",
    "results/reciprocal_pairs/can_mus_one2one_fix.tsv",
    "results/reciprocal_pairs/mac_rat_one2one_fix.tsv",
    "results/reciprocal_pairs/mac_mus_one2one_fix.tsv",
    "results/reciprocal_pairs/rat_mus_one2one_fix.tsv",
]
# Map protein-species
#MAP_PATH = "updated_pipeline/material/sex_experiment/protein_gene_map.tsv"

# Output:

OUT_METRICS = "results/evaluation/vs_b_set/cactus_b_set_metrics.txt"
OUT_MISSED_BY_REF = "results/evaluation/vs_b_set/missed_by_b_set_vs_cactus_new.tsv"
OUT_MISSED_BY_SECOND = "results/evaluation/vs_b_set/missed_by_cactus_vs_b_set.tsv"
OUT_BC_SET_PKL = "results/evaluation/vs_b_set/bc_set.pkl"
OUT_BC_SET_TSV = "results/evaluation/vs_b_set/bc_set.tsv"
# disagreement removed
#OUT_BC_SET_NO_DIS_PKL = "updated_pipeline/results/evaluation/vs_b_set/bc_set_no_disagreement.pkl"
#OUT_BC_SET_NO_DIS_TSV = "updated_pipeline/results/evaluation/vs_b_set/bc_set_no_disagreement.tsv"
# ----------------------------


def load_pairs(files, query_col="Query_Protein", target_col="Target_Protein", sep="\t"):
    """
    Load multiple cactus CSVs and extract only the columns needed for the protein pairs.
    Automatically ignores other columns present in the file.
    """
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep=sep)
        if query_col not in df.columns or target_col not in df.columns:
            raise ValueError(f"Missing required columns in {f}: expected '{query_col}' and '{target_col}'")
        dfs.append(df[[query_col, target_col]])
    return pd.concat(dfs, ignore_index=True)



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

'''
def create_bc_set_with_mapping(ref_df, second_df, map_df):
    """
    Combine BUSCO+Compleasm and Cactus ortholog pairs using species mapping.
    Keeps only pairs where no protein has contradictory partners within the same species combination.
    """

    # --- Attach species info using mapping file ---
    ref_df = ref_df.merge(map_df, left_on="Protein1", right_on="ProteinID", how="left").rename(columns={"Species": "Query_Species"}).drop(columns="ProteinID")
    ref_df = ref_df.merge(map_df, left_on="Protein2", right_on="ProteinID", how="left").rename(columns={"Species": "Target_Species"}).drop(columns="ProteinID")

    second_df = second_df.merge(map_df, left_on="Query_Protein", right_on="ProteinID", how="left").rename(columns={"Species": "Query_Species"}).drop(columns="ProteinID")
    second_df = second_df.merge(map_df, left_on="Target_Protein", right_on="ProteinID", how="left").rename(columns={"Species": "Target_Species"}).drop(columns="ProteinID")

    # --- Convert to species-aware tuples ---
    ref_pairs = set(zip(ref_df["Query_Species"], ref_df["Target_Species"], ref_df["Protein1"], ref_df["Protein2"]))
    sec_pairs = set(zip(second_df["Query_Species"], second_df["Target_Species"], second_df["Query_Protein"], second_df["Target_Protein"]))

    union = ref_pairs | sec_pairs
    final_bc = set()

    # --- Build lookups for partner resolution ---
    def build_lookup(pairs):
        lookup = defaultdict(dict)
        for sp1, sp2, p1, p2 in pairs:
            lookup[(p1, sp2)] = p2
            lookup[(p2, sp1)] = p1
        return lookup

    ref_lookup = build_lookup(ref_pairs)
    sec_lookup = build_lookup(sec_pairs)

    # --- Keep only non-conflicting pairs ---
    removed_by_conflict = []
    for sp1, sp2, p1, p2 in union:
        ref_partner = ref_lookup.get((p1, sp2))
        sec_partner = sec_lookup.get((p1, sp2))
        if (ref_partner in [None, p2]) and (sec_partner in [None, p2]):
            final_bc.add((sp1, sp2, p1, p2))
        else:
            removed_by_conflict.append((sp1, sp2, p1, p2))

    # Check if removed pairs are purely mapping mismatches
    removed_df = pd.DataFrame(removed_by_conflict, columns=["Query_Species", "Target_Species", "Protein1", "Protein2"])
    n_nan = removed_df["Query_Species"].isna().sum() + removed_df["Target_Species"].isna().sum()
    print(f"{n_nan / len(removed_df) * 100:.2f}% of removed pairs have NaN species info")


    # --- Diagnostics ---
    print("\n=== BC-set Disagreement Summary ===")
    print(f"Total pairs in union: {len(union):,}")
    print(f"Pairs kept (no disagreement): {len(final_bc):,}")
    print(f"Removed due to disagreement: {len(removed_by_conflict):,}")

    # Count disagreements per species pair
    if removed_by_conflict:
        sp_counts = Counter(f"{a}-{b}" for a, b, _, _ in removed_by_conflict)
        print("\nDisagreements per species-pair:")
        for sp, n in sorted(sp_counts.items(), key=lambda x: -x[1]):
            print(f"  {sp}: {n}")

    print("===================================\n")

    return union, final_bc
'''


def create_bc_set(ref_set, second_set):
    '''
    # For cactus vs BUSCO+compleasm:
    # Create BC-set = union of both sets. 
    ''' 

    # Index partners for each gene per set
    #partner_ref = {a: b for a, b in ref_set} | {b: a for a, b in ref_set}
    #partner_sec = {a: b for a, b in second_set} | {b: a for a, b in second_set}

    union = ref_set | second_set
    #final_bc = set()

    #for a, b in union:
        # If both sets agree on these two being linked — or one of them doesn't have a partner in the other set
    #    agree = (
    #        partner_ref.get(a) == b or a not in partner_ref
    #    ) and (
    #        partner_sec.get(a) == b or a not in partner_sec
    #    )
        # keep only if no gene has a contradictory mapping
    #    if agree:
    #        final_bc.add(tuple(sorted((a, b))))

    return union       

    

def main():
    # --- Load data ---
    ref_df = pd.read_csv(REF_PATH, sep="\t")
    second_df = load_pairs(SECOND_PATH)
    #map_df = pd.read_csv(MAP_PATH, sep="\t")[["ProteinID", "Species"]].drop_duplicates()

    #map_df["ProteinID"] = map_df["ProteinID"].astype(str)
    ref_df["Protein1"] = ref_df["Protein1"].astype(str)
    ref_df["Protein2"] = ref_df["Protein2"].astype(str)
    second_df["Query_Protein"] = second_df["Query_Protein"].astype(str)
    second_df["Target_Protein"] = second_df["Target_Protein"].astype(str)


    # --- Clean the mapping file ---
    # Remove clearly invalid entries (like containing cofactors, compounds, or brackets)
    '''
    invalid_terms = [
        "ubiquinone", "NAD", "GDP", "GTP", "ADP", "NADP", "NADPH", 
        "acyl", "isomerizing", "quinone", "ammonia", "carboxylating",
        "enoyl", "Mn", "Cu", "flavin", "starch", "protein"
    ]
    pattern = "|".join(invalid_terms)
    map_df = map_df[~map_df["Species"].str.contains(pattern, case=False, na=False)]
    '''

    # Map debugging
    '''
    map_ids = set(map_df["ProteinID"])
    ref_ids = set(ref_df["Protein1"]) | set(ref_df["Protein2"])
    sec_ids = set(second_df["Query_Protein"]) | set(second_df["Target_Protein"])

    print(f"Missing in map (ref): {len(ref_ids - map_ids)}")
    print(f"Missing in map (second): {len(sec_ids - map_ids)}")
    '''

    # --- Convert to sets ---
    ref_set = to_pair_set(ref_df, "Protein1", "Protein2")
    second_set = to_pair_set(second_df, "Query_Protein", "Target_Protein")
  
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

    # --- Create and save BC-set ---
    #bc_set, bc_set_no_dis = create_bc_set_with_mapping(ref_df, second_df, map_df)
    bc_set = create_bc_set(ref_set, second_set)
    #Path(OUT_BC_SET_TSV).parent.mkdir(parents=True, exist_ok=True)

    # check if there are 1:many, appearing after union
    #genes = [g for pair in bc_set for g in pair]
    #counts = Counter(genes)
    #violations = sum(1 for v in counts.values() if v > 1)

    #print(f"Genes that appear in more than one pair: {violations}")

    # save bc-set
    with open(OUT_BC_SET_PKL, "wb") as f:  # Save as .pkl file
        pickle.dump(bc_set, f)
    #pd.DataFrame(bc_set, columns=["Query_Species", "Target_Species", "Protein1", "Protein2"]).to_csv(OUT_BC_SET_TSV, sep="\t", index=False) 
    bc_pairs_df = pd.DataFrame(bc_set, columns=["Protein1", "Protein2"])
    bc_pairs_df.to_csv(OUT_BC_SET_TSV, sep="\t", index=False) # Save as .tsv file
    
    # no disagreement
    #with open(OUT_BC_SET_NO_DIS_PKL, "wb") as f: 
    #    pickle.dump(bc_set_no_dis, f)
    #pd.DataFrame(bc_set_no_dis, columns=["Query_Species", "Target_Species", "Protein1", "Protein2"]).to_csv(OUT_BC_SET_NO_DIS_TSV, sep="\t", index=False)    #bc_pairs_no_dis_df.to_csv(OUT_BC_SET_NO_DIS_TSV, sep="\t", index=False)


if __name__ == "__main__":
    main()
