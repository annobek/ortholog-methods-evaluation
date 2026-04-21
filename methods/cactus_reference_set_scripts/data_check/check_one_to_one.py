import pandas as pd


def check_one_to_one(path):

    df = pd.read_csv(path, sep="\t", dtype=str)

    # --- Basic info ---
    print(f"\n=== Checking {path} ===")
    print(f"Total rows: {len(df):,}")

    # --- 1:1 test ---
    q_dupes = df["Query_Protein"].duplicated().sum()
    t_dupes = df["Target_Protein"].duplicated().sum()
    #q_dupes = df["Protein1"].duplicated().sum()
    #t_dupes = df["Protein2"].duplicated().sum()
   
    if q_dupes == 0 and t_dupes == 0:
        print("All pairs are strictly 1:1 reciprocal orthologs.")
    else:
        print("Not strictly 1:1!")
        print(f"  Duplicate Query_Genes: {q_dupes}")
        print(f"  Duplicate Target_Genes: {t_dupes}")

    # --- Optional deeper check (counts per gene) ---
    if q_dupes or t_dupes:
        print("\nGenes with >1 mapping:")
        q_counts = df["Query_Gene"].value_counts()
        t_counts = df["Target_Gene"].value_counts()
        #q_counts = df["Protein1"].value_counts()
        #t_counts = df["Protein2"].value_counts()
        
        print("  Queries >1:\n", q_counts[q_counts > 1].head())
        print("  Targets >1:\n", t_counts[t_counts > 1].head())


# Path to 1:1 file

# For mammals 

path_list = [
    "results/reciprocal_pairs/can_mac_one2one_fix.tsv",
    "results/reciprocal_pairs/can_rat_one2one_fix.tsv",
    "results/reciprocal_pairs/can_mus_one2one_fix.tsv",
    "results/reciprocal_pairs/mac_rat_one2one_fix.tsv",
    "results/reciprocal_pairs/mac_mus_one2one_fix.tsv",
    "results/reciprocal_pairs/rat_mus_one2one_fix.tsv"
]


# For plants
'''
path_list = [
    "results/reciprocal_pairs/ar_card_one2one.tsv"
]


# For drosophila
path_list = [
    "results/reciprocal_pairs/dmel_dsec_one2one.tsv",
    "results/reciprocal_pairs/dmel_dsim_one2one.tsv",
    "results/reciprocal_pairs/dsec_dsim_one2one.tsv"
]
'''




for path in path_list:
    check_one_to_one(path)