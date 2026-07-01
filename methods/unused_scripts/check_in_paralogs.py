import pandas as pd
#from pathlib import Path

def check_inparalogs(one2one_path, inparalogs_path):
    print(f"\n=== Checking {inparalogs_path} ===")

    # Load both files
    one2one = pd.read_csv(one2one_path, sep="\t", dtype=str)
    inpars = pd.read_csv(inparalogs_path, sep="\t", dtype=str)

    print(f"1:1 pairs: {len(one2one):,}")
    print(f"In-paralog pairs: {len(inpars):,}")

    # Build ortholog gene sets
    one2one_genes = set(one2one["Query_Gene"]) | set(one2one["Target_Gene"])

    # Count how many in-paralog genes are also present in 1:1 set
    overlap_q = inpars["Query_Gene"].isin(one2one_genes).sum()
    overlap_t = inpars["Target_Gene"].isin(one2one_genes).sum()

    print(f"Query genes overlapping 1:1 orthologs: {overlap_q:,}")
    print(f"Target genes overlapping 1:1 orthologs: {overlap_t:,}")

    # Each in-paralog pair should have at least one side in the ortholog set
    both_missing = (~inpars["Query_Gene"].isin(one2one_genes)) & (~inpars["Target_Gene"].isin(one2one_genes))
    print(f"In-paralog pairs with neither gene in 1:1 set: {both_missing.sum():,}")

    # Check duplication frequencies (expected high)
    q_counts = inpars["Query_Gene"].value_counts()
    t_counts = inpars["Target_Gene"].value_counts()
    print(f"Unique Query genes in in-paralogs: {len(q_counts):,}")
    print(f"Unique Target genes in in-paralogs: {len(t_counts):,}")
    print(f"Queries with >1 pairing: {(q_counts > 1).sum():,}")
    print(f"Targets with >1 pairing: {(t_counts > 1).sum():,}")

    # Optional: identify "anchor vs duplicate" roles
    inpars["role"] = inpars.apply(
        lambda r: (
            "target_inparalog" if r["Query_Gene"] in one2one_genes else
            "query_inparalog" if r["Target_Gene"] in one2one_genes else
            "ambiguous"
        ),
        axis=1
    )

    print(inpars["role"].value_counts())
    return inpars


# Example usage:

# For mammals
'''
pairs = [
    ("can", "mac"),
    ("can", "rat"),
    ("can", "mus"),
    ("mac", "rat"),
    ("mac", "mus"),
    ("rat", "mus")
]
'''

# For plants
'''
pairs = [
    ("ar", "card")
]
'''

# For drosophila

pairs = [
    ("dmel", "dsec"),
    ("dmel", "dsim"),
    ("dsec", "dsim")
]


for sp1, sp2 in pairs:
    one2one_path = f"results/reciprocal_pairs/{sp1}_{sp2}_one2one.tsv"
    inpars_path  = f"results/reciprocal_pairs/{sp1}_{sp2}_inparalogs.tsv"
    check_inparalogs(one2one_path, inpars_path)
