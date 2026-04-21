import pandas as pd

def check_one_to_one(path):
    df = pd.read_csv(path, sep="\t", dtype=str)

    print(f"\n=== Checking {path} ===")
    print(f"Total rows: {len(df):,}")

    # Try to auto-detect column names
    if {"Query_Gene", "Target_Gene"}.issubset(df.columns):
        q_col, t_col = "Query_Gene", "Target_Gene"
    elif {"Protein1", "Protein2"}.issubset(df.columns):
        q_col, t_col = "Protein1", "Protein2"
    else:
        raise ValueError(f"Cannot find gene columns in {path}. Columns: {list(df.columns)}")

    q_dupes = df[q_col].duplicated().sum()
    t_dupes = df[t_col].duplicated().sum()

    if q_dupes == 0 and t_dupes == 0:
        print("All pairs are strictly 1:1 reciprocal orthologs.")
    else:
        print("Not strictly 1:1!")
        print(f"  Duplicate {q_col}: {q_dupes}")
        print(f"  Duplicate {t_col}: {t_dupes}")

        # Optional deeper check
        q_counts = df[q_col].value_counts()
        t_counts = df[t_col].value_counts()
        print("\nGenes with >1 mapping:")
        print("  Queries >1:\n", q_counts[q_counts > 1].head())
        print("  Targets >1:\n", t_counts[t_counts > 1].head())


def check_one_to_one_grouped(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    print(f"\n=== Checking merged file {path} ===")
    print(f"Total rows: {len(df):,}")
    
    for (sp1, sp2), sub in df.groupby(["Species1", "Species2"]):
        q_dupes = sub["Protein1"].duplicated().sum()
        t_dupes = sub["Protein2"].duplicated().sum()
        if q_dupes == 0 and t_dupes == 0:
            status = " strictly 1:1"
        else:
            status = f" {q_dupes} / {t_dupes} duplicates"
        print(f"{sp1} <-> {sp2}: {len(sub):,} pairs -> {status}")


check_one_to_one_grouped("results/evaluation/vs_orthofinder/orthofinder_1to1.tsv")
