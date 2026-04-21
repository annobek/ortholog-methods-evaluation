import pandas as pd
#import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib_venn import venn3



ANNOTATIONS = {
    "can": "results/ann_genes/can_id_biotype_fix.tsv",
    "mac": "results/ann_genes/mac_id_biotype_fix.tsv",
    "rat": "results/ann_genes/rat_id_biotype_fix.tsv",
    "mus": "results/ann_genes/mus_id_biotype_fix.tsv",
}

MAP_FILE     = "material/sex_experiment/map_gen_protein_species.tsv"
DIAMOND_FILE = "material/sex_experiment/geometric_mean_sex.txt"

GENE_TABLE_DIR = Path("results/genes_in_blocks/final_output")
ORTHOLOG_DIR   = Path("results/reciprocal_pairs")
OUT_DIR        = Path("results/classification")
OUT_DIR.mkdir(parents=True, exist_ok=True)

SPECIES_PAIRS = [
    ("can", "mac"),
    ("can", "rat"),
    ("can", "mus"),
    ("mac", "rat"),
    ("mac", "mus"),
    ("rat", "mus"),
]



def load_annotation(path):
    df = pd.read_csv(path, sep="\t", dtype=str, header=None, names=["GeneID", "biotype"])
    if not {"GeneID", "biotype"}.issubset(df.columns):
        raise ValueError(f"{path} must have columns 'GeneID' and 'biotype'")
    return dict(zip(df["GeneID"], df["biotype"]))


def load_protein_map(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    if not {"ProteinID", "GeneID"}.issubset(df.columns):
        raise ValueError(f"Missing columns in {path}")
    return dict(zip(df["ProteinID"], df["GeneID"]))


def load_diamond_genes(diamond_file, prot2gene):
    df = pd.read_csv(diamond_file, sep=",", usecols=["qseqid", "sseqid"], engine="python")
    df["Query_Gene"]  = df["qseqid"].map(prot2gene)
    df["Target_Gene"] = df["sseqid"].map(prot2gene)
    return set(df["Query_Gene"].dropna()) | set(df["Target_Gene"].dropna())


def load_syntenic_genes(sp1, sp2):
    p = GENE_TABLE_DIR / f"{sp1}_{sp2}_synteny_tables_genes_fix.tsv"
    if not p.exists():
        return set()
    df = pd.read_csv(p, sep="\t", dtype=str)
    if "Gene-ID" not in df.columns:
        raise ValueError(f"{p} missing Gene-ID column")
    return set(df["Gene-ID"].dropna().unique())


def load_ortholog_genes(sp1, sp2):
    p = ORTHOLOG_DIR / f"{sp1}_{sp2}_one2one_fix.tsv"
    if not p.exists():
        return set()
    df = pd.read_csv(p, sep="\t", dtype=str)
    if not {"Query_Gene", "Target_Gene"}.issubset(df.columns):
        raise ValueError(f"{p} missing required columns")
    return set(df["Query_Gene"].dropna()) | set(df["Target_Gene"].dropna())


def load_inparalog_genes(sp1, sp2):
    p = ORTHOLOG_DIR / f"{sp1}_{sp2}_inparalogs_fix.tsv"
    if not p.exists():
        return set()
    df = pd.read_csv(p, sep="\t", dtype=str)
    return set(df["Query_Gene"].dropna()) | set(df["Target_Gene"].dropna())


def build_gene_to_species():
    g2s = {}
    for sp1, sp2 in SPECIES_PAIRS:
        for p in [f"{sp1}_{sp2}", f"{sp2}_{sp1}"]:
            f = GENE_TABLE_DIR / f"{p}_synteny_tables_genes_fix.tsv"
            if f.exists():
                df = pd.read_csv(f, sep="\t", dtype=str)
                if {"Species", "Gene-ID"}.issubset(df.columns):
                    sub = df[["Species", "Gene-ID"]].dropna().drop_duplicates()
                    g2s.update(dict(zip(sub["Gene-ID"], sub["Species"])))
    print(f"Built gene->species map for {len(g2s):,} genes.")
    return g2s



def classify_for_pair(sp1, sp2, gene2species, diamond_genes):
    print(f"\nProcessing {sp1}-{sp2}")

    biotype1 = load_annotation(ANNOTATIONS[sp1])
    biotype2 = load_annotation(ANNOTATIONS[sp2])

    # Load sets
    syntenic = load_syntenic_genes(sp1, sp2)
    orthologs = load_ortholog_genes(sp1, sp2)
    inparalogs = load_inparalog_genes(sp1, sp2)

    # FIX: remove ortholog anchors from in-paralogs
    inparalogs = {g for g in inparalogs if g not in orthologs}

    # Derived sets
    non_orth_syntenic = syntenic - orthologs
    homologs_not_syntenic = diamond_genes - syntenic

    # Classification
    records = []

    # Orthologs
    for g in orthologs:
        records.append((g, "ortholog"))

    # In-paralogs (syntenic protein-coding)
    for g in inparalogs:
        records.append((g, "in_paralog"))

    # Syntenic but non-coding (non-orthologous, not in-paralog)
    remaining = non_orth_syntenic - inparalogs
    for g in remaining:
        bt = biotype1.get(g) or biotype2.get(g) or "unknown"
        cls = "syntenic_non_coding" if bt != "protein_coding" else "syntenic_protein_paralog"
        records.append((g, cls))

    # Out-paralogs (homologous but non-syntenic)
    for g in homologs_not_syntenic:
        bt = biotype1.get(g) or biotype2.get(g) or "unknown"
        cls = "out_paralog" if bt == "protein_coding" else "non_coding_excluded"
        records.append((g, cls))

    df = pd.DataFrame(records, columns=["GeneID", "classification"])
    df["Species"] = df["GeneID"].map(gene2species).fillna("unknown")

    # Save outputs 
    pair_dir = OUT_DIR / f"{sp1}_{sp2}"
    pair_dir.mkdir(parents=True, exist_ok=True)

    df.to_csv(pair_dir / f"{sp1}_{sp2}_gene_classifications_fix.tsv", sep="\t", index=False)

    # Summaries
    summary = df.groupby("classification")["GeneID"].count().rename("count").reset_index()
    print(summary.to_string(index=False))

    # Sanity checks
    check_sets = {
        "Orthologs": orthologs,
        "In-paralogs": inparalogs,
        "Out-paralogs": {g for g, c in zip(df["GeneID"], df["classification"]) if c == "out_paralog"},
    }

    # Check for overlaps (they should be disjoint)
    overlaps = []
    keys = list(check_sets.keys())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            inter = check_sets[keys[i]] & check_sets[keys[j]]
            if inter:
                overlaps.append((keys[i], keys[j], len(inter)))
    if overlaps:
        print("\nOverlaps detected:")
        for a, b, n in overlaps:
            print(f"  {a} intersect {b} = {n}")
    else:
        print("No overlaps between ortholog/in-/out-paralog sets.")

    # Venn diagram 
    venn_out = pair_dir / f"{sp1}_{sp2}_venn.png"
    plt.figure(figsize=(5,5))
    venn3(subsets=(
        len(orthologs - inparalogs - homologs_not_syntenic),  # only orthologs
        len(inparalogs - orthologs - homologs_not_syntenic),  # only in-paralogs
        len(orthologs & inparalogs),                          # overlap (should be 0)
        len(homologs_not_syntenic - orthologs - inparalogs),  # only out-paralogs
        len(orthologs & homologs_not_syntenic),               # overlap (should be 0)
        len(inparalogs & homologs_not_syntenic),              # overlap (should be 0)
        0                                                     # triple overlap
    ), set_labels=("Orthologs", "In-paralogs", "Out-paralogs"))
    plt.title(f"{sp1}-{sp2} gene relationship classes")
    plt.tight_layout()
    plt.savefig(venn_out)
    plt.close()
    print(f"Venn diagram saved: {venn_out}")

    return df




def main():
    print("Species-aware orthology/paralogy classification\n")

    prot2gene = load_protein_map(MAP_FILE)
    g2s = build_gene_to_species()
    diamond_genes = load_diamond_genes(DIAMOND_FILE, prot2gene)

    all_results = []
    for sp1, sp2 in SPECIES_PAIRS:
        df = classify_for_pair(sp1, sp2, g2s, diamond_genes)
        all_results.append(df)

    merged = pd.concat(all_results, ignore_index=True)

    # Global summaries
    summary = (
        merged.groupby("classification")["GeneID"]
        .count()
        .rename("count")
        .reset_index()
        .sort_values("count", ascending=False)
    )
    summary_by_species = (
        merged.groupby(["Species", "classification"])["GeneID"]
        .count()
        .rename("count")
        .reset_index()
        .sort_values(["Species", "classification"])
    )

    summary.to_csv(OUT_DIR / "global_summary_fix.tsv", sep="\t", index=False)
    summary_by_species.to_csv(OUT_DIR / "global_summary_by_species_fix.tsv", sep="\t", index=False)

    print("\n=== Global summary ===")
    print(summary.to_string(index=False))
    print(f"\nAll outputs saved to {OUT_DIR}/")

if __name__ == "__main__":
    main()
