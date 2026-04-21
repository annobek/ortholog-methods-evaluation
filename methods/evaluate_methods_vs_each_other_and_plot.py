import pandas as pd
import itertools
import numpy as np
import matplotlib.pyplot as plt
import sys

# Evaluate methods vs each other and visualize (jaccard similarity only)

# for mammals
'''
# Input
methods = {
    "Neighborhood":"material/sex_experiment/global_pairs_gm.tsv",
    "Tree-majority":"/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/majority/complete_tree_conserved_orthologs_2.tsv",
    "Tree-standard":"/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/standard/complete_tree_conserved_orthologs_2.tsv",
    "Tree-whitelist":"/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/whitelist/complete_tree_conserved_orthologs_2.tsv",
    "PTP":"/storage/EasyVectorOmics/phylotreepruner/results/ptp_one_to_one_pairs_new.tsv",
    "OrthoFinder":"results/evaluation/vs_orthofinder/orthofinder_1to1.tsv"
}
# Output:
OUT_METRICS = "results/evaluation/all_vs_all/all_vs_all_jaccard.txt"
OUT_PLOT = "results/evaluation/all_vs_all_jaccard_heatmap.png"
'''

# for plants
'''
# Input
methods = {
    "Neighborhood-global-pairs":"/storage/EasyVectorOmics/synteny_algorithm/results/cardamine/global_pairs_gm.tsv",
    "Neighborhood-one-to-one":"/storage/EasyVectorOmics/synteny_algorithm/results/cardamine/one_to_one_gm.tsv",
    "Original-orthologs":"material/cardamine/filtered_orthologs.tsv",
    "PTP":"/storage/EasyVectorOmics/phylotreepruner/results/plants/ptp_one_to_one_pairs.tsv",
    "OrthoFinder":"results/evaluation/plants/vs_orthofinder/orthofinder_1to1.tsv"
}
# Output:
OUT_METRICS = "results/evaluation/plants/all_vs_all/all_vs_all_jaccard.txt"
OUT_PLOT = "results/evaluation/plants/all_vs_all_jaccard_heatmap.png"
'''
# for drosophila
# Input
methods = {
    "Neighborhood-global-pairs":"/storage/EasyVectorOmics/synteny_algorithm/results/drosophila/global_pairs_gm.tsv",
    "Neighborhood-one-to-one":"/storage/EasyVectorOmics/synteny_algorithm/results/drosophila/one_to_one_gm.tsv",
    "OrthoFinder":"results/evaluation/drosophila/vs_orthofinder/orthofinder_1to1.tsv"
}
# Output:
OUT_METRICS = "results/evaluation/drosophila/all_vs_all/all_vs_all_jaccard.txt"
OUT_PLOT = "results/evaluation/drosophila/all_vs_all_jaccard_heatmap.png"

def normalize_id(x: str):
    """Make IDs comparable across files."""
    if pd.isna(x):
        return None
    x = str(x).strip()
    x = x.replace("protein", "")        # remove literal 'protein' if present
    x = x.replace(" ", "_")             # unify spaces/underscores
    return x

def load_ortholog_pairs(path, sep="\t"):
    """
    Load ortholog pairs from any supported file format.
    Returns standardized columns: ['Species1', 'Protein1', 'Species2', 'Protein2'].
    Prefers Protein columns over Gene columns.
    """
    df = pd.read_csv(path, sep=sep, dtype=str)
    df.columns = [c.strip() for c in df.columns]  # remove stray spaces
    cols = {c.lower(): c for c in df.columns}

    # --- Detect columns ---
    species_cols = [cols[c] for c in cols if "species" in c]
    prot_cols = [cols[c] for c in cols if "protein" in c]
    gene_cols = [cols[c] for c in cols if "gene" in c]

    # --- Prefer Protein columns if available ---
    if len(prot_cols) >= 2:
        use_prot = prot_cols
    elif len(gene_cols) >= 2:
        print(f"No Protein columns found in {path}, using Gene columns instead.")
        use_prot = gene_cols
    else:
        raise ValueError(f"Could not detect suitable ID columns in {path}. Found: {df.columns.tolist()}")

    # --- Create output ---
    df_out = df[[species_cols[0], use_prot[0], species_cols[1], use_prot[1]]].copy()
    df_out.columns = ["Species1", "Protein1", "Species2", "Protein2"]

    # --- Normalize ---
    df_out["Protein1"] = df_out["Protein1"].astype(str).str.strip()
    df_out["Protein2"] = df_out["Protein2"].astype(str).str.strip()

    # --- Drop invalid / self pairs ---
    df_out = df_out[df_out["Protein1"] != df_out["Protein2"]].drop_duplicates()

    return df_out



def to_pair_set(path, a="Protein1", b="Protein2", sep="\t"):
    df = pd.read_csv(path, sep=sep)
    return {tuple(sorted((x, y))) for x, y in zip(df[a], df[b]) if pd.notna(x) and pd.notna(y)}

# --- Load all sets ---
pair_sets = {}

for name, path in methods.items():
    df = load_ortholog_pairs(path)
    pair_sets[name] = {tuple(sorted((p1, p2))) for p1, p2 in zip(df["Protein1"], df["Protein2"])}
    print(f"{name}: loaded {len(pair_sets[name]):,} pairs")

for name, s in pair_sets.items():
    ids = sorted({x for pair in s for x in pair})
    print(f"\n{name} sample IDs:")
    print(ids[:10])

# --- Compute pairwise Jaccard similarities ---
names = list(pair_sets.keys())
n = len(names)
matrix = np.zeros((n, n))

for i, j in itertools.combinations(range(n), 2):
    a, b = pair_sets[names[i]], pair_sets[names[j]]
    inter = len(a & b)
    union = len(a | b)
    matrix[i, j] = inter / union if union else 0
    matrix[j, i] = matrix[i, j]

np.fill_diagonal(matrix, 1.0)
df_jaccard = pd.DataFrame(matrix, index=names, columns=names)

sys.stdout = open(OUT_METRICS, "w")

print("\n=== Jaccard Similarity Matrix ===")
print(df_jaccard.round(3))

sys.stdout.close()
sys.stdout = sys.__stdout__

# --- Visualize as heatmap (upper triangle) ---
fig, ax = plt.subplots(figsize=(8,6))
im = ax.imshow(matrix, cmap="coolwarm", vmin=0, vmax=1)

ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels(names, rotation=45, ha="right")
ax.set_yticklabels(names)

# Hide lower triangle
for i in range(n):
    for j in range(n):
        if i > j:
            ax.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, color='white', ec='none'))

# Annotate with numbers
for i in range(n):
    for j in range(n):
        if i <= j:
            ax.text(j, i, f"{matrix[i,j]:.2f}", ha="center", va="center", color="black")

plt.colorbar(im, ax=ax, label="Jaccard Similarity")
ax.set_title("Pairwise Jaccard Similarity Between Orthology Methods")
plt.tight_layout()
plt.savefig(OUT_PLOT)
plt.show()
