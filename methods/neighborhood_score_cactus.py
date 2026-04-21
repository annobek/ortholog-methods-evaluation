import pandas as pd
import pickle
import networkx as nx
from collections import defaultdict
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np



# Evaluation Cactus vs Synteny neighborhood algorithm - closer look at the results
#---------------------------------------------------------------------------------


# False negatives - the pairs that are found by Cactus, but not by Neighborhood
#-----------------

# We need to see if there is some neighborhood detected around the Cactus-only genes. 
    # -> Hypothesis: there is none, that's why it wasn't caught by neighborhood algorithm

# We need to count number of shared neighbors (Neighborhood Score - N-score) between genes for each Cactus-only pair.
# We will use count_shared_genes() function from synteny neighborhood algorithm.


# Function from the synteny algorithm (adapted for the input dicts)
def count_shared_genes(neighborhood_1, neighborhood_2, relationship_dict):
    """
    Counts the maximum number of unique gene pairs between two neighborhoods
    using the Hopcroft-Karp algorithm for bipartite matching.

    neighborhood_1: list of genes around gene in species1
    neighborhood_2: list of genes around gene_related in species2
    relationship_dict: dict gene_in_species1 -> [genes_in_species2]
    """
    neighborhood_set_2 = set(neighborhood_2)

    B = nx.Graph()
    B.add_nodes_from(neighborhood_1, bipartite=0)
    B.add_nodes_from(neighborhood_2, bipartite=1)

    for gene1 in neighborhood_1:
        related_genes = relationship_dict.get(gene1, [])
        for related_gene in related_genes:   # already a string ID
            if related_gene in neighborhood_set_2:
                B.add_edge(gene1, related_gene)

    if len(B.edges) == 0:
        return 0, 0

    try:
        # top_nodes can be either side; we choose neighborhood_1
        matching = nx.bipartite.maximum_matching(B, top_nodes=neighborhood_1)
    except Exception:
        matching = {}

    max_unique_pairs = len(matching) // 2
    return max_unique_pairs, len(B.edges)


# -------------------------------------------------------------
# 1. Load input data and define output
# -------------------------------------------------------------

# Dictionary of neighborhood for each gene
neighborhood_dict_pkl = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/neighborhoods.pkl"

# Dictionary of how close the genes are (geometric mean)
geomean_dict_pkl = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/geometric_mean_pairs_protein.pkl"

# Cactus-only pairs - false negatives (obtained from evaluation of Neighborhood vs Cactus)
false_negatives_path  = "results/evaluation/mammalia/vs_neighborhood/missed_by_nbh_gm.tsv"

# Neighborhood-only pairs 

# Map protein-gene
map_protein_genes_path = "/storage/EasyVectorOmics/synteny_algorithm/material/mammalian/gene_to_protein_map.tsv"

# Output:
# Cactus-only pairs with N-scores
out_fn_tsv = "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/missed_by_nbh_gm_with_scores.tsv"

# Plot of disctibution of N-scores
plot_distr_n_scores = "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/n_scores_distribution_cactus_vs_nbh_gm.png"

# Load all files and dicts
with open(neighborhood_dict_pkl, "rb") as f:
    nbh_dict = pickle.load(f)          # expected: {GeneID(str): [neighbor_gene1, neighbor_gene2, ...]}
with open(geomean_dict_pkl, "rb") as f:
    geomean_dict = pickle.load(f)      # expected: {(prot1, prot2): score, ...}
fn_df = pd.read_csv(false_negatives_path, sep="\t")
map_df = pd.read_csv(map_protein_genes_path, sep="\t")

# -------------------------------------------------------------
# 2. Build Protein -> Gene mapping
# -------------------------------------------------------------

# Adjust column names if needed
# expected columns: 'GeneID', 'ProteinID'
map_df["GeneID"] = map_df["GeneID"].astype(str)
map_df["ProteinID"] = map_df["ProteinID"].astype(str)

prot2gene = dict(zip(map_df["ProteinID"], map_df["GeneID"]))

# -------------------------------------------------------------
# 3. Build relationship_dict from geometric-mean protein pairs
# -------------------------------------------------------------

relationship_dict = defaultdict(list)

for (p1, p2), score in geomean_dict.items():
    p1 = str(p1)
    p2 = str(p2)
    g1 = prot2gene.get(p1)
    g2 = prot2gene.get(p2)
    if g1 is None or g2 is None:
        continue
    # make relationships symmetric
    relationship_dict[g1].append(g2)
    relationship_dict[g2].append(g1)

# quick check
print("Example entry in relationship_dict:")
some_gene = next(iter(relationship_dict))
print(some_gene, "->", relationship_dict[some_gene][:5])

# -------------------------------------------------------------
# 4. Figure out which columns hold the protein pairs (false negatives)
# -------------------------------------------------------------

if {"Query_Protein", "Target_Protein"}.issubset(fn_df.columns):
    p1_col, p2_col = "Query_Protein", "Target_Protein"
elif {"Protein1", "Protein2"}.issubset(fn_df.columns):
    p1_col, p2_col = "Protein1", "Protein2"
else:
    raise ValueError(f"Don't know which columns contain protein IDs in {false_negatives_path}; columns are {fn_df.columns.tolist()}")

fn_df[p1_col] = fn_df[p1_col].astype(str)
fn_df[p2_col] = fn_df[p2_col].astype(str)

# -------------------------------------------------------------
# 5. Compute neighborhood scores for each Cactus-only pair
# -------------------------------------------------------------

scores = []
num_edges_list = []
missing_anything = 0

for _, row in fn_df.iterrows():
    p1 = row[p1_col]
    p2 = row[p2_col]

    g1 = prot2gene.get(p1)
    g2 = prot2gene.get(p2)

    if g1 is None or g2 is None:
        # cannot map protein -> gene
        scores.append(None)
        num_edges_list.append(None)
        missing_anything += 1
        continue

    # neighborhoods: directly from nbh_dict {gene: [neighbors]}
    neigh1 = nbh_dict.get(g1, [])
    neigh2 = nbh_dict.get(g2, [])

    s, e = count_shared_genes(neigh1, neigh2, relationship_dict)
    scores.append(s)
    num_edges_list.append(e)

fn_df["NeighborhoodScore"] = scores
fn_df["NeighborhoodEdges"] = num_edges_list

print(f"Pairs where we could not map or find data: {missing_anything}")

# Save detailed FN table with scores
fn_df.to_csv(out_fn_tsv, sep="\t", index=False)
print(f"Saved scores to {out_fn_tsv}")

# -------------------------------------------------------------
# 6. Plot distribution of neighborhood scores
# -------------------------------------------------------------

valid_scores = [s for s in scores if isinstance(s, int)]

if valid_scores:
    plt.figure(figsize=(6,4))
    max_score = max(valid_scores)
    bins = range(0, max_score + 2)  # one bin per integer score

    plt.hist(valid_scores, bins=bins, align="left", edgecolor="black")
    plt.xlabel("Neighborhood score (max # of shared neighbors)")
    plt.ylabel("Number of False Negatives")
    plt.title("Neighborhood scores for False Negatives (Cactus-only pairs)")
    plt.xticks(bins)
    plt.tight_layout()
    plt.savefig(plot_distr_n_scores)

else:
    print("No valid scores to plot.")



# ----------------------------------------------------------------
# Look at contradictory results - Where Cactus and Neighborhood 
# assign different orthologs for the SAME gene
# ----------------------------------------------------------------

# We need to get N-scores for Cactus and Neighborhood orthologs

# 0) small helper: gene -> protein (for reporting)
gene2prot = dict(zip(map_df["GeneID"], map_df["ProteinID"]))

# -------------------------------------------------------------
# 1. Load Cactus orthologs
# -------------------------------------------------------------

cactus_orth_files = [
    "results/reciprocal_pairs/can_mac_one2one_fix.tsv",
    "results/reciprocal_pairs/can_rat_one2one_fix.tsv",
    "results/reciprocal_pairs/can_mus_one2one_fix.tsv",
    "results/reciprocal_pairs/mac_rat_one2one_fix.tsv",
    "results/reciprocal_pairs/mac_mus_one2one_fix.tsv",
    "results/reciprocal_pairs/rat_mus_one2one_fix.tsv",
]


# -------------------------------------------------------------
# 2. Load Neighborhood one-to-one orthologs (gene-based)
# -------------------------------------------------------------

neighborhood_orth_path = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_gm.tsv"

nbh_df = pd.read_csv(neighborhood_orth_path, sep="\t", dtype=str)

# Expect Gene1 / Gene2 columns
if not {"Gene1", "Gene2"}.issubset(nbh_df.columns):
    raise ValueError(f"Neighborhood file missing Gene1/Gene2, has {nbh_df.columns.tolist()}")

nbh_df = nbh_df.dropna(subset=["Gene1", "Gene2"])
nbh_df["Gene1"] = nbh_df["Gene1"].astype(str)
nbh_df["Gene2"] = nbh_df["Gene2"].astype(str)
nbh_df = nbh_df.drop_duplicates()

print(f"Neighborhood pairs (gene space): {len(nbh_df):,}")



# --- helper: load ONE cactus file -> Gene1/Gene2
def load_one_cactus_file_to_genes(f, map_df):
    df = pd.read_csv(f, sep="\t", dtype=str)
    required = {"Query_Protein", "Target_Protein"}
    if not required.issubset(df.columns):
        raise ValueError(f"{f} missing columns {required}, has {df.columns.tolist()}")

    df = df.merge(map_df[["ProteinID", "GeneID"]], left_on="Query_Protein",
                  right_on="ProteinID", how="left") \
           .rename(columns={"GeneID": "Gene1"}) \
           .drop(columns=["ProteinID"])

    df = df.merge(map_df[["ProteinID", "GeneID"]], left_on="Target_Protein",
                  right_on="ProteinID", how="left") \
           .rename(columns={"GeneID": "Gene2"}) \
           .drop(columns=["ProteinID"])

    df = df.dropna(subset=["Gene1", "Gene2"])
    df["Gene1"] = df["Gene1"].astype(str)
    df["Gene2"] = df["Gene2"].astype(str)
    return df[["Gene1", "Gene2"]].drop_duplicates()


# --- helper: build oriented map A->B 
def oriented_map_A_to_B(df_pairs, A_set, B_set):
    m = {}
    for g1, g2 in zip(df_pairs["Gene1"], df_pairs["Gene2"]):
        if g1 in A_set and g2 in B_set:
            m[g1] = g2
        elif g2 in A_set and g1 in B_set:
            m[g2] = g1
    return m



# -------------------------------------------------------------
# 3. Build bidirectional maps: gene -> ortholog for neighborhood
# -------------------------------------------------------------

def build_bidirectional_map(df, col1="Gene1", col2="Gene2"):
    mapping = {}
    for g1, g2 in zip(df[col1], df[col2]):
        mapping[g1] = g2
        mapping[g2] = g1
    return mapping

nbh_map = build_bidirectional_map(nbh_df, "Gene1", "Gene2")

print(f"Genes with Neighborhood assignments: {len(nbh_map):,}")

# -------------------------------------------------------------
# 4. Identify contradictory assignments
# -------------------------------------------------------------


OUT_PLOTS_PREFIX = "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/contradictory_cactus_vs_nbh_gm"
out_contr_tsv = "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/contradictory_results_gm.tsv"

results = []

for cactus_file in cactus_orth_files:
    pair_name = Path(cactus_file).stem.replace("_one2one_fix", "")

    cactus_pair = load_one_cactus_file_to_genes(cactus_file, map_df)

    # Sets of genes of both species (how they are in cactus_file)
    A = set(cactus_pair["Gene1"])
    B = set(cactus_pair["Gene2"])

    # Neighborhood pairs only between the same two sets (so between same species)
    nbh_pair = nbh_df[
        (nbh_df["Gene1"].isin(A) & nbh_df["Gene2"].isin(B)) |
        (nbh_df["Gene1"].isin(B) & nbh_df["Gene2"].isin(A))
    ].drop_duplicates().copy()

    # check: inside species-pair there should be close to 0
    print(f"[{pair_name}] cactus Gene1>1 partner:",
          (cactus_pair.groupby("Gene1")["Gene2"].nunique() > 1).sum())

    # Build mappings
    cactus_map_A = cactus_pair.drop_duplicates("Gene1").set_index("Gene1")["Gene2"].to_dict()
    nbh_map_A    = oriented_map_A_to_B(nbh_pair, A, B)

    common = set(cactus_map_A.keys()) & set(nbh_map_A.keys())
    contradictory = [(g, cactus_map_A[g], nbh_map_A[g]) for g in common if cactus_map_A[g] != nbh_map_A[g]]

    print(f"[{pair_name}] contradictory genes: {len(contradictory)}")

    # Count N-score for contradictory cases
    for gene, cactus_t, neigh_t in contradictory:
        if gene not in nbh_dict or cactus_t not in nbh_dict or neigh_t not in nbh_dict:
            continue

        nbh_gene   = nbh_dict[gene]
        nbh_cactus = nbh_dict[cactus_t]
        nbh_neigh  = nbh_dict[neigh_t]

        score_C, edges_C = count_shared_genes(nbh_gene, nbh_cactus, relationship_dict)
        score_N, edges_N = count_shared_genes(nbh_gene, nbh_neigh, relationship_dict)

        max_C = min(len(nbh_gene), len(nbh_cactus))
        max_N = min(len(nbh_gene), len(nbh_neigh))
        rel_C = score_C / max_C if max_C > 0 else None
        rel_N = score_N / max_N if max_N > 0 else None

        # Important:
        
        # Caclulate ratio (can be >1)
        ratio = score_C / score_N if score_N > 0 else None

        # Calculate fraction (is in interval [0..1])
        frac01 = score_C / (score_C + score_N) if (score_C + score_N) > 0 else None

        results.append({
            "Pair": pair_name,
            "Gene": gene,
            "Cactus_Gene": cactus_t,
            "Neighborhood_Gene": neigh_t,
            "Cactus_Protein": gene2prot.get(cactus_t, ""),
            "Neighborhood_Protein": gene2prot.get(neigh_t, ""),
            "Nscore_C": score_C,
            "Edges_C": edges_C,
            "Nscore_N": score_N,
            "Edges_N": edges_N,
            "RelScore_C": rel_C,
            "RelScore_N": rel_N,
            "Ratio_C_over_N": ratio,
            "Fraction_C": frac01,
            "MaxPossible_C": max_C,
            "MaxPossible_N": max_N,
        })

print(f"Contradictory cases with full neighborhood info: {len(results):,}")

res_df = pd.DataFrame(results)
res_df.to_csv(out_contr_tsv, sep="\t", index=False)
print(f"Saved detailed contradictory table to: {out_contr_tsv}")


# -----------------------------------------------------------
# 6. Plots: distributions of Nscore_C, Nscore_N, and Fraction
#------------------------------------------------------------

valid_C = res_df["Nscore_C"].dropna()
valid_N = res_df["Nscore_N"].dropna()
valid_frac = res_df["Fraction_C"].dropna()

plt.figure(figsize=(10,4))
plt.hist(valid_C, bins=20)
plt.xlabel("Neighborhood score for Cactus ortholog (Nscore_C)")
plt.ylabel("Number of genes")
plt.title("N-score for Cactus ortholog in contradictory cases")
plt.tight_layout()
plt.savefig(OUT_PLOTS_PREFIX + "_Nscore_C.png")
plt.close()

plt.figure(figsize=(10,4))
plt.hist(valid_N, bins=20)
plt.xlabel("Neighborhood score for Neighborhood ortholog (Nscore_N)")
plt.ylabel("Number of genes")
plt.title("N-score for Neighborhood ortholog in contradictory cases")
plt.tight_layout()
plt.savefig(OUT_PLOTS_PREFIX + "_Nscore_N.png")
plt.close()

plt.figure(figsize=(10,4))

# bins exactly on [0,1]
bin_width = 0.02
bins = np.arange(-bin_width/2, 1 + bin_width/2 + 1e-9, bin_width)

plt.hist(valid_frac, bins=bins)
plt.xlim(0, 1)
plt.xticks(np.arange(0, 1.01, 0.1))

plt.xlabel("Fraction = score_C / (score_C + score_N)")
plt.title("Fraction in [0..1] in contradictory cases")
plt.ylabel("Number of cases")
plt.tight_layout()
plt.savefig(OUT_PLOTS_PREFIX + "_Fraction.png")
plt.close()


print("Plots saved with prefix:", OUT_PLOTS_PREFIX)

print("Some final checks:")
eq = (res_df["Nscore_C"] == res_df["Nscore_N"]).mean()
c_better = (res_df["Nscore_C"] > res_df["Nscore_N"]).mean()
n_better = (res_df["Nscore_N"] > res_df["Nscore_C"]).mean()
print("ties:", eq, "C better:", c_better, "N better:", n_better)


