import pandas as pd
import pickle
import networkx as nx
from collections import defaultdict
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns


# =============================================================================
# Evaluation Cactus vs prot-syn - False Negatives Analysis
# =============================================================================

# False negatives are pairs found by Cactus but not by prot-syn
# We classify them as:
#   - PS-no-results: prot-syn found no ortholog for this query gene
#   - PS-contradictory: prot-syn assigned a different target gene

# =============================================================================
# File Paths
# =============================================================================

# Input files
neighborhood_dict_pkl = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/neighborhoods.pkl"
geomean_dict_pkl = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/geometric_mean.pkl"
map_protein_genes_path = "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/gene_to_protein_map_FIXED.tsv"

# Prot-syn ortholog file
protsyn_orth_path = "material/sex_experiment/one_to_one_sum.tsv"

# Cactus ortholog files (one per species pair)
cactus_orth_files = [
    "results/reciprocal_pairs/can_mac_one2one_fix.tsv",
    "results/reciprocal_pairs/can_rat_one2one_fix.tsv",
    "results/reciprocal_pairs/can_mus_one2one_fix.tsv",
    "results/reciprocal_pairs/mac_rat_one2one_fix.tsv",
    "results/reciprocal_pairs/mac_mus_one2one_fix.tsv",
    "results/reciprocal_pairs/rat_mus_one2one_fix.tsv",
]

# False negatives files (one per species pair)
false_negatives_files = [
    "results/evaluation/mammalia/vs_prot_syn/sum_approach/false_negatives_can_mac.tsv",
    "results/evaluation/mammalia/vs_prot_syn/sum_approach/false_negatives_can_rat.tsv",
    "results/evaluation/mammalia/vs_prot_syn/sum_approach/false_negatives_can_mus.tsv",
    "results/evaluation/mammalia/vs_prot_syn/sum_approach/false_negatives_mac_rat.tsv",
    "results/evaluation/mammalia/vs_prot_syn/sum_approach/false_negatives_mac_mus.tsv",
    "results/evaluation/mammalia/vs_prot_syn/sum_approach/false_negatives_rat_mus.tsv",
]

# Output directory
out_dir = Path("results/evaluation/mammalia/vs_prot_syn/sum_approach/disagreement_investigation/false_negatives")
out_dir.mkdir(parents=True, exist_ok=True)

# Output files
out_fn_with_scores = out_dir / "false_negatives_with_scores_all.tsv"
fn_classified_out = out_dir / "false_negatives_classified_with_nscores_all.tsv"
out_path_better = out_dir / "cases_cactus_picked_better_n_score.tsv"

# Plot files
plot_distr_n_scores = out_dir / "n_scores_distribution_cactus_vs_prot_syn.png"
OUT_PLOTS_PREFIX = str(out_dir / "contradictory_cactus_vs_prot_syn")


# =============================================================================
# Helper Functions
# =============================================================================

def count_shared_genes(gene, gene_related, neighborhood_1, neighborhood_2, relationship_dict):
    """
    Count maximum number of unique gene pairs between two neighborhoods
    using bipartite matching.
    
    Args:
        gene: Query gene (for compatibility, not used)
        gene_related: Related gene (for compatibility, not used)
        neighborhood_1: List of genes around gene in species1
        neighborhood_2: List of genes around gene_related in species2
        relationship_dict: Dict of gene relationships {gene: [(related_gene, score), ...]}
    
    Returns:
        (max_unique_pairs, num_edges): Tuple of matching size and edge count
    """
    neighborhood_set_2 = set(neighborhood_2)
    
    # Create bipartite graph
    B = nx.Graph()
    B.add_nodes_from(neighborhood_1, bipartite=0)
    B.add_nodes_from(neighborhood_2, bipartite=1)
    
    # Add edges based on relationship dictionary
    for gene1 in neighborhood_1:
        if gene1 in neighborhood_set_2:
            continue
            
        related_genes = relationship_dict.get(gene1)
        if not related_genes:
            continue
            
        for gene2, _ in related_genes:  # Unpack tuple (gene, score)
            if gene2 in neighborhood_set_2:
                B.add_edge(gene1, gene2)
    
    if B.number_of_edges() == 0:
        return 0, 0
    
    try:
        matching = nx.bipartite.maximum_matching(B, top_nodes=neighborhood_2)
    except Exception:
        matching = {}
    
    max_unique_pairs = len(matching) // 2
    return max_unique_pairs, B.number_of_edges()


def load_cactus_file_to_genes(filepath, map_df):
    """
    Load Cactus ortholog file and convert proteins to genes.
    
    Returns:
        DataFrame with Gene1, Gene2 columns
    """
    df = pd.read_csv(filepath, sep="\t", dtype=str)
    
    if not {"Query_Protein", "Target_Protein"}.issubset(df.columns):
        raise ValueError(f"File {filepath} missing Query_Protein/Target_Protein columns")
    
    # Map proteins to genes
    df = df.merge(map_df[["ProteinID", "GeneID"]], 
                  left_on="Query_Protein", right_on="ProteinID", how="left") \
           .rename(columns={"GeneID": "Gene1"}) \
           .drop(columns=["ProteinID"])
    
    df = df.merge(map_df[["ProteinID", "GeneID"]], 
                  left_on="Target_Protein", right_on="ProteinID", how="left") \
           .rename(columns={"GeneID": "Gene2"}) \
           .drop(columns=["ProteinID"])
    
    df = df.dropna(subset=["Gene1", "Gene2"])
    df["Gene1"] = df["Gene1"].astype(str)
    df["Gene2"] = df["Gene2"].astype(str)
    
    return df[["Gene1", "Gene2"]].drop_duplicates()


def get_species_sets(map_df, species_pair_name):
    """
    Extract species sets A and B from the map based on species pair name.
    
    Args:
        map_df: Mapping dataframe with GeneID and Species columns
        species_pair_name: e.g., 'can_mac'
    
    Returns:
        (set_A, set_B, code_A, code_B): Two sets of gene IDs and their species codes
    """
    species_codes = species_pair_name.split("_")
    if len(species_codes) != 2:
        raise ValueError(f"Cannot parse species from {species_pair_name}")
    
    code_to_prefix = {
        "can": "canis lupus",
        "mac": "macaca fascicularis",
        "mus": "mus musculus",
        "rat": "rattus norvegicus",
    }
    
    sp1_key = code_to_prefix.get(species_codes[0])
    sp2_key = code_to_prefix.get(species_codes[1])
    
    if sp1_key is None or sp2_key is None:
        raise ValueError(f"Unknown species code in {species_pair_name}")
    
    def matches_species(species_str, target):
        if pd.isna(species_str):
            return False
        return species_str.lower().startswith(target.lower())
    
    A = set(map_df[map_df["Species"].apply(
        lambda x: matches_species(x, sp1_key))]["GeneID"].astype(str))
    B = set(map_df[map_df["Species"].apply(
        lambda x: matches_species(x, sp2_key))]["GeneID"].astype(str))
    
    # Verify no overlap
    overlap = A & B
    if overlap:
        raise ValueError(f"Species sets {species_pair_name} overlap! {len(overlap)} genes in both.")
    
    return A, B, species_codes[0], species_codes[1]


def build_oriented_map(pairs_df, set_A, set_B):
    """
    Build a mapping from genes in set_A to genes in set_B.
    
    Args:
        pairs_df: DataFrame with Gene1, Gene2 columns
        set_A: Set of genes in species A
        set_B: Set of genes in species B
    
    Returns:
        Dict mapping gene_A -> gene_B
    """
    mapping = {}
    for g1, g2 in zip(pairs_df["Gene1"], pairs_df["Gene2"]):
        if g1 in set_A and g2 in set_B:
            mapping[g1] = g2
        elif g2 in set_A and g1 in set_B:
            mapping[g2] = g1
    return mapping


# =============================================================================
# Load Data
# =============================================================================

print("Loading data...")

# Load dictionaries
with open(neighborhood_dict_pkl, "rb") as f:
    nbh_dict = pickle.load(f)  # {GeneID: [neighbor_gene1, neighbor_gene2, ...]}

with open(geomean_dict_pkl, "rb") as f:
    geomean_dict = pickle.load(f)  # {(prot1, prot2): score, ...}

# Load mapping table
map_df = pd.read_csv(map_protein_genes_path, sep="\t")
map_df["GeneID"] = map_df["GeneID"].astype(str)
map_df["ProteinID"] = map_df["ProteinID"].astype(str)

# Load prot-syn orthologs
protsyn_df = pd.read_csv(protsyn_orth_path, sep="\t", dtype=str)

print(f"Loaded {len(map_df)} protein-gene mappings")
print(f"Loaded {len(protsyn_df)} prot-syn ortholog pairs")


# =============================================================================
# Build Mappings
# =============================================================================

# Protein -> Gene mapping
prot2gene = dict(zip(map_df["ProteinID"], map_df["GeneID"]))
gene2prot = dict(zip(map_df["GeneID"], map_df["ProteinID"]))

# The geomean_dict is ALREADY in the format we need!
# Structure: {gene: [(related_gene, score), ...]}
# Just ensure keys are strings
relationship_dict = {}
for gene, related_list in geomean_dict.items():
    gene_str = str(gene)
    # Ensure all related genes are also strings
    relationship_dict[gene_str] = [(str(rel_gene), score) for rel_gene, score in related_list]

print(f"Loaded relationship dict with {len(relationship_dict)} genes")

# Quick verification
sample_gene = next(iter(relationship_dict))
print(f"Sample entry: {sample_gene} -> {relationship_dict[sample_gene][:3]}")


# =============================================================================
# Process Prot-syn Orthologs - Build Maps per Species Pair
# =============================================================================

if not {"Gene1", "Gene2"}.issubset(protsyn_df.columns):
    raise ValueError(f"Prot-syn file missing Gene1/Gene2 columns")

protsyn_df = protsyn_df.dropna(subset=["Gene1", "Gene2"])
protsyn_df["Gene1"] = protsyn_df["Gene1"].astype(str)
protsyn_df["Gene2"] = protsyn_df["Gene2"].astype(str)
protsyn_df = protsyn_df.drop_duplicates()

print(f"\nBuilding prot-syn maps per species pair...")

protsyn_maps_by_pair = {}

for cactus_file in cactus_orth_files:
    # Extract species pair name from filename
    pair_name_raw = Path(cactus_file).stem.replace("_one2one_fix", "")
    
    # Normalize pair name (sort species codes)
    species_codes = pair_name_raw.split("_")
    pair_name = "_".join(sorted(species_codes))
    
    print(f"  Processing {pair_name_raw} → {pair_name}")
    
    # Get species sets (use raw name for lookup, it has original order)
    A, B, sp_code_A, sp_code_B = get_species_sets(map_df, pair_name_raw)
    print(f"    Species A ({sp_code_A}): {len(A)} genes")
    print(f"    Species B ({sp_code_B}): {len(B)} genes")
    
    # Filter prot-syn for this species pair
    protsyn_pair = protsyn_df[
        (protsyn_df["Gene1"].isin(A) & protsyn_df["Gene2"].isin(B)) |
        (protsyn_df["Gene1"].isin(B) & protsyn_df["Gene2"].isin(A))
    ].drop_duplicates().copy()
    
    print(f"    Prot-syn orthologs: {len(protsyn_pair)}")
    
    # Build oriented map A -> B
    protsyn_maps_by_pair[pair_name] = build_oriented_map(protsyn_pair, A, B)
    print(f"    Map size: {len(protsyn_maps_by_pair[pair_name])}")


# =============================================================================
# Process All FN Files
# =============================================================================

print(f"\n{'='*60}")
print("Processing False Negatives Files")
print(f"{'='*60}")

all_fn_with_scores = []
all_classified = []

for fn_file in false_negatives_files:
    pair_name_raw = Path(fn_file).stem.replace("false_negatives_", "")
    species_codes = pair_name_raw.split("_")
    pair_name = "_".join(sorted(species_codes))
    
    print(f"\n--- Processing {pair_name} ---")
    
    # Load FN file
    fn_df = pd.read_csv(fn_file, sep="\t")
    print(f"Loaded {len(fn_df)} false negatives")
    
    # Check required columns
    required_cols = {"Query_Protein", "Target_Protein", "Query_Gene", "Target_Gene", 
                     "Query_Species", "Target_Species"}
    if not required_cols.issubset(fn_df.columns):
        print(f"WARNING: Missing columns in {fn_file}, skipping")
        continue
    
    fn_df["Query_Protein"] = fn_df["Query_Protein"].astype(str)
    fn_df["Target_Protein"] = fn_df["Target_Protein"].astype(str)
    fn_df["Query_Gene"] = fn_df["Query_Gene"].astype(str)
    fn_df["Target_Gene"] = fn_df["Target_Gene"].astype(str)
    
    # -------------------------------------------------------------------------
    # Compute N-scores for all FN pairs
    # -------------------------------------------------------------------------
    
    scores = []
    num_edges_list = []
    missing_count = 0
    
    for _, row in fn_df.iterrows():
        g1 = row["Query_Gene"]
        g2 = row["Target_Gene"]
        
        neigh1 = nbh_dict.get(g1, [])
        neigh2 = nbh_dict.get(g2, [])
        
        if not neigh1 or not neigh2:
            scores.append(None)
            num_edges_list.append(None)
            missing_count += 1
        else:
            s, e = count_shared_genes(g1, g2, neigh1, neigh2, relationship_dict)
            scores.append(s)
            num_edges_list.append(e)
    
    fn_df["NeighborhoodScore"] = scores
    fn_df["NeighborhoodEdges"] = num_edges_list
    fn_df["Pair"] = pair_name
    
    print(f"Could not compute scores for {missing_count} pairs")
    
    all_fn_with_scores.append(fn_df)
    
    # -------------------------------------------------------------------------
    # Classify FN: No-Results vs Contradictory
    # -------------------------------------------------------------------------
    
    if pair_name not in protsyn_maps_by_pair:
        print(f"WARNING: No prot-syn map for {pair_name}, skipping classification")
        continue
    
    protsyn_map = protsyn_maps_by_pair[pair_name]
    
    for _, row in fn_df.iterrows():
        gA = row["Query_Gene"]
        gB = row["Target_Gene"]
        pA = row["Query_Protein"]
        pB = row["Target_Protein"]
        
        # Check if prot-syn has an assignment for gA
        gX = protsyn_map.get(gA)
        
        # Classify
        if gX is None:
            fn_type = "PS-no-results"
        else:
            fn_type = "PS-contradictory"
        
        # Compute N-scores
        nscore_cactus = None
        nscore_protsyn = None
        
        # N-score for Cactus pair (gA -> gB)
        if gA in nbh_dict and gB in nbh_dict:
            nscore_cactus, _ = count_shared_genes(
                gA, gB, nbh_dict[gA], nbh_dict[gB], relationship_dict
            )
        
        # N-score for prot-syn pair (gA -> gX), if exists
        if gX is not None and gA in nbh_dict and gX in nbh_dict:
            nscore_protsyn, _ = count_shared_genes(
                gA, gX, nbh_dict[gA], nbh_dict[gX], relationship_dict
            )
        
        # Compute fraction
        if nscore_cactus is not None and nscore_protsyn is not None and nscore_protsyn > 0:
            fraction = nscore_cactus / nscore_protsyn
        else:
            fraction = None
        
        pX = gene2prot.get(gX) if gX is not None else None
        
        all_classified.append({
            "Pair": pair_name,
            "Query_Protein": pA,
            "Target_Protein": pB,
            "Query_Gene": gA,
            "Cactus_Target_Gene": gB,
            "ProtSyn_Target_Gene": gX if gX else "",
            "ProtSyn_Target_Protein": pX if pX else "",
            "Query_Species": row["Query_Species"],
            "Target_Species": row["Target_Species"],
            "FN_Type": fn_type,
            "Nscore_Cactus": nscore_cactus,
            "Nscore_ProtSyn": nscore_protsyn,
            "Fraction_Cactus_ProtSyn": fraction,
        })


# =============================================================================
# Combine All Results
# =============================================================================

print(f"\n{'='*60}")
print("Combining Results")
print(f"{'='*60}")

# Combine all FN with scores
all_fn_df = pd.concat(all_fn_with_scores, ignore_index=True)
all_fn_df.to_csv(out_fn_with_scores, sep="\t", index=False)
print(f"✓ Saved all FN with scores: {out_fn_with_scores}")
print(f"  Total FN pairs: {len(all_fn_df)}")

# Combine all classified FN
classified_df = pd.DataFrame(all_classified)
classified_df.to_csv(fn_classified_out, sep="\t", index=False)
print(f"✓ Saved all classified FN: {fn_classified_out}")
print(f"  Total classified: {len(classified_df)}")

print("\nClassification summary:")
print(classified_df["FN_Type"].value_counts())

print("\nClassification by species pair:")
print(classified_df.groupby("Pair")["FN_Type"].value_counts())


# =============================================================================
# Plot N-score Distribution (All FN)
# =============================================================================

valid_scores = all_fn_df["NeighborhoodScore"].dropna()

if not valid_scores.empty:
    plt.figure(figsize=(8, 5))
    max_score = int(valid_scores.max())
    bins = range(0, max_score + 2)
    
    plt.hist(valid_scores, bins=bins, align="left", edgecolor="black", color="steelblue")
    plt.xlabel("Neighborhood score (max # of shared neighbors)", fontsize=12)
    plt.ylabel("Number of False Negatives", fontsize=12)
    plt.title("Neighborhood scores for False Negatives (Cactus-only pairs)", fontsize=13)
    plt.xticks(bins)
    plt.grid(axis='y', alpha=0.3, linestyle='--')
    plt.tight_layout()
    plt.savefig(plot_distr_n_scores, dpi=300)
    plt.close()
    print(f"\n✓ Saved N-score distribution plot: {plot_distr_n_scores}")
else:
    print("\nNo valid scores to plot.")


# =============================================================================
# Find Cases Where Cactus Picked Better N-score
# =============================================================================

print(f"\n{'='*60}")
print("Finding Cases Where Cactus Picked Better N-score")
print(f"{'='*60}")

cases_better = classified_df[classified_df["FN_Type"] == "PS-contradictory"].copy()

# Ensure numeric
cases_better["Nscore_Cactus"] = pd.to_numeric(cases_better["Nscore_Cactus"], errors="coerce")
cases_better["Nscore_ProtSyn"] = pd.to_numeric(cases_better["Nscore_ProtSyn"], errors="coerce")

# Keep only valid cases
cases_better = cases_better.dropna(subset=["Nscore_Cactus", "Nscore_ProtSyn"])
cases_better = cases_better[(cases_better["Nscore_Cactus"] > 0) & (cases_better["Nscore_ProtSyn"] > 0)]

# Recalculate fraction
cases_better["Fraction_Cactus_ProtSyn"] = cases_better["Nscore_Cactus"] / cases_better["Nscore_ProtSyn"]

# Filter where Cactus is better
cases_better = cases_better[cases_better["Fraction_Cactus_ProtSyn"] > 1].copy()
cases_better = cases_better.sort_values("Fraction_Cactus_ProtSyn", ascending=False)

# Save
cases_better.to_csv(out_path_better, sep="\t", index=False)

contradictory_total = len(classified_df[classified_df["FN_Type"] == "PS-contradictory"])
pct = (len(cases_better) / contradictory_total * 100) if contradictory_total > 0 else 0

print(f"✓ Saved cases where Cactus picked better: {out_path_better}")
print(f"  Total cases: {len(cases_better)}")
print(f"  Percentage of contradictory: {pct:.1f}%")

if not cases_better.empty:
    print(f"\nTop 10 cases where Cactus picked much better ortholog:")
    print(cases_better[[
        "Pair", "Query_Gene", "Cactus_Target_Gene", "ProtSyn_Target_Gene",
        "Nscore_Cactus", "Nscore_ProtSyn", "Fraction_Cactus_ProtSyn"
    ]].head(10).to_string(index=False))

print(f"\nCases where Cactus picked better, by species pair:")
print(cases_better["Pair"].value_counts())


# =============================================================================
# Density Plots for Contradictory Cases
# =============================================================================

print(f"\n{'='*60}")
print("Generating Density Plots")
print(f"{'='*60}")

contradictory = classified_df[classified_df["FN_Type"] == "PS-contradictory"]

valid_C = contradictory["Nscore_Cactus"].dropna()
valid_N = contradictory["Nscore_ProtSyn"].dropna()
valid_frac = contradictory["Fraction_Cactus_ProtSyn"].dropna()

# Plot 1: Cactus N-scores
if not valid_C.empty:
    plt.figure(figsize=(12, 5))
    sns.kdeplot(data=valid_C, fill=True, color='steelblue', alpha=0.6, linewidth=2)
    plt.xlabel("Neighborhood score for Cactus ortholog (Nscore_Cactus)", fontsize=12)
    plt.ylabel("Density", fontsize=12)
    plt.title("N-score distribution for Cactus orthologs in contradictory cases", fontsize=13)
    plt.xlim(left=0)
    plt.grid(axis='both', alpha=0.3, linestyle='--')
    plt.tight_layout()
    plt.savefig(OUT_PLOTS_PREFIX + "_Nscore_Cactus_density.png", dpi=300)
    plt.close()

# Plot 2: Prot-syn N-scores
if not valid_N.empty:
    plt.figure(figsize=(12, 5))
    sns.kdeplot(data=valid_N, fill=True, color='coral', alpha=0.6, linewidth=2)
    plt.xlabel("Neighborhood score for prot-syn ortholog (Nscore_ProtSyn)", fontsize=12)
    plt.ylabel("Density", fontsize=12)
    plt.title("N-score distribution for prot-syn orthologs in contradictory cases", fontsize=13)
    plt.xlim(left=0)
    plt.grid(axis='both', alpha=0.3, linestyle='--')
    plt.tight_layout()
    plt.savefig(OUT_PLOTS_PREFIX + "_Nscore_ProtSyn_density.png", dpi=300)
    plt.close()

# Plot 3: Fraction
if not valid_frac.empty:
    plt.figure(figsize=(12, 5))
    sns.kdeplot(data=valid_frac, fill=True, color='mediumseagreen', alpha=0.6, linewidth=2)
    plt.xlabel("Fraction = Nscore_Cactus / Nscore_ProtSyn", fontsize=12)
    plt.ylabel("Density", fontsize=12)
    plt.title("Fraction distribution in contradictory cases", fontsize=13)
    plt.axvline(x=1.0, color='red', linestyle='--', linewidth=2.5,
                label='Equal scores (Fraction=1)', zorder=10)
    plt.xlim(left=0)
    plt.legend(fontsize=11, loc='upper right')
    plt.grid(axis='both', alpha=0.3, linestyle='--')
    plt.tight_layout()
    plt.savefig(OUT_PLOTS_PREFIX + "_Fraction_density.png", dpi=300)
    plt.close()

print(f"✓ Saved 3 density plots to {OUT_PLOTS_PREFIX}_*_density.png")


# =============================================================================
# Summary Statistics
# =============================================================================

print(f"\n{'='*60}")
print("SUMMARY STATISTICS")
print(f"{'='*60}")

if not valid_C.empty:
    print(f"\nCactus N-scores:")
    print(f"  Mean:      {valid_C.mean():.2f}")
    print(f"  Median:    {valid_C.median():.2f}")
    print(f"  N-score=0: {(valid_C == 0).sum()} ({(valid_C == 0).sum()/len(valid_C)*100:.1f}%)")

if not valid_N.empty:
    print(f"\nProt-syn N-scores:")
    print(f"  Mean:      {valid_N.mean():.2f}")
    print(f"  Median:    {valid_N.median():.2f}")
    print(f"  N-score=0: {(valid_N == 0).sum()} ({(valid_N == 0).sum()/len(valid_N)*100:.1f}%)")

if not valid_frac.empty:
    print(f"\nFraction (Cactus/ProtSyn):")
    print(f"  Mean:              {valid_frac.mean():.2f}")
    print(f"  Median:            {valid_frac.median():.2f}")
    print(f"  Cactus better (>1):  {(valid_frac > 1).sum()} ({(valid_frac > 1).sum()/len(valid_frac)*100:.1f}%)")
    print(f"  ProtSyn better (<1): {(valid_frac < 1).sum()} ({(valid_frac < 1).sum()/len(valid_frac)*100:.1f}%)")
    print(f"  Equal (=1):          {(valid_frac == 1).sum()} ({(valid_frac == 1).sum()/len(valid_frac)*100:.1f}%)")

print(f"\n{'='*60}")
print("✓ Analysis Complete!")
print(f"{'='*60}")