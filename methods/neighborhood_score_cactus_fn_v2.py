import pandas as pd
import pickle
import networkx as nx
from collections import defaultdict
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np


# =============================================================================
# File Paths Configuration
# =============================================================================

# Input files
NEIGHBORHOOD_DICT_PKL = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/neighborhoods.pkl"
GEOMEAN_DICT_PKL = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/geometric_mean.pkl"
MAP_PROTEIN_GENES_PATH = "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/gene_to_protein_map_FIXED.tsv"
PROTSYN_ORTH_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_nscore.tsv"
FALSE_NEGATIVES_PATH = "results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/all_false_negatives_nscore.tsv"

# Output directory
OUT_DIR = Path("results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/false_negatives")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Output files
FN_CLASSIFIED_OUT = OUT_DIR / "false_negatives_classified.tsv"
FN_WITH_SCORES_OUT = OUT_DIR / "false_negatives_with_all_scores.tsv"
OUT_PATH_CACTUS_BETTER = OUT_DIR / "contradictory_cactus_higher_nscore_than_protsyn.tsv"
PLOT_DISTRIBUTIONS = OUT_DIR / "fn_distributions_2x2.png"
OUT_PATH_ANOMALOUS = OUT_DIR / "anomalous_protsyn_nscore_zero.tsv"


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


def create_species_mapper():
    """
    Create a mapping from various species name formats to a canonical form.
    """
    species_map = {
        # Dog variants
        "canis lupus": "dog",
        "canis lupus familiaris": "dog",
        "canis_lupus": "dog",
        "canis": "dog",
        
        # Mouse variants
        "mus musculus": "mouse",
        "mus": "mouse",
        
        # Rat variants
        "rattus norvegicus": "rat",
        "rattus": "rat",
        
        # Macaque variants
        "macaca fascicularis": "macaque",
        "macaca": "macaque",
    }
    return species_map


def normalize_species_robust(species_str):
    """
    Normalize species names to canonical form, handling multiple formats.
    
    Args:
        species_str: Species name in any format
    
    Returns:
        Canonical species name (e.g., "dog", "mouse")
    """
    if pd.isna(species_str):
        return ""
    
    # Convert to lowercase and replace underscores with spaces
    normalized = species_str.lower().replace("_", " ").strip()
    
    # Try to map to canonical form
    species_map = create_species_mapper()
    
    # Try exact match first
    if normalized in species_map:
        return species_map[normalized]
    
    # Try first two words
    parts = normalized.split()
    if len(parts) >= 2:
        two_words = " ".join(parts[:2])
        if two_words in species_map:
            return species_map[two_words]
    
    # Try first word only (genus)
    if len(parts) >= 1:
        first_word = parts[0]
        if first_word in species_map:
            return species_map[first_word]
    
    # If no match, return the normalized form
    return normalized

def build_protsyn_index(protsyn_df, gene_to_species):
    """
    Build an index for fast prot-syn lookup.
    Structure: {(gene, target_species_canonical): partner_gene}
    
    This pre-computes all possible lookups so we don't search repeatedly.
    """
    print("\nBuilding prot-syn index for fast lookup...")
    
    index = {}
    
    for _, row in protsyn_df.iterrows():
        g1 = row["Gene1"]
        g2 = row["Gene2"]
        sp1_norm = normalize_species_robust(row.get("Species1", ""))
        sp2_norm = normalize_species_robust(row.get("Species2", ""))
        
        # For each gene, store which partner it has in which species
        # Gene1 -> Gene2 (species1 -> species2)
        key1 = (g1, sp2_norm)
        if key1 not in index:
            index[key1] = g2
        
        # Gene2 -> Gene1 (species2 -> species1)
        key2 = (g2, sp1_norm)
        if key2 not in index:
            index[key2] = g1
    
    print(f"✓ Built index with {len(index)} gene-species pairs")
    return index


def find_protsyn_partner_fast(gene, target_species, protsyn_index, gene_to_species):
    """
    Find prot-syn partner using pre-built index (FAST).
    
    Args:
        gene: Gene ID to search for
        target_species: Species we're looking for a partner in
        protsyn_index: Pre-built index {(gene, target_species): partner}
        gene_to_species: Dict mapping gene -> species
    
    Returns:
        partner_gene (str or None)
    """
    gene_species = gene_to_species.get(gene)
    if gene_species is None:
        return None
    
    # Normalize target species
    target_species_norm = normalize_species_robust(target_species)
    
    # Lookup in index
    key = (gene, target_species_norm)
    return protsyn_index.get(key)



# =============================================================================
# Data Loading Functions
# =============================================================================

def load_data():
    """
    Load all required data files.
    
    Returns:
        Tuple of (nbh_dict, geomean_dict, map_df, protsyn_df, fn_df)
    """
    print("="*60)
    print("STEP 1: Loading Data")
    print("="*60)
    
    # Load dictionaries
    with open(NEIGHBORHOOD_DICT_PKL, "rb") as f:
        nbh_dict = pickle.load(f)
    
    with open(GEOMEAN_DICT_PKL, "rb") as f:
        geomean_dict = pickle.load(f)
    
    # Load tables
    map_df = pd.read_csv(MAP_PROTEIN_GENES_PATH, sep="\t")
    map_df["GeneID"] = map_df["GeneID"].astype(str)
    map_df["ProteinID"] = map_df["ProteinID"].astype(str)
    
    protsyn_df = pd.read_csv(PROTSYN_ORTH_PATH, sep="\t", dtype=str)
    fn_df = pd.read_csv(FALSE_NEGATIVES_PATH, sep="\t")
    
    print(f"✓ Loaded {len(map_df)} protein-gene mappings")
    print(f"✓ Loaded {len(protsyn_df)} prot-syn ortholog pairs")
    print(f"✓ Loaded {len(fn_df)} false negatives")
    
    return nbh_dict, geomean_dict, map_df, protsyn_df, fn_df


def build_mappings(map_df, geomean_dict):
    """
    Build mapping dictionaries from loaded data.
    
    Args:
        map_df: Mapping dataframe with GeneID, ProteinID, Species
        geomean_dict: Dictionary of gene relationships (ALREADY in correct format)
    
    Returns:
        Tuple of (prot2gene, gene2prot, gene_to_species, relationship_dict)
    """
    print("\nBuilding mappings...")
    
    # Protein <-> Gene mapping
    prot2gene = dict(zip(map_df["ProteinID"], map_df["GeneID"]))
    gene2prot = dict(zip(map_df["GeneID"], map_df["ProteinID"]))
    
    # Gene -> Species mapping
    gene_to_species = dict(zip(map_df["GeneID"], map_df["Species"]))
    
    # Use geomean_dict directly as relationship_dict (just ensure string keys)
    # The dict is already in format: {gene: [(related_gene, score), ...]}
    relationship_dict = {str(gene): related_list for gene, related_list in geomean_dict.items()}
    
    print(f"✓ Using relationship dict with {len(relationship_dict)} genes (no rebuilding needed)")
    
    return prot2gene, gene2prot, gene_to_species, relationship_dict



def prepare_protsyn_data(protsyn_df, gene_to_species):
    """
    Prepare prot-syn dataframe by adding species info if missing.
    
    Args:
        protsyn_df: Prot-syn ortholog dataframe
        gene_to_species: Dict mapping gene -> species
    
    Returns:
        Cleaned prot-syn dataframe
    """
    # Add species info if missing
    if "Species1" not in protsyn_df.columns or "Species2" not in protsyn_df.columns:
        print("Adding species info to prot-syn...")
        protsyn_df["Species1"] = protsyn_df["Gene1"].map(gene_to_species)
        protsyn_df["Species2"] = protsyn_df["Gene2"].map(gene_to_species)
    
    # Clean data
    protsyn_df = protsyn_df.dropna(subset=["Gene1", "Gene2"])
    protsyn_df["Gene1"] = protsyn_df["Gene1"].astype(str)
    protsyn_df["Gene2"] = protsyn_df["Gene2"].astype(str)
    
    return protsyn_df


# =============================================================================
# Classification Functions
# =============================================================================

def classify_fn_pairs_fast(fn_df, protsyn_index, gene_to_species, gene2prot):
    """
    Classify false negative pairs using fast index lookup.
    
    Args:
        fn_df: False negatives dataframe
        protsyn_index: Pre-built prot-syn index
        gene_to_species: Dict mapping gene -> species
        gene2prot: Dict mapping gene -> protein
    
    Returns:
        DataFrame with classified FN pairs
    """
    print("\n" + "="*60)
    print("STEP 2: Classifying False Negatives (OPTIMIZED)")
    print("="*60)
    
    # Validate columns
    required_cols = {"Protein1", "Protein2", "Gene1", "Gene2", "Species1", "Species2"}
    if not required_cols.issubset(fn_df.columns):
        raise ValueError(f"FN file missing required columns. Has: {fn_df.columns.tolist()}")
    
    # Ensure correct types
    fn_df = fn_df.copy()
    fn_df["Protein1"] = fn_df["Protein1"].astype(str)
    fn_df["Protein2"] = fn_df["Protein2"].astype(str)
    fn_df["Gene1"] = fn_df["Gene1"].astype(str)
    fn_df["Gene2"] = fn_df["Gene2"].astype(str)
    
    print(f"\nClassifying {len(fn_df)} FN pairs...")
    
    # Vectorize the lookup (much faster than iterrows)
    results = []
    
    for idx in range(len(fn_df)):
        row = fn_df.iloc[idx]
        
        g1 = row["Gene1"]
        g2 = row["Gene2"]
        p1 = row["Protein1"]
        p2 = row["Protein2"]
        sp1 = row["Species1"]
        sp2 = row["Species2"]
        
        # Fast index lookup
        partner_g1 = find_protsyn_partner_fast(g1, sp2, protsyn_index, gene_to_species)
        partner_g2 = find_protsyn_partner_fast(g2, sp1, protsyn_index, gene_to_species)
        
        # Determine status
        if partner_g1 is None:
            status_g1 = "no_result"
        elif partner_g1 == g2:
            status_g1 = "agreement"
        else:
            status_g1 = "contradictory"
        
        if partner_g2 is None:
            status_g2 = "no_result"
        elif partner_g2 == g1:
            status_g2 = "agreement"
        else:
            status_g2 = "contradictory"
        
        # Overall category
        if status_g1 == "no_result" and status_g2 == "no_result":
            category = "both_no_result"
        elif status_g1 == "contradictory" and status_g2 == "contradictory":
            category = "both_contradictory"
        elif status_g1 == "contradictory":
            category = "gene1_only"
        elif status_g2 == "contradictory":
            category = "gene2_only"
        else:
            category = "other"
        
        # Get protein IDs
        partner_p1 = gene2prot.get(partner_g1) if partner_g1 else None
        partner_p2 = gene2prot.get(partner_g2) if partner_g2 else None
        
        results.append({
            "Pair": row.get("Pair", ""),
            "Protein1": p1,
            "Protein2": p2,
            "Gene1": g1,
            "Gene2": g2,
            "Species1": sp1,
            "Species2": sp2,
            "PS_Partner_Gene1": partner_g1 if partner_g1 else "",
            "PS_Partner_Prot1": partner_p1 if partner_p1 else "",
            "PS_Partner_Gene2": partner_g2 if partner_g2 else "",
            "PS_Partner_Prot2": partner_p2 if partner_p2 else "",
            "Gene1_Status": status_g1,
            "Gene2_Status": status_g2,
            "FN_Category": category,
        })
        
        if (idx + 1) % 1000 == 0:
            print(f"  Processed {idx + 1}/{len(fn_df)} rows...")
    
    classified_df = pd.DataFrame(results)
    
    # Save
    classified_df.to_csv(FN_CLASSIFIED_OUT, sep="\t", index=False)
    print(f"\n✓ Saved classified FN: {FN_CLASSIFIED_OUT}")
    
    print("\nClassification summary:")
    print(classified_df["FN_Category"].value_counts())
    
    if "Pair" in classified_df.columns and classified_df["Pair"].notna().any():
        print("\nBy species pair:")
        print(classified_df.groupby("Pair")["FN_Category"].value_counts())
    
    return classified_df

# =============================================================================
# N-score Computation Functions
# =============================================================================

def compute_nscore_for_pair(g1, g2, nbh_dict, relationship_dict):
    """
    Compute N-score for a single gene pair.
    
    Args:
        g1: Gene 1 ID
        g2: Gene 2 ID
        nbh_dict: Neighborhood dictionary
        relationship_dict: Relationship dictionary
    
    Returns:
        N-score (int or None)
    """
    neigh1 = nbh_dict.get(g1, [])
    neigh2 = nbh_dict.get(g2, [])
    
    if not neigh1 or not neigh2:
        return None
    
    score, _ = count_shared_genes(g1, g2, neigh1, neigh2, relationship_dict)
    return score


def compute_all_nscores(classified_df, nbh_dict, relationship_dict):
    """
    Compute N-scores for all pairs (Cactus and prot-syn).
    
    Args:
        classified_df: Classified FN dataframe
        nbh_dict: Neighborhood dictionary
        relationship_dict: Relationship dictionary
    
    Returns:
        DataFrame with all N-scores added
    """
    print("\n" + "="*60)
    print("STEP 3: Computing N-scores")
    print("="*60)
    
    # Compute Cactus N-scores
    print("\nComputing N-scores for Cactus pairs (Gene1-Gene2)...")
    cactus_scores = []
    missing_cactus = 0
    
    for _, row in classified_df.iterrows():
        score = compute_nscore_for_pair(row["Gene1"], row["Gene2"], nbh_dict, relationship_dict)
        if score is None:
            missing_cactus += 1
        cactus_scores.append(score)
    
    classified_df["Nscore_Cactus"] = cactus_scores
    print(f"  Could not compute for {missing_cactus} pairs")
    
    # Compute prot-syn N-scores for Gene1
    print("\nComputing N-scores for prot-syn pairs (Gene1-Partner1)...")
    ps_scores_g1 = []
    missing_ps_g1 = 0
    
    for _, row in classified_df.iterrows():
        partner_g1 = row["PS_Partner_Gene1"]
        
        if not partner_g1 or partner_g1 == "":
            ps_scores_g1.append(None)
            continue
        
        score = compute_nscore_for_pair(row["Gene1"], partner_g1, nbh_dict, relationship_dict)
        if score is None:
            missing_ps_g1 += 1
        ps_scores_g1.append(score)
    
    classified_df["Nscore_ProtSyn_Gene1"] = ps_scores_g1
    print(f"  Could not compute for {missing_ps_g1} pairs (with partners)")
    
    # Compute prot-syn N-scores for Gene2
    print("\nComputing N-scores for prot-syn pairs (Gene2-Partner2)...")
    ps_scores_g2 = []
    missing_ps_g2 = 0
    
    for _, row in classified_df.iterrows():
        partner_g2 = row["PS_Partner_Gene2"]
        
        if not partner_g2 or partner_g2 == "":
            ps_scores_g2.append(None)
            continue
        
        score = compute_nscore_for_pair(row["Gene2"], partner_g2, nbh_dict, relationship_dict)
        if score is None:
            missing_ps_g2 += 1
        ps_scores_g2.append(score)
    
    classified_df["Nscore_ProtSyn_Gene2"] = ps_scores_g2
    print(f"  Could not compute for {missing_ps_g2} pairs (with partners)")
    
    # Compute fractions
    print("\nComputing fractions...")
    fractions_g1 = []
    fractions_g2 = []

    for _, row in classified_df.iterrows():
        nscore_c = row["Nscore_Cactus"]
        nscore_ps_g1 = row["Nscore_ProtSyn_Gene1"]
        nscore_ps_g2 = row["Nscore_ProtSyn_Gene2"]
        
        # Fraction for Gene1 - UPDATED to allow Cactus N-score = 0
        if pd.notna(nscore_c) and pd.notna(nscore_ps_g1) and nscore_ps_g1 > 0:
            fractions_g1.append(nscore_c / nscore_ps_g1)
        else:
            fractions_g1.append(None)
        
        # Fraction for Gene2 - UPDATED to allow Cactus N-score = 0
        if pd.notna(nscore_c) and pd.notna(nscore_ps_g2) and nscore_ps_g2 > 0:
            fractions_g2.append(nscore_c / nscore_ps_g2)
        else:
            fractions_g2.append(None)

    classified_df["Fraction_Gene1"] = fractions_g1
    classified_df["Fraction_Gene2"] = fractions_g2

    # Save
    classified_df.to_csv(FN_WITH_SCORES_OUT, sep="\t", index=False)
    print(f"\n✓ Saved FN with all scores: {FN_WITH_SCORES_OUT}")

    # REPORT UNEXPECTED CASES: Prot-syn N-score = 0
    print("\n--- Checking for Unexpected Cases ---")

    ps_zero_g1 = classified_df[
        (classified_df["PS_Partner_Gene1"].fillna("").str.strip() != "") &
        (classified_df["Nscore_ProtSyn_Gene1"] == 0)
    ]
    ps_zero_g2 = classified_df[
        (classified_df["PS_Partner_Gene2"].fillna("").str.strip() != "") &
        (classified_df["Nscore_ProtSyn_Gene2"] == 0)
    ]

    zero_sides = len(ps_zero_g1) + len(ps_zero_g2)

    # Unique FN pairs where at least one side has prot-syn Nscore==0
    pairs_zero_any_side = set(zip(ps_zero_g1["Gene1"].astype(str), ps_zero_g1["Gene2"].astype(str))) | \
                        set(zip(ps_zero_g2["Gene1"].astype(str), ps_zero_g2["Gene2"].astype(str)))

    zero_pairs = len(pairs_zero_any_side)

    if zero_sides > 0:
        print("⚠ WARNING: Found prot-syn partner(s) with N-score = 0")
        print(f"  Zero-sides (Gene1-side + Gene2-side): {zero_sides}")
        print(f"    - Gene1-side zero: {len(ps_zero_g1)}")
        print(f"    - Gene2-side zero: {len(ps_zero_g2)}")
        print(f"  Unique FN pairs affected (>=1 side zero): {zero_pairs}")
    else:
        print("✓ All prot-syn pairs have N-score > 0 (as expected)")


    return classified_df

# =============================================================================
# Analysis Functions
# =============================================================================

def find_cactus_better_cases(classified_df):
    """
    Find cases where Cactus picked higher N-score than prot-syn.
    
    Args:
        classified_df: Classified FN dataframe with N-scores
    
    Returns:
        DataFrame with cases where Cactus is better
    """
    print("\n" + "="*60)
    print("STEP 4: Finding Cases Where Cactus Picked Better N-score")
    print("="*60)
    
    better_cases = []

    gene1_better_sides = 0
    gene2_better_sides = 0
    gene1_total_sides = 0
    gene2_total_sides = 0

    # Initialize to avoid UnboundLocalError when there are 0 contradictory cases
    gene1_better = pd.DataFrame()
    gene2_better = pd.DataFrame()

    
    # Check Gene1 contradictory cases
    gene1_contradictory = classified_df[
        (classified_df["Gene1_Status"] == "contradictory") &
        (classified_df["Nscore_Cactus"].notna()) &
        (classified_df["Nscore_ProtSyn_Gene1"].notna()) &
        (classified_df["Nscore_Cactus"] > 0) &
        (classified_df["Nscore_ProtSyn_Gene1"] > 0)
    ].copy()
    
    if len(gene1_contradictory) > 0:
        gene1_better = gene1_contradictory[gene1_contradictory["Fraction_Gene1"] > 1]
        pct = (len(gene1_better) / len(gene1_contradictory) * 100)
        print(f"\nGene1-side: Cactus better in {len(gene1_better)}/{len(gene1_contradictory)} sides ({pct:.1f}%)")
        
        for _, row in gene1_better.iterrows():
            better_cases.append({
                **row.to_dict(),
                "Better_For": "gene1",
                "Relevant_Fraction": row["Fraction_Gene1"],
                "Relevant_PS_Nscore": row["Nscore_ProtSyn_Gene1"],
                "Relevant_PS_Partner_Gene": row["PS_Partner_Gene1"],
                "Relevant_PS_Partner_Prot": row["PS_Partner_Prot1"],
            })
    else:
        print(f"\nGene1: No contradictory cases found (0 cases)")

    gene1_total_sides = len(gene1_contradictory)
    gene1_better_sides = len(gene1_better)

    
    # Check Gene2 contradictory cases
    gene2_contradictory = classified_df[
        (classified_df["Gene2_Status"] == "contradictory") &
        (classified_df["Nscore_Cactus"].notna()) &
        (classified_df["Nscore_ProtSyn_Gene2"].notna()) &
        (classified_df["Nscore_Cactus"] > 0) &
        (classified_df["Nscore_ProtSyn_Gene2"] > 0)
    ].copy()
    
    if len(gene2_contradictory) > 0:
        gene2_better = gene2_contradictory[gene2_contradictory["Fraction_Gene2"] > 1]
        pct = (len(gene2_better) / len(gene2_contradictory) * 100)
        print(f"Gene2-side: Cactus better in {len(gene2_better)}/{len(gene2_contradictory)} sides ({pct:.1f}%)")
        
        for _, row in gene2_better.iterrows():
            # Check if this case was already added (both genes better)
            existing = [c for c in better_cases if str(c["Gene1"]) == str(row["Gene1"]) and str(c["Gene2"]) == str(row["Gene2"])]
            
            if existing:
                # Update existing entry to mark as "both"
                existing[0]["Better_For"] = "both"
                existing[0]["Gene2_Fraction"] = row["Fraction_Gene2"]
                existing[0]["Gene2_PS_Nscore"] = row["Nscore_ProtSyn_Gene2"]
            else:
                # Add new entry for gene2 only
                better_cases.append({
                    **row.to_dict(),
                    "Better_For": "gene2",
                    "Relevant_Fraction": row["Fraction_Gene2"],
                    "Relevant_PS_Nscore": row["Nscore_ProtSyn_Gene2"],
                    "Relevant_PS_Partner_Gene": row["PS_Partner_Gene2"],
                    "Relevant_PS_Partner_Prot": row["PS_Partner_Prot2"],
                })
    else:
        print(f"Gene2: No contradictory cases found (0 cases)")

    gene2_total_sides = len(gene2_contradictory)
    gene2_better_sides = len(gene2_better)


    if better_cases:
        better_df = pd.DataFrame(better_cases)
        
        # Sort by Better_For and then by Relevant_Fraction (descending)
        better_df = better_df.sort_values(["Better_For", "Relevant_Fraction"], ascending=[True, False])
        
        better_df.to_csv(OUT_PATH_CACTUS_BETTER, sep="\t", index=False)
        
        print(f"\n✓ Saved cases where Cactus picked better: {OUT_PATH_CACTUS_BETTER}")
        print(f"  Total cases: {len(better_df)}")
        print(f"  Breakdown:")
        print(f"    - Gene1 only: {(better_df['Better_For'] == 'gene1').sum()}")
        print(f"    - Gene2 only: {(better_df['Better_For'] == 'gene2').sum()}")
        print(f"    - Both genes: {(better_df['Better_For'] == 'both').sum()}")
        
        if "Pair" in better_df.columns:
            print(f"\n  By species pair:")
            print(better_df.groupby("Pair")["Better_For"].value_counts())
        
        better_sides_total = gene1_better_sides + gene2_better_sides

        print("\n--- Clarification of units ---")
        print(f"  Better-sides total (Gene1-side + Gene2-side): {better_sides_total}")
        print(f"  Better-pairs saved to file (unique FN pairs): {len(better_df)}")
        print(f"  (Pairs can contribute 1 or 2 better-sides, hence totals differ.)")

        return better_df
    else:
        print("\n✓ No cases found where Cactus picked better N-score")
        print("   All contradictory cases show prot-syn with equal or better N-scores.")
        return pd.DataFrame()





def save_protsyn_zero_nscore_cases(classified_df):
    """
    Save cases where prot-syn picked a partner with N-score = 0.
    Option A: one row per FN pair, with BOTH prot-syn alternatives (if present).
    """
    print("\n" + "="*60)
    print("STEP 4.5: Finding Anomalous Cases (Prot-syn N-score = 0)")
    print("="*60)

    # Build a dict keyed by (Gene1, Gene2) to merge gene1-side and gene2-side anomalies
    anomalous = {}

    # --- Gene1-side anomalies ---
    ps_zero_g1 = classified_df[
        (classified_df["PS_Partner_Gene1"].fillna("").str.strip() != "") &
        (classified_df["Nscore_ProtSyn_Gene1"].notna()) &
        (classified_df["Nscore_ProtSyn_Gene1"] == 0)
    ].copy()

    print(f"\nGene1-side: Found {len(ps_zero_g1)} cases with prot-syn N-score = 0")

    for _, row in ps_zero_g1.iterrows():
        key = (str(row["Gene1"]), str(row["Gene2"]))

        if key not in anomalous:
            anomalous[key] = {
                "Pair": row.get("Pair", ""),
                "Protein1": row["Protein1"],
                "Protein2": row["Protein2"],
                "Gene1": row["Gene1"],
                "Gene2": row["Gene2"],
                "Species1": row["Species1"],
                "Species2": row["Species2"],
                "Nscore_Cactus": row["Nscore_Cactus"],
                # initialize both sides as empty
                "ProtSyn_Partner_Gene1": "",
                "ProtSyn_Partner_Protein1": "",
                "Nscore_ProtSyn_Gene1": np.nan,
                "ProtSyn_Partner_Gene2": "",
                "ProtSyn_Partner_Protein2": "",
                "Nscore_ProtSyn_Gene2": np.nan,
            }

        anomalous[key]["ProtSyn_Partner_Gene1"] = row["PS_Partner_Gene1"]
        anomalous[key]["ProtSyn_Partner_Protein1"] = row["PS_Partner_Prot1"]
        anomalous[key]["Nscore_ProtSyn_Gene1"] = row["Nscore_ProtSyn_Gene1"]

    # --- Gene2-side anomalies ---
    ps_zero_g2 = classified_df[
        (classified_df["PS_Partner_Gene2"].fillna("").str.strip() != "") &
        (classified_df["Nscore_ProtSyn_Gene2"].notna()) &
        (classified_df["Nscore_ProtSyn_Gene2"] == 0)
    ].copy()

    print(f"Gene2-side: Found {len(ps_zero_g2)} cases with prot-syn N-score = 0")

    for _, row in ps_zero_g2.iterrows():
        key = (str(row["Gene1"]), str(row["Gene2"]))

        if key not in anomalous:
            anomalous[key] = {
                "Pair": row.get("Pair", ""),
                "Protein1": row["Protein1"],
                "Protein2": row["Protein2"],
                "Gene1": row["Gene1"],
                "Gene2": row["Gene2"],
                "Species1": row["Species1"],
                "Species2": row["Species2"],
                "Nscore_Cactus": row["Nscore_Cactus"],
                # initialize both sides as empty
                "ProtSyn_Partner_Gene1": "",
                "ProtSyn_Partner_Protein1": "",
                "Nscore_ProtSyn_Gene1": np.nan,
                "ProtSyn_Partner_Gene2": "",
                "ProtSyn_Partner_Protein2": "",
                "Nscore_ProtSyn_Gene2": np.nan,
            }

        anomalous[key]["ProtSyn_Partner_Gene2"] = row["PS_Partner_Gene2"]
        anomalous[key]["ProtSyn_Partner_Protein2"] = row["PS_Partner_Prot2"]
        anomalous[key]["Nscore_ProtSyn_Gene2"] = row["Nscore_ProtSyn_Gene2"]

    if not anomalous:
        print("\n✓ No anomalous cases found (all prot-syn N-scores > 0, as expected)")
        return pd.DataFrame()

    anomalous_df = pd.DataFrame(list(anomalous.values()))

    # Determine anomaly type based on which side(s) are present
    has_g1 = anomalous_df["ProtSyn_Partner_Gene1"].fillna("").str.strip() != ""
    has_g2 = anomalous_df["ProtSyn_Partner_Gene2"].fillna("").str.strip() != ""

    anomalous_df["Anomaly_Type"] = np.where(
        has_g1 & has_g2, "both_ps_nscore_zero",
        np.where(has_g1, "gene1_ps_nscore_zero", "gene2_ps_nscore_zero")
    )

    # Optional: enforce a consistent, unambiguous column order (recommended)
    desired_cols = [
        "Pair",
        "Protein1", "Protein2",
        "Gene1", "Gene2",
        "Species1", "Species2",
        "Nscore_Cactus",
        "ProtSyn_Partner_Gene1", "ProtSyn_Partner_Protein1", "Nscore_ProtSyn_Gene1",
        "ProtSyn_Partner_Gene2", "ProtSyn_Partner_Protein2", "Nscore_ProtSyn_Gene2",
        "Anomaly_Type",
    ]
    # Keep only columns that exist (safe if you tweak upstream)
    anomalous_df = anomalous_df[[c for c in desired_cols if c in anomalous_df.columns]]

    # Save
    out_path = OUT_DIR / "anomalous_protsyn_nscore_zero.tsv"
    anomalous_df.to_csv(out_path, sep="\t", index=False)


    print(f"\n WARNING: Saved {len(anomalous_df)} anomalous cases: {out_path}")
    print("  Breakdown:")
    print(f"    - Gene1 only: {(anomalous_df['Anomaly_Type'] == 'gene1_ps_nscore_zero').sum()}")
    print(f"    - Gene2 only: {(anomalous_df['Anomaly_Type'] == 'gene2_ps_nscore_zero').sum()}")
    print(f"    - Both genes: {(anomalous_df['Anomaly_Type'] == 'both_ps_nscore_zero').sum()}")

    gene1_only = (anomalous_df['Anomaly_Type'] == 'gene1_ps_nscore_zero').sum()
    gene2_only = (anomalous_df['Anomaly_Type'] == 'gene2_ps_nscore_zero').sum()
    both = (anomalous_df['Anomaly_Type'] == 'both_ps_nscore_zero').sum()

    # This converts "pair rows" back into "side events"
    zero_sides_implied = gene1_only + gene2_only + 2 * both

    print(f"  Note: File is ONE ROW PER FN PAIR.")
    print(f"        Pair rows = {len(anomalous_df)}; implied zero-sides = {zero_sides_implied}")


    return anomalous_df


# =============================================================================
# Plotting Functions
# =============================================================================

def plot_distributions(classified_df):
    """
    Create 2x2 grid plot of N-score distributions.
    
    Args:
        classified_df: Classified FN dataframe with N-scores
    """
    print("\n" + "="*60)
    print("STEP 5: Plotting Distributions (2x2 Grid)")
    print("="*60)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('False Negatives N-score Distributions', fontsize=16, fontweight='bold')
    
    # ---- Plot 1: All FN (Cactus N-scores) ----
    ax1 = axes[0, 0]
    all_cactus = classified_df["Nscore_Cactus"].dropna()
    
    if not all_cactus.empty:
        max_score = int(all_cactus.max())
        bins = np.arange(0, max_score + 2, 1)
        
        ax1.hist(all_cactus, bins=bins, edgecolor='black', color='steelblue', alpha=0.7)
        ax1.set_xlabel('N-score (Cactus)', fontsize=11)
        ax1.set_ylabel('Number of FN pairs', fontsize=11)
        ax1.set_title('All False Negatives\n(Cactus pairs)', fontsize=12, fontweight='bold')
        ax1.grid(axis='y', alpha=0.3, linestyle='--')
        # NOTE: removed n= label to avoid confusion (pairs vs sides)
    
    # ---- Plot 2: No Result Cases ----
    ax2 = axes[0, 1]
    no_result = classified_df[classified_df["FN_Category"] == "both_no_result"]["Nscore_Cactus"].dropna()
    
    if not no_result.empty:
        max_score = int(no_result.max())
        bins = np.arange(0, max_score + 2, 1)
        
        ax2.hist(no_result, bins=bins, edgecolor='black', color='coral', alpha=0.7)
        ax2.set_xlabel('N-score (Cactus)', fontsize=11)
        ax2.set_ylabel('Number of FN pairs', fontsize=11)
        ax2.set_title('No Result Cases\n(neither gene in prot-syn)', fontsize=12, fontweight='bold')
        ax2.grid(axis='y', alpha=0.3, linestyle='--')
        # NOTE: removed n= label to avoid confusion (pairs vs sides)
    else:
        ax2.text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14)
        ax2.set_title('No Result Cases', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Number of FN pairs', fontsize=11)
    
    # ---- Plot 3: Contradictory Cases (Prot-syn N-scores) ----
    ax3 = axes[1, 0]
    contradictory_mask = classified_df["FN_Category"].isin([
        "both_contradictory", "gene1_only", "gene2_only"
    ])
    contradictory = classified_df[contradictory_mask]
    
    # Combine both prot-syn N-scores (NOTE: this counts SIDES, not pairs)
    ps_combined = pd.concat([
        contradictory["Nscore_ProtSyn_Gene1"].dropna(),
        contradictory["Nscore_ProtSyn_Gene2"].dropna()
    ])
    
    if not ps_combined.empty:
        max_score = int(ps_combined.max())
        bins = np.arange(0, max_score + 2, 1)
        
        ax3.hist(ps_combined, bins=bins, edgecolor='black', color='mediumseagreen', alpha=0.7)
        ax3.set_xlabel('N-score (Prot-syn)', fontsize=11)
        ax3.set_ylabel('Number of prot-syn sides', fontsize=11)
        ax3.set_title('Contradictory Cases\n(prot-syn partners)', fontsize=12, fontweight='bold')
        ax3.grid(axis='y', alpha=0.3, linestyle='--')
        # NOTE: removed n= label to avoid confusion (pairs vs sides)
    else:
        ax3.text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14)
        ax3.set_title('Contradictory Cases', fontsize=12, fontweight='bold')
        ax3.set_ylabel('Number of prot-syn sides', fontsize=11)
    
    # ---- Plot 4: Fraction Distribution ----
    ax4 = axes[1, 1]
    fractions_combined = pd.concat([
        contradictory["Fraction_Gene1"].dropna(),
        contradictory["Fraction_Gene2"].dropna()
    ])

    if not fractions_combined.empty:
        # Use more bins for fraction since it's continuous
        bins = np.linspace(0, fractions_combined.max(), 30)
        
        ax4.hist(fractions_combined, bins=bins, edgecolor='black', color='mediumpurple', alpha=0.7)
        ax4.axvline(x=1.0, color='red', linestyle='--', linewidth=2.5, 
                    label='Equal (Fraction=1)', zorder=10)
        ax4.set_xlabel('Fraction (Cactus / Prot-syn)', fontsize=11)
        ax4.set_ylabel('Number of prot-syn sides', fontsize=11)
        ax4.set_title('Fraction Distribution\n(contradictory cases)', fontsize=12, fontweight='bold')
        
        # Move legend to middle right - away from both the peak and n= text
        ax4.legend(loc='center right', fontsize=10, framealpha=0.95)
        
        ax4.grid(axis='y', alpha=0.3, linestyle='--')
        # NOTE: removed n= label to avoid confusion (pairs vs sides)
    else:
        ax4.text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14)
        ax4.set_title('Fraction Distribution', fontsize=12, fontweight='bold')
        ax4.set_ylabel('Number of prot-syn sides', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(PLOT_DISTRIBUTIONS, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Saved 2x2 distribution plot: {PLOT_DISTRIBUTIONS}")



def print_summary_statistics(classified_df):
    """
    Print summary statistics for all distributions.
    
    Args:
        classified_df: Classified FN dataframe with N-scores
    """
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    
    # All Cactus N-scores
    all_cactus = classified_df["Nscore_Cactus"].dropna()
    print(f"\n1. All FN (Cactus N-scores):")
    if not all_cactus.empty:
        print(f"   Count (pairs):     {len(all_cactus)}")
        print(f"   Mean:      {all_cactus.mean():.2f}")
        print(f"   Median:    {all_cactus.median():.2f}")
        print(f"   N-score=0: {(all_cactus == 0).sum()} ({(all_cactus == 0).sum()/len(all_cactus)*100:.1f}%)")
    
    # No Result cases
    no_result = classified_df[classified_df["FN_Category"] == "both_no_result"]["Nscore_Cactus"].dropna()
    print(f"\n2. No Result Cases (Cactus N-scores):")
    if not no_result.empty:
        print(f"   Count (pairs):     {len(no_result)}")
        print(f"   Mean:      {no_result.mean():.2f}")
        print(f"   Median:    {no_result.median():.2f}")
    
    # Contradictory cases - Prot-syn N-scores
    contradictory_mask = classified_df["FN_Category"].isin([
        "both_contradictory", "gene1_only", "gene2_only"
    ])
    contradictory = classified_df[contradictory_mask]
    ps_combined = pd.concat([
        contradictory["Nscore_ProtSyn_Gene1"].dropna(),
        contradictory["Nscore_ProtSyn_Gene2"].dropna()
    ])
    
    print(f"\n3. Contradictory Cases (Prot-syn N-scores):")
    if not ps_combined.empty:
        print(f"   Count (sides):     {len(ps_combined)}")
        print(f"   Mean:      {ps_combined.mean():.2f}")
        print(f"   Median:    {ps_combined.median():.2f}")
        print(f"   N-score=0: {(ps_combined == 0).sum()} ({(ps_combined == 0).sum()/len(ps_combined)*100:.1f}%)")
    
    # Fraction distribution
    fractions_combined = pd.concat([
        contradictory["Fraction_Gene1"].dropna(),
        contradictory["Fraction_Gene2"].dropna()
    ])
    
    print(f"\n4. Fraction Distribution:")
    if not fractions_combined.empty:
        print(f"   Count (sides):              {len(fractions_combined)}")
        print(f"   Mean:               {fractions_combined.mean():.2f}")
        print(f"   Median:             {fractions_combined.median():.2f}")
        print(f"   Cactus better (>1):  {(fractions_combined > 1).sum()} ({(fractions_combined > 1).sum()/len(fractions_combined)*100:.1f}%)")
        print(f"   ProtSyn better (<1): {(fractions_combined < 1).sum()} ({(fractions_combined < 1).sum()/len(fractions_combined)*100:.1f}%)")
        print(f"   Equal (=1):          {(fractions_combined == 1).sum()} ({(fractions_combined == 1).sum()/len(fractions_combined)*100:.1f}%)")


# =============================================================================
# Main Workflow
# =============================================================================
def main():
    """
    Main workflow for false negatives analysis (OPTIMIZED).
    """
    # Step 1: Load data
    nbh_dict, geomean_dict, map_df, protsyn_df, fn_df = load_data()
    
    # Step 2: Build mappings
    prot2gene, gene2prot, gene_to_species, relationship_dict = build_mappings(map_df, geomean_dict)
    
    # Step 3: Prepare prot-syn data
    protsyn_df = prepare_protsyn_data(protsyn_df, gene_to_species)
    
    # Step 3.5: BUILD INDEX
    protsyn_index = build_protsyn_index(protsyn_df, gene_to_species)
    
    # Step 4: Classify FN pairs
    classified_df = classify_fn_pairs_fast(fn_df, protsyn_index, gene_to_species, gene2prot)
    
    # Step 5: Compute N-scores (now includes reporting of prot-syn N-score = 0)
    classified_df = compute_all_nscores(classified_df, nbh_dict, relationship_dict)
    
    # Step 6: Save anomalous cases (prot-syn N-score = 0) - NEW!
    anomalous_df = save_protsyn_zero_nscore_cases(classified_df)
    
    # Step 7: Find cases where Cactus is better
    better_df = find_cactus_better_cases(classified_df)
    
    # Step 8: Plot distributions
    plot_distributions(classified_df)
    
    # Step 9: Print summary statistics
    print_summary_statistics(classified_df)
    
    print("\n" + "="*60)
    print("✓ Analysis Complete!")
    print("="*60)


if __name__ == "__main__":
    main()

