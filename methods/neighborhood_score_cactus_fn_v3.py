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
MAP_PROTEIN_GENES_PATH = "material/sex_experiment/gene_to_protein_map_FIXED.tsv"
PROTSYN_ORTH_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_sum.tsv"
FALSE_NEGATIVES_PATH = "results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/all_false_negatives_sum.tsv"

# Output directory
OUT_DIR = Path("results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_negatives")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Output files
FN_WITH_IDS_OUT = OUT_DIR / "false_negatives_with_ids.tsv"
FN_CLASSIFIED_OUT = OUT_DIR / "false_negatives_classified.tsv"
OUT_PATH_CACTUS_BETTER = OUT_DIR / "contradictory_cactus_higher_nscore_than_protsyn.tsv"
OUT_PATH_ANOMALOUS = OUT_DIR / "anomalous_protsyn_nscore_zero.tsv"
PLOT_DISTRIBUTIONS = OUT_DIR / "fn_distributions_2x2.png"


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
    Build a SYMMETRIC index for fast prot-syn lookup.
    Structure: {(gene, target_species_canonical): partner_gene}
    
    This pre-computes all possible lookups so we don't search repeatedly.
    Both directions are indexed, making it order-independent.
    """
    print("\nBuilding symmetric prot-syn index for fast lookup...")
    
    index = {}
    
    for _, row in protsyn_df.iterrows():
        g1 = row["Gene1"]
        g2 = row["Gene2"]
        sp1_norm = normalize_species_robust(row.get("Species1", ""))
        sp2_norm = normalize_species_robust(row.get("Species2", ""))
        
        # Gene1 -> Gene2 (species1 -> species2)
        key1 = (g1, sp2_norm)
        if key1 not in index:
            index[key1] = g2
        
        # Gene2 -> Gene1 (species2 -> species1) - SYMMETRIC
        key2 = (g2, sp1_norm)
        if key2 not in index:
            index[key2] = g1
    
    print(f"✓ Built symmetric index with {len(index)} gene-species pairs")
    return index


def find_protsyn_partner(gene, target_species, protsyn_index):
    """
    Find prot-syn partner using pre-built index.
    
    Args:
        gene: Gene ID to search for
        target_species: Species we're looking for a partner in (already normalized)
        protsyn_index: Pre-built index {(gene, target_species): partner}
    
    Returns:
        partner_gene (str or None)
    """
    key = (gene, target_species)
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
    print(f"✓ Loaded {len(fn_df)} false negative pairs")
    
    return nbh_dict, geomean_dict, map_df, protsyn_df, fn_df


def build_mappings(map_df, geomean_dict):
    """
    Build mapping dictionaries from loaded data.
    
    Args:
        map_df: Mapping dataframe with GeneID, ProteinID, Species
        geomean_dict: Dictionary of gene relationships
    
    Returns:
        Tuple of (prot2gene, gene2prot, gene_to_species, relationship_dict)
    """
    print("\nBuilding mappings...")
    
    # Protein <-> Gene mapping
    prot2gene = dict(zip(map_df["ProteinID"], map_df["GeneID"]))
    gene2prot = dict(zip(map_df["GeneID"], map_df["ProteinID"]))
    
    # Gene -> Species mapping
    gene_to_species = dict(zip(map_df["GeneID"], map_df["Species"]))
    
    # Use geomean_dict directly as relationship_dict
    relationship_dict = {str(gene): related_list for gene, related_list in geomean_dict.items()}
    
    print(f"✓ Built mappings for {len(gene_to_species)} genes")
    
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


def add_fn_ids(fn_df):
    """
    Add stable FN_ID to false negatives dataframe.
    
    Args:
        fn_df: False negatives dataframe
    
    Returns:
        DataFrame with FN_ID column added
    """
    print("\n" + "="*60)
    print("STEP 2: Adding FN_IDs to Pair-Level Table")
    print("="*60)
    
    # Validate columns
    required_cols = {"Protein1", "Protein2", "Gene1", "Gene2", "Species1", "Species2"}
    if not required_cols.issubset(fn_df.columns):
        raise ValueError(f"FN file missing required columns. Has: {fn_df.columns.tolist()}")
    
    fn_df = fn_df.copy()
    fn_df["FN_ID"] = range(1, len(fn_df) + 1)
    
    # Ensure correct types
    fn_df["Protein1"] = fn_df["Protein1"].astype(str)
    fn_df["Protein2"] = fn_df["Protein2"].astype(str)
    fn_df["Gene1"] = fn_df["Gene1"].astype(str)
    fn_df["Gene2"] = fn_df["Gene2"].astype(str)
    
    # Save pair-level table
    fn_df.to_csv(FN_WITH_IDS_OUT, sep="\t", index=False)
    print(f"✓ Saved pair-level table with {len(fn_df)} FN pairs: {FN_WITH_IDS_OUT}")
    
    return fn_df


# =============================================================================
# Gene-Level Classification Functions
# =============================================================================

def create_gene_level_rows(fn_df, protsyn_index, gene_to_species, gene2prot):
    """
    Convert FN pairs to gene-level analysis rows.
    Each FN pair produces 2 rows (one per gene).
    
    Args:
        fn_df: False negatives dataframe with FN_IDs
        protsyn_index: Pre-built prot-syn index
        gene_to_species: Dict mapping gene -> species
        gene2prot: Dict mapping gene -> protein
    
    Returns:
        DataFrame with gene-level classification
    """
    print("\n" + "="*60)
    print("STEP 3: Creating Gene-Level Classification Table")
    print("="*60)
    
    gene_rows = []
    
    for _, row in fn_df.iterrows():
        fn_id = row["FN_ID"]
        g1 = row["Gene1"]
        g2 = row["Gene2"]
        p1 = row["Protein1"]
        p2 = row["Protein2"]
        sp1_raw = row["Species1"]
        sp2_raw = row["Species2"]
        
        # Normalize species
        sp1 = normalize_species_robust(sp1_raw)
        sp2 = normalize_species_robust(sp2_raw)
        
        # Get pair identifier if it exists
        pair_id = row.get("Pair", "")
        
        # ----- ROW 1: Focus on Gene1 -----
        # Cactus says: Gene1 (sp1) <-> Gene2 (sp2)
        # ProtSyn lookup: Gene1 in sp2 (target species)
        ps_partner_g1 = find_protsyn_partner(g1, sp2, protsyn_index)
        
        if ps_partner_g1 is None:
            fn_class_g1 = "no_result"
        elif ps_partner_g1 == g2:
            fn_class_g1 = "agreement"  # Won't be saved, but good for debugging
        else:
            fn_class_g1 = "contradictory"
        
        gene_rows.append({
            "FN_ID": fn_id,
            "Pair": pair_id,
            "Focus_Gene": g1,
            "Focus_Protein": p1,
            "Focus_Species": sp1_raw,
            "FN_Gene_Position": "Gene1",
            "Cactus_Ortholog_Gene": g2,
            "Cactus_Ortholog_Protein": p2,
            "Cactus_Ortholog_Species": sp2_raw,
            "ProtSyn_Ortholog_Gene": ps_partner_g1 if ps_partner_g1 else "",
            "ProtSyn_Ortholog_Protein": gene2prot.get(ps_partner_g1, "") if ps_partner_g1 else "",
            "FN_Class": fn_class_g1,
        })
        
        # ----- ROW 2: Focus on Gene2 -----
        # Cactus says: Gene2 (sp2) <-> Gene1 (sp1)
        # ProtSyn lookup: Gene2 in sp1 (target species)
        ps_partner_g2 = find_protsyn_partner(g2, sp1, protsyn_index)
        
        if ps_partner_g2 is None:
            fn_class_g2 = "no_result"
        elif ps_partner_g2 == g1:
            fn_class_g2 = "agreement"
        else:
            fn_class_g2 = "contradictory"
        
        gene_rows.append({
            "FN_ID": fn_id,
            "Pair": pair_id,
            "Focus_Gene": g2,
            "Focus_Protein": p2,
            "Focus_Species": sp2_raw,
            "FN_Gene_Position": "Gene2",
            "Cactus_Ortholog_Gene": g1,
            "Cactus_Ortholog_Protein": p1,
            "Cactus_Ortholog_Species": sp1_raw,
            "ProtSyn_Ortholog_Gene": ps_partner_g2 if ps_partner_g2 else "",
            "ProtSyn_Ortholog_Protein": gene2prot.get(ps_partner_g2, "") if ps_partner_g2 else "",
            "FN_Class": fn_class_g2,
        })
    
    gene_df = pd.DataFrame(gene_rows)
    
    print(f"✓ Created {len(gene_df)} gene-level rows from {len(fn_df)} FN pairs")
    print(f"\nClassification summary:")
    print(gene_df["FN_Class"].value_counts())
    
    if "Pair" in gene_df.columns and gene_df["Pair"].notna().any():
        print(f"\nBy species pair:")
        print(gene_df.groupby("Pair")["FN_Class"].value_counts())
    
    return gene_df


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


def compute_nscores(gene_df, nbh_dict, relationship_dict):
    """
    Compute N-scores for all gene-level rows.
    
    Args:
        gene_df: Gene-level classification dataframe
        nbh_dict: Neighborhood dictionary
        relationship_dict: Relationship dictionary
    
    Returns:
        DataFrame with N-scores added
    """
    print("\n" + "="*60)
    print("STEP 4: Computing N-scores")
    print("="*60)
    
    # Compute Cactus N-scores
    print("\nComputing N-scores for Cactus pairs...")
    cactus_scores = []
    missing_cactus = 0
    
    for _, row in gene_df.iterrows():
        score = compute_nscore_for_pair(
            row["Focus_Gene"], 
            row["Cactus_Ortholog_Gene"], 
            nbh_dict, 
            relationship_dict
        )
        if score is None:
            missing_cactus += 1
        cactus_scores.append(score)
    
    gene_df["Nscore_Cactus"] = cactus_scores
    print(f"  Could not compute for {missing_cactus} gene rows")
    
    # Compute ProtSyn N-scores
    print("\nComputing N-scores for ProtSyn pairs...")
    ps_scores = []
    missing_ps = 0
    
    for _, row in gene_df.iterrows():
        ps_partner = row["ProtSyn_Ortholog_Gene"]
        
        if not ps_partner or ps_partner == "":
            ps_scores.append(None)
            continue
        
        score = compute_nscore_for_pair(
            row["Focus_Gene"], 
            ps_partner, 
            nbh_dict, 
            relationship_dict
        )
        if score is None:
            missing_ps += 1
        ps_scores.append(score)
    
    gene_df["Nscore_ProtSyn"] = ps_scores
    print(f"  Could not compute for {missing_ps} gene rows (with ProtSyn partners)")
    
    # Compute fractions
    print("\nComputing score fractions...")
    fractions = []
    
    for _, row in gene_df.iterrows():
        nscore_c = row["Nscore_Cactus"]
        nscore_ps = row["Nscore_ProtSyn"]
        
        # Only compute fraction if both scores exist and ProtSyn > 0
        if pd.notna(nscore_c) and pd.notna(nscore_ps) and nscore_ps > 0:
            fractions.append(nscore_c / nscore_ps)
        else:
            fractions.append(None)
    
    gene_df["Score_Fraction"] = fractions
    
    # Save classified gene-level table
    gene_df.to_csv(FN_CLASSIFIED_OUT, sep="\t", index=False)
    print(f"\n✓ Saved gene-level classification: {FN_CLASSIFIED_OUT}")
    
    # Check for anomalous cases
    print("\n--- Checking for Unexpected Cases ---")
    ps_zero = gene_df[
        (gene_df["ProtSyn_Ortholog_Gene"].fillna("").str.strip() != "") &
        (gene_df["Nscore_ProtSyn"] == 0)
    ]
    
    if len(ps_zero) > 0:
        print(f" WARNING: Found {len(ps_zero)} gene rows with ProtSyn N-score = 0")
        unique_fn_ids = ps_zero["FN_ID"].nunique()
        print(f"  Affecting {unique_fn_ids} unique FN pairs")
    else:
        print("✓ All ProtSyn pairs have N-score > 0 (as expected)")
    
    return gene_df


# =============================================================================
# Analysis Functions
# =============================================================================

def extract_cactus_better_cases(gene_df):
    """
    Extract cases where Cactus picked a higher N-score than ProtSyn.
    
    Args:
        gene_df: Gene-level classification dataframe with N-scores
    
    Returns:
        DataFrame with cases where Cactus is better
    """
    print("\n" + "="*60)
    print("STEP 5: Extracting Cases Where Cactus Picked Better N-score")
    print("="*60)
    
    # Filter for contradictory cases where Cactus is better
    better = gene_df[
        (gene_df["FN_Class"] == "contradictory") &
        (gene_df["Nscore_Cactus"].notna()) &
        (gene_df["Nscore_ProtSyn"].notna()) &
        (gene_df["Nscore_Cactus"] > 0) &
        (gene_df["Nscore_ProtSyn"] > 0) &
        (gene_df["Score_Fraction"] > 1)
    ].copy()
    
    if len(better) > 0:
        # Sort by fraction (descending)
        better = better.sort_values("Score_Fraction", ascending=False)
        
        better.to_csv(OUT_PATH_CACTUS_BETTER, sep="\t", index=False)
        
        print(f"✓ Saved {len(better)} gene rows where Cactus picked better: {OUT_PATH_CACTUS_BETTER}")
        print(f"  Affecting {better['FN_ID'].nunique()} unique FN pairs")
        
        if "Pair" in better.columns:
            print(f"\n  By species pair:")
            print(better["Pair"].value_counts())
        
        # Calculate percentage
        contradictory_total = len(gene_df[
            (gene_df["FN_Class"] == "contradictory") &
            (gene_df["Nscore_Cactus"].notna()) &
            (gene_df["Nscore_ProtSyn"].notna()) &
            (gene_df["Nscore_Cactus"] > 0) &
            (gene_df["Nscore_ProtSyn"] > 0)
        ])
        
        if contradictory_total > 0:
            pct = (len(better) / contradictory_total * 100)
            print(f"\n  Cactus better in {len(better)}/{contradictory_total} contradictory gene rows ({pct:.1f}%)")
        
        return better
    else:
        print("✓ No cases found where Cactus picked better N-score")
        print("  All contradictory cases show ProtSyn with equal or better N-scores.")
        return pd.DataFrame()


def extract_anomalous_cases(gene_df):
    """
    Extract cases where ProtSyn picked a partner with N-score = 0.
    
    Args:
        gene_df: Gene-level classification dataframe with N-scores
    
    Returns:
        DataFrame with anomalous cases
    """
    print("\n" + "="*60)
    print("STEP 6: Extracting Anomalous Cases (ProtSyn N-score = 0)")
    print("="*60)
    
    anomalous = gene_df[
        (gene_df["ProtSyn_Ortholog_Gene"].fillna("").str.strip() != "") &
        (gene_df["Nscore_ProtSyn"].notna()) &
        (gene_df["Nscore_ProtSyn"] == 0)
    ].copy()
    
    if len(anomalous) > 0:
        # Sort by FN_ID for easier viewing
        anomalous = anomalous.sort_values("FN_ID")
        
        anomalous.to_csv(OUT_PATH_ANOMALOUS, sep="\t", index=False)
        
        print(f" WARNING: Saved {len(anomalous)} anomalous gene rows: {OUT_PATH_ANOMALOUS}")
        print(f"  Affecting {anomalous['FN_ID'].nunique()} unique FN pairs")
        
        if "Pair" in anomalous.columns:
            print(f"\n  By species pair:")
            print(anomalous["Pair"].value_counts())
        
        return anomalous
    else:
        print("✓ No anomalous cases found (all ProtSyn N-scores > 0, as expected)")
        return pd.DataFrame()


# =============================================================================
# Plotting Functions
# =============================================================================

def plot_distributions(gene_df):
    """
    Create 2x2 grid plot of N-score distributions (gene-level).
    
    Args:
        gene_df: Gene-level classification dataframe with N-scores
    """
    print("\n" + "="*60)
    print("STEP 7: Plotting Distributions (2x2 Grid)")
    print("="*60)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('False Negatives N-score Distributions (Gene-Level)', fontsize=16, fontweight='bold')
    
    # ---- Plot 1: All FN (Cactus N-scores) ----
    ax1 = axes[0, 0]
    all_cactus = gene_df["Nscore_Cactus"].dropna()
    
    if not all_cactus.empty:
        max_score = int(all_cactus.max())
        bins = np.arange(0, max_score + 2, 1)
        
        ax1.hist(all_cactus, bins=bins, edgecolor='black', color='steelblue', alpha=0.7)
        ax1.set_xlabel('N-score (Cactus)', fontsize=11)
        ax1.set_ylabel('Frequency', fontsize=11)
        ax1.set_title('All False Negatives\n(Cactus pairs)', fontsize=12, fontweight='bold')
        ax1.grid(axis='y', alpha=0.3, linestyle='--')
    
    # ---- Plot 2: No Result Cases ----
    ax2 = axes[0, 1]
    no_result = gene_df[gene_df["FN_Class"] == "no_result"]["Nscore_Cactus"].dropna()
    
    if not no_result.empty:
        max_score = int(no_result.max())
        bins = np.arange(0, max_score + 2, 1)
        
        ax2.hist(no_result, bins=bins, edgecolor='black', color='coral', alpha=0.7)
        ax2.set_xlabel('N-score (Cactus)', fontsize=11)
        ax2.set_ylabel('Frequency', fontsize=11)
        ax2.set_title('No Result Cases\n(gene not in prot-syn)', fontsize=12, fontweight='bold')
        ax2.grid(axis='y', alpha=0.3, linestyle='--')
    else:
        ax2.text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14)
        ax2.set_title('No Result Cases', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Frequency', fontsize=11)
    
    # ---- Plot 3: Contradictory Cases (ProtSyn N-scores) ----
    ax3 = axes[1, 0]
    contradictory = gene_df[gene_df["FN_Class"] == "contradictory"]
    ps_scores = contradictory["Nscore_ProtSyn"].dropna()
    
    if not ps_scores.empty:
        max_score = int(ps_scores.max())
        bins = np.arange(0, max_score + 2, 1)
        
        ax3.hist(ps_scores, bins=bins, edgecolor='black', color='mediumseagreen', alpha=0.7)
        ax3.set_xlabel('N-score (ProtSyn)', fontsize=11)
        ax3.set_ylabel('Frequency', fontsize=11)
        ax3.set_title('Contradictory Cases\n(prot-syn partners)', fontsize=12, fontweight='bold')
        ax3.grid(axis='y', alpha=0.3, linestyle='--')
    else:
        ax3.text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14)
        ax3.set_title('Contradictory Cases', fontsize=12, fontweight='bold')
        ax3.set_ylabel('Frequency', fontsize=11)
    
    # ---- Plot 4: Fraction Distribution ----
    ax4 = axes[1, 1]
    fractions = contradictory["Score_Fraction"].dropna()
    
    if not fractions.empty:
        bins = np.linspace(0, fractions.max(), 30)
        
        ax4.hist(fractions, bins=bins, edgecolor='black', color='mediumpurple', alpha=0.7)
        ax4.axvline(x=1.0, color='red', linestyle='--', linewidth=2.5, 
                    label='Equal (Fraction=1)', zorder=10)
        ax4.set_xlabel('Fraction (Cactus / ProtSyn)', fontsize=11)
        ax4.set_ylabel('Frequency', fontsize=11)
        ax4.set_title('Fraction Distribution\n(contradictory cases)', fontsize=12, fontweight='bold')
        ax4.legend(loc='center right', fontsize=10, framealpha=0.95)
        ax4.grid(axis='y', alpha=0.3, linestyle='--')
    else:
        ax4.text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14)
        ax4.set_title('Fraction Distribution', fontsize=12, fontweight='bold')
        ax4.set_ylabel('Frequency', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(PLOT_DISTRIBUTIONS, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Saved 2x2 distribution plot: {PLOT_DISTRIBUTIONS}")


def print_summary_statistics(gene_df):
    """
    Print summary statistics for all distributions (gene-level).
    
    Args:
        gene_df: Gene-level classification dataframe with N-scores
    """
    print("\n" + "="*60)
    print("SUMMARY STATISTICS (Gene-Level)")
    print("="*60)
    
    # All Cactus N-scores
    all_cactus = gene_df["Nscore_Cactus"].dropna()
    print(f"\n1. All FN (Cactus N-scores):")
    if not all_cactus.empty:
        print(f"   Count:     {len(all_cactus)} gene rows")
        print(f"   Mean:      {all_cactus.mean():.2f}")
        print(f"   Median:    {all_cactus.median():.2f}")
        print(f"   N-score=0: {(all_cactus == 0).sum()} ({(all_cactus == 0).sum()/len(all_cactus)*100:.1f}%)")
    
    # No Result cases
    no_result = gene_df[gene_df["FN_Class"] == "no_result"]["Nscore_Cactus"].dropna()
    print(f"\n2. No Result Cases (Cactus N-scores):")
    if not no_result.empty:
        print(f"   Count:     {len(no_result)} gene rows")
        print(f"   Mean:      {no_result.mean():.2f}")
        print(f"   Median:    {no_result.median():.2f}")
    
    # Contradictory cases - ProtSyn N-scores
    contradictory = gene_df[gene_df["FN_Class"] == "contradictory"]
    ps_scores = contradictory["Nscore_ProtSyn"].dropna()
    
    print(f"\n3. Contradictory Cases (ProtSyn N-scores):")
    if not ps_scores.empty:
        print(f"   Count:     {len(ps_scores)} gene rows")
        print(f"   Mean:      {ps_scores.mean():.2f}")
        print(f"   Median:    {ps_scores.median():.2f}")
        print(f"   N-score=0: {(ps_scores == 0).sum()} ({(ps_scores == 0).sum()/len(ps_scores)*100:.1f}%)")
    
    # Fraction distribution
    fractions = contradictory["Score_Fraction"].dropna()
    
    print(f"\n4. Fraction Distribution:")
    if not fractions.empty:
        print(f"   Count:              {len(fractions)} gene rows")
        print(f"   Mean:               {fractions.mean():.2f}")
        print(f"   Median:             {fractions.median():.2f}")
        print(f"   Cactus better (>1):  {(fractions > 1).sum()} ({(fractions > 1).sum()/len(fractions)*100:.1f}%)")
        print(f"   ProtSyn better (<1): {(fractions < 1).sum()} ({(fractions < 1).sum()/len(fractions)*100:.1f}%)")
        print(f"   Equal (=1):          {(fractions == 1).sum()} ({(fractions == 1).sum()/len(fractions)*100:.1f}%)")


# =============================================================================
# Main Workflow
# =============================================================================

def main():
    """
    Main workflow for gene-level false negatives analysis.
    """
    # Step 1: Load data
    nbh_dict, geomean_dict, map_df, protsyn_df, fn_df = load_data()
    
    # Step 2: Build mappings
    prot2gene, gene2prot, gene_to_species, relationship_dict = build_mappings(map_df, geomean_dict)
    
    # Step 3: Prepare prot-syn data and build index
    protsyn_df = prepare_protsyn_data(protsyn_df, gene_to_species)
    protsyn_index = build_protsyn_index(protsyn_df, gene_to_species)
    
    # Step 4: Add FN_IDs to pair-level table
    fn_df = add_fn_ids(fn_df)
    
    # Step 5: Create gene-level classification
    gene_df = create_gene_level_rows(fn_df, protsyn_index, gene_to_species, gene2prot)
    
    # Step 6: Compute N-scores
    gene_df = compute_nscores(gene_df, nbh_dict, relationship_dict)
    
    # Step 7: Extract special cases
    better_df = extract_cactus_better_cases(gene_df)
    anomalous_df = extract_anomalous_cases(gene_df)
    
    # Step 8: Plot distributions
    plot_distributions(gene_df)
    
    # Step 9: Print summary statistics
    print_summary_statistics(gene_df)
    
    print("\n" + "="*60)
    print("✓ Gene-Level Analysis Complete!")
    print("="*60)
    print(f"\nOutputs:")
    print(f"  1. Pair-level table:     {FN_WITH_IDS_OUT}")
    print(f"  2. Gene-level table:     {FN_CLASSIFIED_OUT}")
    print(f"  3. Cactus better cases:  {OUT_PATH_CACTUS_BETTER}")
    print(f"  4. Anomalous cases:      {OUT_PATH_ANOMALOUS}")
    print(f"  5. Distribution plots:   {PLOT_DISTRIBUTIONS}")


if __name__ == "__main__":
    main()
