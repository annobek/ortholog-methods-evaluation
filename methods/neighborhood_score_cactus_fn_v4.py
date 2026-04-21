import pandas as pd
import pickle
import networkx as nx
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
# Uncomment for corresponding prot-syn approach
# prot-syn nscore approach
#PROTSYN_ORTH_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_nscore.tsv"
#FALSE_NEGATIVES_PATH = "results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/all_false_negatives_nscore.tsv"
# prot-syn sum approach
PROTSYN_ORTH_PATH = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_sum.tsv"
FALSE_NEGATIVES_PATH = "results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/all_false_negatives_sum.tsv"

# Output directory
# prot-syn nscore
#OUT_DIR = Path(
#    "results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/"
#    "disagreement_investigation/false_negatives"
#)
# prot-syn sum
OUT_DIR = Path(
    "results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/"
    "disagreement_investigation/false_negatives"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Output files
FN_WITH_IDS_OUT = OUT_DIR / "false_negatives_with_ids_pair_level.tsv"
FN_CLASSIFIED_OUT = OUT_DIR / "false_negatives_classified_gene_level.tsv"
OUT_PATH_CACTUS_BETTER = OUT_DIR / "contradictory_cactus_higher_nscore_than_protsyn_gene_level.tsv"
OUT_PATH_ANOMALOUS = OUT_DIR / "anomalous_protsyn_nscore_zero_gene_level.tsv"
PLOT_DISTRIBUTIONS = OUT_DIR / "fn_distributions_2x2_gene_level.png"
FN_PAIR_LEVEL_OUT = OUT_DIR / "false_negatives_classified_pair_level.tsv"
PLOT_PAIR_DISTRIBUTIONS_PAIR = OUT_DIR / "fn_distributions_2x2_pair_level.png"
EDGE_AFFECTED_GENE_ROWS_OUT = OUT_DIR / "false_negatives_edge_affected_gene_rows.tsv"

# =============================================================================
# Helper Functions
# =============================================================================

def count_shared_genes(gene, gene_related, neighborhood_1, neighborhood_2, relationship_dict):
    """
    Compute the maximum number of unique cross-neighborhood gene pairs using bipartite matching.

    Notes:
        `gene` and `gene_related` are kept for compatibility with older code but are not used
        in the calculation.

    Args:
        gene (str): Focus gene (unused; kept for compatibility).
        gene_related (str): Partner gene (unused; kept for compatibility).
        neighborhood_1 (list[str]): Neighborhood genes around gene in species 1.
        neighborhood_2 (list[str]): Neighborhood genes around gene_related in species 2.
        relationship_dict (dict[str, list[tuple[str, float]]]):
            Mapping gene -> list of (related_gene, score).

    Returns:
        tuple[int, int]:
            (max_unique_pairs, num_edges)
            - max_unique_pairs: size of maximum matching (N-score component)
            - num_edges: number of edges in bipartite graph (debug/diagnostic)
    """
    neighborhood_set_2 = set(neighborhood_2)

    B = nx.Graph()
    B.add_nodes_from(neighborhood_1, bipartite=0)
    B.add_nodes_from(neighborhood_2, bipartite=1)

    for gene1 in neighborhood_1:
        if gene1 in neighborhood_set_2:
            continue

        related_genes = relationship_dict.get(gene1)
        if not related_genes:
            continue

        for gene2, _ in related_genes:
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
    Create a mapping from common/raw species strings to canonical short names.

    Returns:
        dict[str, str]: Mapping of normalized species strings -> canonical species labels.
    """
    return {
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


def normalize_species_robust(species_str):
    """
    Normalize species names to a canonical form.

    Args:
        species_str (str): Species name in any format (underscores/spaces/case variations).

    Returns:
        str: Canonical species name (e.g., "dog", "mouse") or normalized raw string if unknown.
    """
    if pd.isna(species_str):
        return ""

    normalized = species_str.lower().replace("_", " ").strip()
    species_map = create_species_mapper()

    if normalized in species_map:
        return species_map[normalized]

    parts = normalized.split()
    if len(parts) >= 2:
        two_words = " ".join(parts[:2])
        if two_words in species_map:
            return species_map[two_words]

    if len(parts) >= 1:
        first_word = parts[0]
        if first_word in species_map:
            return species_map[first_word]

    return normalized


def build_protsyn_index(protsyn_df):
    """
    Build a symmetric lookup index for prot-syn ortholog assignment.

    The index allows fast lookup:
        (gene, target_species_canonical) -> partner_gene

    Both directions are indexed, so order does not matter.

    Args:
        protsyn_df (pd.DataFrame): Must contain columns Gene1, Gene2, and optionally Species1, Species2.

    Returns:
        dict[tuple[str, str], str]: Symmetric index for prot-syn partner lookup.
    """
    print("\nBuilding symmetric prot-syn index for fast lookup...")

    index = {}
    for _, row in protsyn_df.iterrows():
        g1 = row["Gene1"]
        g2 = row["Gene2"]
        sp1_norm = normalize_species_robust(row.get("Species1", ""))
        sp2_norm = normalize_species_robust(row.get("Species2", ""))

        key1 = (g1, sp2_norm)
        if key1 not in index:
            index[key1] = g2

        key2 = (g2, sp1_norm)
        if key2 not in index:
            index[key2] = g1

    print(f"✓ Built symmetric index with {len(index)} gene-species pairs")
    return index


def find_protsyn_partner(gene, target_species, protsyn_index):
    """
    Lookup the prot-syn partner gene for a given (gene, target_species).

    Args:
        gene (str): Focus gene.
        target_species (str): Canonical target species name.
        protsyn_index (dict[tuple[str, str], str]): Output of build_protsyn_index().

    Returns:
        str | None: Partner gene if found, otherwise None.
    """
    return protsyn_index.get((gene, target_species))


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_data():
    """
    Load dictionaries and TSV inputs.

    Returns:
        tuple:
            nbh_dict (dict): gene -> neighborhood gene list
            geomean_dict (dict): gene -> list of related (gene, score)
            map_df (pd.DataFrame): protein-gene mapping
            protsyn_df (pd.DataFrame): prot-syn one-to-one ortholog pairs
            fn_df (pd.DataFrame): false negative pairs (Cactus reference)
    """
    print("=" * 60)
    print("STEP 1: Loading Data")
    print("=" * 60)

    with open(NEIGHBORHOOD_DICT_PKL, "rb") as f:
        nbh_dict = pickle.load(f)

    with open(GEOMEAN_DICT_PKL, "rb") as f:
        geomean_dict = pickle.load(f)

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
    Build mapping dictionaries needed downstream.

    Args:
        map_df (pd.DataFrame): Must contain GeneID, ProteinID, Species columns.
        geomean_dict (dict): gene -> list of related (gene, score)

    Returns:
        tuple:
            prot2gene (dict[str, str]): ProteinID -> GeneID
            gene2prot (dict[str, str]): GeneID -> ProteinID
            gene_to_species (dict[str, str]): GeneID -> Species
            relationship_dict (dict[str, list[tuple[str, float]]]): gene -> related list
    """
    print("\nBuilding mappings...")

    prot2gene = dict(zip(map_df["ProteinID"], map_df["GeneID"]))
    gene2prot = dict(zip(map_df["GeneID"], map_df["ProteinID"]))
    gene_to_species = dict(zip(map_df["GeneID"], map_df["Species"]))

    relationship_dict = {str(gene): related_list for gene, related_list in geomean_dict.items()}

    print(f"✓ Built mappings for {len(gene_to_species)} genes")
    return prot2gene, gene2prot, gene_to_species, relationship_dict


def prepare_protsyn_data(protsyn_df, gene_to_species):
    """
    Ensure prot-syn table has species columns and correct dtypes.

    Args:
        protsyn_df (pd.DataFrame): Prot-syn ortholog pairs with Gene1, Gene2, optional Species1/Species2.
        gene_to_species (dict[str, str]): GeneID -> Species (raw).

    Returns:
        pd.DataFrame: Cleaned prot-syn dataframe with Gene1/Gene2 as strings.
    """
    if "Species1" not in protsyn_df.columns or "Species2" not in protsyn_df.columns:
        print("Adding species info to prot-syn...")
        protsyn_df["Species1"] = protsyn_df["Gene1"].map(gene_to_species)
        protsyn_df["Species2"] = protsyn_df["Gene2"].map(gene_to_species)

    protsyn_df = protsyn_df.dropna(subset=["Gene1", "Gene2"])
    protsyn_df["Gene1"] = protsyn_df["Gene1"].astype(str)
    protsyn_df["Gene2"] = protsyn_df["Gene2"].astype(str)

    return protsyn_df


def add_fn_ids(fn_df):
    """
    Add a stable FN_ID to each FN pair row and write the pair-level FN-with-IDs table.

    Args:
        fn_df (pd.DataFrame): False negatives dataframe (one row per pair).

    Returns:
        pd.DataFrame: Copy of fn_df with FN_ID column added.

    Output:
        Writes FN_WITH_IDS_OUT (TSV) with FN_ID.
    """
    print("\n" + "=" * 60)
    print("STEP 2: Adding FN_IDs to Pair-Level Table")
    print("=" * 60)

    required_cols = {"Protein1", "Protein2", "Gene1", "Gene2", "Species1", "Species2"}
    if not required_cols.issubset(fn_df.columns):
        raise ValueError(f"FN file missing required columns. Has: {fn_df.columns.tolist()}")

    fn_df = fn_df.copy()
    fn_df["FN_ID"] = range(1, len(fn_df) + 1)

    fn_df["Protein1"] = fn_df["Protein1"].astype(str)
    fn_df["Protein2"] = fn_df["Protein2"].astype(str)
    fn_df["Gene1"] = fn_df["Gene1"].astype(str)
    fn_df["Gene2"] = fn_df["Gene2"].astype(str)

    fn_df.to_csv(FN_WITH_IDS_OUT, sep="\t", index=False)
    print(f"✓ Saved pair-level table with {len(fn_df)} FN pairs: {FN_WITH_IDS_OUT}")

    return fn_df


# =============================================================================
# Gene-Level Classification Functions
# =============================================================================

def create_gene_level_rows(fn_df, protsyn_index, gene2prot):
    """
    Expand FN pairs into gene-level rows (2 rows per FN pair) and classify each gene-side.

    For each FN pair (Gene1, Gene2):
        - Gene1-row: query Gene1 against target species of Gene2
        - Gene2-row: query Gene2 against target species of Gene1

    Args:
        fn_df (pd.DataFrame): FN pairs with FN_ID.
        protsyn_index (dict): Output of build_protsyn_index().
        gene2prot (dict[str, str]): GeneID -> ProteinID.

    Returns:
        pd.DataFrame: Gene-level table with columns describing:
            - Focus gene side (Gene1 or Gene2)
            - Cactus reference ortholog
            - ProtSyn ortholog (if any)
            - FN_Class: no_result / contradictory / agreement
    """
    print("\n" + "=" * 60)
    print("STEP 3: Creating Gene-Level Classification Table")
    print("=" * 60)

    gene_rows = []

    for _, row in fn_df.iterrows():
        fn_id = row["FN_ID"]
        g1, g2 = row["Gene1"], row["Gene2"]
        p1, p2 = row["Protein1"], row["Protein2"]
        sp1_raw, sp2_raw = row["Species1"], row["Species2"]

        sp1 = normalize_species_robust(sp1_raw)
        sp2 = normalize_species_robust(sp2_raw)
        pair_id = row.get("Pair", "")

        # Focus Gene1
        ps_partner_g1 = find_protsyn_partner(g1, sp2, protsyn_index)
        if ps_partner_g1 is None:
            fn_class_g1 = "no_result"
        elif ps_partner_g1 == g2:
            fn_class_g1 = "agreement"
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

        # Focus Gene2
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
    print("\nClassification summary (gene-level rows):")
    print(gene_df["FN_Class"].value_counts(dropna=False))

    if "Pair" in gene_df.columns and gene_df["Pair"].notna().any():
        print("\nBy species pair (gene-level rows):")
        print(gene_df.groupby("Pair")["FN_Class"].value_counts())

    return gene_df


# =============================================================================
# N-score Computation Functions
# =============================================================================

def compute_nscore_for_pair(g1, g2, nbh_dict, relationship_dict):
    """
    Compute the N-score for a gene pair (g1, g2) using their neighborhood lists.

    Args:
        g1 (str): Gene ID 1.
        g2 (str): Gene ID 2.
        nbh_dict (dict[str, list[str]]): Gene -> neighborhood genes.
        relationship_dict (dict[str, list[tuple[str, float]]]): Gene -> related gene list.

    Returns:
        int | None:
            N-score (>=0) if both neighborhoods exist, else None.
    """
    neigh1 = nbh_dict.get(g1, [])
    neigh2 = nbh_dict.get(g2, [])
    if not neigh1 or not neigh2:
        return None
    score, _ = count_shared_genes(g1, g2, neigh1, neigh2, relationship_dict)
    return score


def compute_nscores(gene_df, nbh_dict, relationship_dict):
    """
    Compute N-scores for both Cactus reference pairs and ProtSyn-picked pairs per gene-row,
    and compute gene-level fraction (Cactus/ProtSyn) where defined.

    Args:
        gene_df (pd.DataFrame): Output of create_gene_level_rows().
        nbh_dict (dict[str, list[str]]): Neighborhood dictionary.
        relationship_dict (dict[str, list[tuple[str, float]]]): Relationship dictionary.

    Returns:
        pd.DataFrame: gene_df with columns added:
            - Nscore_Cactus
            - Nscore_ProtSyn
            - Score_Fraction (Cactus/ProtSyn, only if ProtSyn > 0)
    """
    print("\n" + "=" * 60)
    print("STEP 4: Computing N-scores")
    print("=" * 60)

    print("\nComputing N-scores for Cactus pairs...")
    cactus_scores = []
    missing_cactus = 0
    for _, row in gene_df.iterrows():
        score = compute_nscore_for_pair(
            row["Focus_Gene"],
            row["Cactus_Ortholog_Gene"],
            nbh_dict,
            relationship_dict,
        )
        if score is None:
            missing_cactus += 1
        cactus_scores.append(score)
    gene_df["Nscore_Cactus"] = cactus_scores
    print(f"  Could not compute for {missing_cactus} gene rows")

    print("\nComputing N-scores for ProtSyn pairs...")
    ps_scores = []
    missing_ps = 0
    for _, row in gene_df.iterrows():
        ps_partner = row["ProtSyn_Ortholog_Gene"]
        if not ps_partner or ps_partner.strip() == "":
            ps_scores.append(None)
            continue

        score = compute_nscore_for_pair(
            row["Focus_Gene"],
            ps_partner,
            nbh_dict,
            relationship_dict,
        )
        if score is None:
            missing_ps += 1
        ps_scores.append(score)

    gene_df["Nscore_ProtSyn"] = ps_scores
    print(f"  Could not compute for {missing_ps} gene rows (with ProtSyn partners)")

    print("\nComputing score fractions (gene-level)...")
    gene_df["Score_Fraction"] = gene_df.apply(
        lambda r: (r["Nscore_Cactus"] / r["Nscore_ProtSyn"])
        if (pd.notna(r["Nscore_Cactus"]) and pd.notna(r["Nscore_ProtSyn"]) and r["Nscore_ProtSyn"] > 0)
        else np.nan,
        axis=1
    )

    print("\n--- Checking for Unexpected Cases ---")
    ps_zero = gene_df[
        (gene_df["ProtSyn_Ortholog_Gene"].fillna("").str.strip() != "") &
        (gene_df["Nscore_ProtSyn"] == 0)
    ]
    if len(ps_zero) > 0:
        print(f" WARNING: Found {len(ps_zero)} gene rows with ProtSyn N-score = 0")
        print(f"  Affecting {ps_zero['FN_ID'].nunique()} unique FN pairs")
    else:
        print("✓ All ProtSyn pairs have N-score > 0 (as expected)")

    return gene_df


# =============================================================================
# Relative N-score (edge-affected genes)
# =============================================================================

def add_theoretical_max_and_relative_scores(gene_df, nbh_dict, window_max=20):
    """
    Add theoretical maximum neighborhood score per Focus_Gene and relative scores.

    Theory:
        The matching N-score is bounded by the number of available neighbors.
        Near scaffold/chromosome edges, fewer neighbors exist, so the maximum possible
        N-score can be < window_max (default 20).

    Args:
        gene_df (pd.DataFrame): Gene-level dataframe containing Focus_Gene and N-scores.
        nbh_dict (dict[str, list[str]]): Neighborhood dictionary.
        window_max (int): Maximum neighborhood size considered by the algorithm (default 20).

    Returns:
        pd.DataFrame: gene_df with columns:
            - Nscore_TheoreticalMax
            - Nscore_Cactus_Relative
            - Nscore_ProtSyn_Relative

    Output:
        Prints counts of affected genes/rows/pairs (theoretical max < window_max).
    """
    print("\n" + "=" * 60)
    print("STEP 10: Relative N-score (theoretical max < 20 near scaffold edges)")
    print("=" * 60)

    def theo_max(g):
        neigh = nbh_dict.get(str(g), [])
        return min(len(neigh), window_max) if neigh is not None else 0

    gene_df["Nscore_TheoreticalMax"] = gene_df["Focus_Gene"].astype(str).apply(theo_max)

    def safe_rel(actual, tmax):
        if pd.isna(actual) or pd.isna(tmax) or tmax <= 0:
            return np.nan
        return actual / tmax

    gene_df["Nscore_Cactus_Relative"] = gene_df.apply(
        lambda r: safe_rel(r["Nscore_Cactus"], r["Nscore_TheoreticalMax"]), axis=1
    )
    gene_df["Nscore_ProtSyn_Relative"] = gene_df.apply(
        lambda r: safe_rel(r["Nscore_ProtSyn"], r["Nscore_TheoreticalMax"]), axis=1
    )

    affected_mask = gene_df["Nscore_TheoreticalMax"].notna() & (gene_df["Nscore_TheoreticalMax"] < window_max)
    affected_gene_rows = int(affected_mask.sum())
    total_gene_rows = len(gene_df)

    affected_unique_genes = gene_df.loc[affected_mask, "Focus_Gene"].nunique()
    affected_unique_pairs = gene_df.loc[affected_mask, "FN_ID"].nunique()

    print(f"\nTheoretical max window: {window_max}")
    print(
        f"Affected gene rows (theoretical max < {window_max}): {affected_gene_rows}/{total_gene_rows} "
        f"({(affected_gene_rows / total_gene_rows * 100) if total_gene_rows else 0:.2f}%)"
    )
    print(f"Affected unique genes: {affected_unique_genes}")
    print(f"Affected unique FN pairs: {affected_unique_pairs}")

    if affected_gene_rows > 0:
        print("\nReduced theoretical maxima (value counts):")
        print(gene_df.loc[affected_mask, "Nscore_TheoreticalMax"].value_counts().sort_index())

    return gene_df



def save_edge_affected_gene_rows(gene_df, window_max=20):
    """
    Save a filtered TSV containing only gene-level rows where the theoretical maximum
    neighborhood score is reduced (i.e., < window_max). These are the scaffold/contig-edge
    affected cases.

    Args:
        gene_df (pd.DataFrame):
            Gene-level FN table. Must contain:
              - Nscore_TheoreticalMax
        window_max (int):
            Max neighborhood size used by the algorithm (default 20).

    Returns:
        pd.DataFrame:
            Filtered dataframe (edge-affected rows only).

    Output:
        Writes EDGE_AFFECTED_GENE_ROWS_OUT (TSV).
        If no rows are affected, writes an empty TSV with headers.
    """
    if "Nscore_TheoreticalMax" not in gene_df.columns:
        raise ValueError("Expected column 'Nscore_TheoreticalMax' not found. "
                         "Run add_theoretical_max_and_relative_scores() first.")

    affected = gene_df[
        gene_df["Nscore_TheoreticalMax"].notna() &
        (gene_df["Nscore_TheoreticalMax"] < window_max)
    ].copy()

    affected.to_csv(EDGE_AFFECTED_GENE_ROWS_OUT, sep="\t", index=False)

    print(f"\n✓ Saved edge-affected gene rows (theoretical max < {window_max}): {EDGE_AFFECTED_GENE_ROWS_OUT}")
    print(f"  Rows saved: {len(affected)}")
    print(f"  Unique genes: {affected['Focus_Gene'].nunique() if len(affected) else 0}")
    print(f"  Unique FN pairs: {affected['FN_ID'].nunique() if len(affected) else 0}")

    return affected


# =============================================================================
# Analysis (special-case extractions)
# =============================================================================

def extract_cactus_better_cases(gene_df):
    """
    Extract gene-rows where Cactus reference has higher N-score than ProtSyn-picked ortholog.

    Args:
        gene_df (pd.DataFrame): Gene-level dataframe with Nscore_Cactus, Nscore_ProtSyn, Score_Fraction.

    Returns:
        pd.DataFrame: Filtered rows ("better") sorted by Score_Fraction desc.

    Output:
        Writes OUT_PATH_CACTUS_BETTER if any cases exist.
    """
    print("\n" + "=" * 60)
    print("STEP 5: Extracting Cases Where Cactus Picked Better N-score")
    print("=" * 60)

    better = gene_df[
        (gene_df["FN_Class"] == "contradictory") &
        (gene_df["Nscore_Cactus"].notna()) &
        (gene_df["Nscore_ProtSyn"].notna()) &
        (gene_df["Nscore_Cactus"] > 0) &
        (gene_df["Nscore_ProtSyn"] > 0) &
        (gene_df["Score_Fraction"] > 1)
    ].copy()

    if len(better) == 0:
        print("✓ No cases found where Cactus picked better N-score")
        print("  All contradictory cases show ProtSyn with equal or better N-scores.")
        return pd.DataFrame()

    better = better.sort_values("Score_Fraction", ascending=False)
    better.to_csv(OUT_PATH_CACTUS_BETTER, sep="\t", index=False)

    print(f"✓ Saved {len(better)} gene rows where Cactus picked better: {OUT_PATH_CACTUS_BETTER}")
    print(f"  Affecting {better['FN_ID'].nunique()} unique FN pairs")

    if "Pair" in better.columns:
        print("\n  By species pair:")
        print(better["Pair"].value_counts())

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


def extract_anomalous_cases(gene_df):
    """
    Extract gene-rows where ProtSyn-picked ortholog has N-score == 0.

    Args:
        gene_df (pd.DataFrame): Gene-level dataframe with ProtSyn_Ortholog_Gene and Nscore_ProtSyn.

    Returns:
        pd.DataFrame: Filtered anomalous rows.

    Output:
        Writes OUT_PATH_ANOMALOUS if any cases exist.
    """
    print("\n" + "=" * 60)
    print("STEP 6: Extracting Anomalous Cases (ProtSyn N-score = 0)")
    print("=" * 60)

    anomalous = gene_df[
        (gene_df["ProtSyn_Ortholog_Gene"].fillna("").str.strip() != "") &
        (gene_df["Nscore_ProtSyn"].notna()) &
        (gene_df["Nscore_ProtSyn"] == 0)
    ].copy()

    if len(anomalous) == 0:
        print("✓ No anomalous cases found (all ProtSyn N-scores > 0, as expected)")
        return pd.DataFrame()

    anomalous = anomalous.sort_values("FN_ID")
    anomalous.to_csv(OUT_PATH_ANOMALOUS, sep="\t", index=False)

    print(f" WARNING: Saved {len(anomalous)} anomalous gene rows: {OUT_PATH_ANOMALOUS}")
    print(f"  Affecting {anomalous['FN_ID'].nunique()} unique FN pairs")

    if "Pair" in anomalous.columns:
        print("\n  By species pair:")
        print(anomalous["Pair"].value_counts())

    return anomalous


# =============================================================================
# Plotting Functions
# =============================================================================

def plot_distributions(gene_df):
    """
    Save a 2x2 histogram grid of gene-level distributions.

    Args:
        gene_df (pd.DataFrame): Gene-level dataframe with Nscore_Cactus, Nscore_ProtSyn, Score_Fraction.

    Output:
        Writes PLOT_DISTRIBUTIONS (PNG).
    """
    print("\n" + "=" * 60)
    print("STEP 7: Plotting Distributions (2x2 Grid)")
    print("=" * 60)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("False Negatives N-score Distributions (Gene-Level)", fontsize=16, fontweight="bold")

    ax1 = axes[0, 0]
    all_cactus = gene_df["Nscore_Cactus"].dropna()
    if not all_cactus.empty:
        max_score = int(all_cactus.max())
        bins = np.arange(0, max_score + 2, 1)
        ax1.hist(all_cactus, bins=bins, edgecolor="black", color="steelblue", alpha=0.7)
        ax1.set_xlabel("N-score (Cactus)", fontsize=11)
        ax1.set_ylabel("Frequency", fontsize=11)
        ax1.set_title("All False Negatives\n(Cactus pairs)", fontsize=12, fontweight="bold")
        ax1.grid(axis="y", alpha=0.3, linestyle="--")

    ax2 = axes[0, 1]
    no_result = gene_df[gene_df["FN_Class"] == "no_result"]["Nscore_Cactus"].dropna()
    if not no_result.empty:
        max_score = int(no_result.max())
        bins = np.arange(0, max_score + 2, 1)
        ax2.hist(no_result, bins=bins, edgecolor="black", color="coral", alpha=0.7)
        ax2.set_xlabel("N-score (Cactus)", fontsize=11)
        ax2.set_ylabel("Frequency", fontsize=11)
        ax2.set_title("No Result Cases\n(gene not in prot-syn)", fontsize=12, fontweight="bold")
        ax2.grid(axis="y", alpha=0.3, linestyle="--")
    else:
        ax2.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=14)
        ax2.set_title("No Result Cases", fontsize=12, fontweight="bold")
        ax2.set_ylabel("Frequency", fontsize=11)

    ax3 = axes[1, 0]
    contradictory = gene_df[gene_df["FN_Class"] == "contradictory"]
    ps_scores = contradictory["Nscore_ProtSyn"].dropna()
    if not ps_scores.empty:
        max_score = int(ps_scores.max())
        bins = np.arange(0, max_score + 2, 1)
        ax3.hist(ps_scores, bins=bins, edgecolor="black", color="mediumseagreen", alpha=0.7)
        ax3.set_xlabel("N-score (ProtSyn)", fontsize=11)
        ax3.set_ylabel("Frequency", fontsize=11)
        ax3.set_title("Contradictory Cases\n(prot-syn partners)", fontsize=12, fontweight="bold")
        ax3.grid(axis="y", alpha=0.3, linestyle="--")
    else:
        ax3.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=14)
        ax3.set_title("Contradictory Cases", fontsize=12, fontweight="bold")
        ax3.set_ylabel("Frequency", fontsize=11)

    ax4 = axes[1, 1]
    fractions = contradictory["Score_Fraction"].dropna()
    if not fractions.empty:
        bins = np.linspace(0, fractions.max(), 30)
        ax4.hist(fractions, bins=bins, edgecolor="black", color="mediumpurple", alpha=0.7)
        ax4.axvline(x=1.0, color="red", linestyle="--", linewidth=2.5, label="Equal (Fraction=1)", zorder=10)
        ax4.set_xlabel("Fraction (Cactus / ProtSyn)", fontsize=11)
        ax4.set_ylabel("Frequency", fontsize=11)
        ax4.set_title("Fraction Distribution\n(contradictory cases)", fontsize=12, fontweight="bold")
        ax4.legend(loc="center right", fontsize=10, framealpha=0.95)
        ax4.grid(axis="y", alpha=0.3, linestyle="--")
    else:
        ax4.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=14)
        ax4.set_title("Fraction Distribution", fontsize=12, fontweight="bold")
        ax4.set_ylabel("Frequency", fontsize=11)

    plt.tight_layout()
    plt.savefig(PLOT_DISTRIBUTIONS, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"\n✓ Saved 2x2 distribution plot: {PLOT_DISTRIBUTIONS}")


def build_pair_level_table(fn_df, gene_df):
    """
    Aggregate gene-level results into a pair-level FN table (one row per FN_ID).

    Pair_Class:
        - no_result: both gene-sides are no_result
        - contradictory: otherwise

    Aggregations:
        - Nscore_Cactus_pair: max over the two gene-rows (should be identical)
        - ProtSyn_Nscore_mean: mean of available ProtSyn scores (ignores NA)
        - ProtSyn_Nscore_max: max of available ProtSyn scores
        - Fraction_mean/max: Cactus_pair / ProtSyn_(mean/max), only if denom > 0

    Args:
        fn_df (pd.DataFrame): Original FN pair table with FN_ID and metadata.
        gene_df (pd.DataFrame): Gene-level table with Nscore_Cactus/Nscore_ProtSyn.

    Returns:
        pd.DataFrame: Pair-level table.

    Output:
        Writes FN_PAIR_LEVEL_OUT (TSV).
    """
    print("\n" + "=" * 60)
    print("STEP 8: Building Pair-Level FN Table")
    print("=" * 60)

    cls = gene_df.groupby("FN_ID")["FN_Class"].apply(list)

    def pair_class(classes):
        classes = [c for c in classes if pd.notna(c)]
        if len(classes) == 0:
            return "unknown"
        if all(c == "no_result" for c in classes):
            return "no_result"
        return "contradictory"

    pair_class_series = cls.apply(pair_class).rename("Pair_Class")
    cactus_pair = gene_df.groupby("FN_ID")["Nscore_Cactus"].max().rename("Nscore_Cactus_pair")

    cactus_nunique = gene_df.groupby("FN_ID")["Nscore_Cactus"].nunique(dropna=True)
    bad = cactus_nunique[cactus_nunique > 1]
    if len(bad) > 0:
        print(f"WARNING: {len(bad)} FN_IDs have >1 distinct Nscore_Cactus across gene-rows (unexpected).")

    ps_mean = gene_df.groupby("FN_ID")["Nscore_ProtSyn"].mean(numeric_only=True).rename("ProtSyn_Nscore_mean")
    ps_max = gene_df.groupby("FN_ID")["Nscore_ProtSyn"].max().rename("ProtSyn_Nscore_max")

    pair_df = pd.concat([pair_class_series, cactus_pair, ps_mean, ps_max], axis=1).reset_index()

    def safe_div(num, den):
        if pd.isna(num) or pd.isna(den) or den <= 0:
            return np.nan
        return num / den

    pair_df["Fraction_mean"] = pair_df.apply(
        lambda r: safe_div(r["Nscore_Cactus_pair"], r["ProtSyn_Nscore_mean"]), axis=1
    )
    pair_df["Fraction_max"] = pair_df.apply(
        lambda r: safe_div(r["Nscore_Cactus_pair"], r["ProtSyn_Nscore_max"]), axis=1
    )

    keep_cols = ["FN_ID", "Pair", "Protein1", "Protein2", "Gene1", "Gene2", "Species1", "Species2"]
    meta = fn_df.copy()
    keep_cols = [c for c in keep_cols if c in meta.columns]
    meta = meta[keep_cols]

    pair_df = meta.merge(pair_df, on="FN_ID", how="left")

    side_cls = gene_df.pivot_table(index="FN_ID", columns="FN_Gene_Position", values="FN_Class", aggfunc="first")
    for side in ["Gene1", "Gene2"]:
        if side in side_cls.columns:
            pair_df[f"NoResult_{side}"] = (side_cls[side] == "no_result").values
            pair_df[f"Contradictory_{side}"] = (side_cls[side] == "contradictory").values

    pair_df.to_csv(FN_PAIR_LEVEL_OUT, sep="\t", index=False)
    print(f"✓ Saved pair-level table: {FN_PAIR_LEVEL_OUT}")
    print("\nPair-level class summary:")
    print(pair_df["Pair_Class"].value_counts(dropna=False))

    return pair_df


def plot_pair_distributions(pair_df):
    """
    Save a 2x2 histogram grid of pair-level distributions.

    Args:
        pair_df (pd.DataFrame): Pair-level table with Nscore_Cactus_pair, ProtSyn_Nscore_mean/max, Fraction_mean/max.

    Output:
        Writes PLOT_PAIR_DISTRIBUTIONS_PAIR (PNG).
    """
    print("\n" + "=" * 60)
    print("STEP 9: Plotting Distributions (Pair-Level 2x2 Grid)")
    print("=" * 60)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("False Negatives N-score Distributions (Pair-Level)", fontsize=16, fontweight="bold")

    ax1 = axes[0, 0]
    all_cactus = pair_df["Nscore_Cactus_pair"].dropna()
    if not all_cactus.empty:
        max_score = int(all_cactus.max())
        bins = np.arange(0, max_score + 2, 1)
        ax1.hist(all_cactus, bins=bins, edgecolor="black", alpha=0.7)
        ax1.set_xlabel("N-score (Cactus, pair)", fontsize=11)
        ax1.set_ylabel("Frequency", fontsize=11)
        ax1.set_title("All False Negative Pairs\n(Cactus reference)", fontsize=12, fontweight="bold")
        ax1.grid(axis="y", alpha=0.3, linestyle="--")

    ax2 = axes[0, 1]
    nores = pair_df[pair_df["Pair_Class"] == "no_result"]["Nscore_Cactus_pair"].dropna()
    if not nores.empty:
        max_score = int(nores.max())
        bins = np.arange(0, max_score + 2, 1)
        ax2.hist(nores, bins=bins, edgecolor="black", alpha=0.7)
        ax2.set_xlabel("N-score (Cactus, pair)", fontsize=11)
        ax2.set_ylabel("Frequency", fontsize=11)
        ax2.set_title("No-result Pairs\n(both genes missing in prot-syn)", fontsize=12, fontweight="bold")
        ax2.grid(axis="y", alpha=0.3, linestyle="--")
    else:
        ax2.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=14)
        ax2.set_title("No-result Pairs", fontsize=12, fontweight="bold")

    ax3 = axes[1, 0]
    contr = pair_df[pair_df["Pair_Class"] == "contradictory"]
    ps_max = contr["ProtSyn_Nscore_max"].dropna()
    ps_mean = contr["ProtSyn_Nscore_mean"].dropna()

    if not ps_max.empty or not ps_mean.empty:
        max_val = 0
        if not ps_max.empty:
            max_val = max(max_val, int(ps_max.max()))
        if not ps_mean.empty:
            max_val = max(max_val, int(ps_mean.max()))
        bins = np.arange(0, max_val + 2, 1)

        if not ps_max.empty:
            ax3.hist(ps_max, bins=bins, edgecolor="black", alpha=0.5, label="ProtSyn max")
        if not ps_mean.empty:
            ax3.hist(ps_mean, bins=bins, edgecolor="black", alpha=0.5, label="ProtSyn mean")

        ax3.set_xlabel("N-score (ProtSyn, aggregated per pair)", fontsize=11)
        ax3.set_ylabel("Frequency", fontsize=11)
        ax3.set_title("Contradictory Pairs\nProtSyn mean vs max", fontsize=12, fontweight="bold")
        ax3.legend()
        ax3.grid(axis="y", alpha=0.3, linestyle="--")
    else:
        ax3.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=14)
        ax3.set_title("Contradictory Pairs", fontsize=12, fontweight="bold")

    ax4 = axes[1, 1]
    fr_max = contr["Fraction_max"].dropna()
    fr_mean = contr["Fraction_mean"].dropna()

    if not fr_max.empty or not fr_mean.empty:
        max_f = 0
        if not fr_max.empty:
            max_f = max(max_f, float(fr_max.max()))
        if not fr_mean.empty:
            max_f = max(max_f, float(fr_mean.max()))

        bins = np.linspace(0, max_f, 30) if max_f > 0 else np.linspace(0, 1, 30)

        if not fr_max.empty:
            ax4.hist(fr_max, bins=bins, edgecolor="black", alpha=0.5, label="Fraction max")
        if not fr_mean.empty:
            ax4.hist(fr_mean, bins=bins, edgecolor="black", alpha=0.5, label="Fraction mean")

        ax4.axvline(x=1.0, color="red", linestyle="--", linewidth=2.5, label="Equal (=1)", zorder=10)
        ax4.set_xlabel("Fraction (Cactus / ProtSyn)", fontsize=11)
        ax4.set_ylabel("Frequency", fontsize=11)
        ax4.set_title("Contradictory Pairs\nFraction mean vs max", fontsize=12, fontweight="bold")
        ax4.legend(loc="center right", fontsize=10, framealpha=0.95)
        ax4.grid(axis="y", alpha=0.3, linestyle="--")
    else:
        ax4.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=14)
        ax4.set_title("Fraction Distribution", fontsize=12, fontweight="bold")

    plt.tight_layout()
    plt.savefig(PLOT_PAIR_DISTRIBUTIONS_PAIR, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"✓ Saved pair-level 2x2 distribution plot: {PLOT_PAIR_DISTRIBUTIONS_PAIR}")


def print_summary_statistics(gene_df):
    """
    Print gene-level summary statistics used for quick sanity-checking.

    Args:
        gene_df (pd.DataFrame): Gene-level dataframe with Nscore_Cactus, Nscore_ProtSyn, Score_Fraction.

    Returns:
        None (prints to stdout).
    """
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS (Gene-Level)")
    print("=" * 60)

    all_cactus = gene_df["Nscore_Cactus"].dropna()
    print("\n1. All FN (Cactus N-scores):")
    if not all_cactus.empty:
        print(f"   Count:     {len(all_cactus)} gene rows")
        print(f"   Mean:      {all_cactus.mean():.2f}")
        print(f"   Median:    {all_cactus.median():.2f}")
        print(f"   N-score=0: {(all_cactus == 0).sum()} ({(all_cactus == 0).sum()/len(all_cactus)*100:.1f}%)")

    no_result = gene_df[gene_df["FN_Class"] == "no_result"]["Nscore_Cactus"].dropna()
    print("\n2. No Result Cases (Cactus N-scores):")
    if not no_result.empty:
        print(f"   Count:     {len(no_result)} gene rows")
        print(f"   Mean:      {no_result.mean():.2f}")
        print(f"   Median:    {no_result.median():.2f}")

    contradictory = gene_df[gene_df["FN_Class"] == "contradictory"]
    ps_scores = contradictory["Nscore_ProtSyn"].dropna()
    print("\n3. Contradictory Cases (ProtSyn N-scores):")
    if not ps_scores.empty:
        print(f"   Count:     {len(ps_scores)} gene rows")
        print(f"   Mean:      {ps_scores.mean():.2f}")
        print(f"   Median:    {ps_scores.median():.2f}")
        print(f"   N-score=0: {(ps_scores == 0).sum()} ({(ps_scores == 0).sum()/len(ps_scores)*100:.1f}%)")

    fractions = contradictory["Score_Fraction"].dropna()
    print("\n4. Fraction Distribution:")
    if not fractions.empty:
        print(f"   Count:              {len(fractions)} gene rows")
        print(f"   Mean:               {fractions.mean():.2f}")
        print(f"   Median:             {fractions.median():.2f}")
        print(f"   Cactus better (>1):  {(fractions > 1).sum()} ({(fractions > 1).sum()/len(fractions)*100:.1f}%)")
        print(f"   ProtSyn better (<1): {(fractions < 1).sum()} ({(fractions < 1).sum()/len(fractions)*100:.1f}%)")
        print(f"   Equal (=1):          {(fractions == 1).sum()} ({(fractions == 1).sum()/len(fractions)*100:.1f}%)")


def print_pair_level_summary(pair_df):
    """
    Print key pair-level counts.

    Reported:
        - count of no_result vs contradictory FN pairs
        - count of contradictory pairs where Cactus is better (Fraction > 1), for mean/max aggregation
        - count of contradictory pairs where ProtSyn aggregated score is zero

    Args:
        pair_df (pd.DataFrame): Pair-level dataframe.

    Returns:
        None (prints to stdout).
    """
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS (Pair-Level)")
    print("=" * 60)

    print("\n1. Pair-level class counts:")
    cls_counts = pair_df["Pair_Class"].value_counts(dropna=False)
    for k, v in cls_counts.items():
        print(f"   {k:>13}: {v}")

    contr = pair_df[pair_df["Pair_Class"] == "contradictory"].copy()

    valid_mean = (
        contr["Nscore_Cactus_pair"].notna() &
        contr["ProtSyn_Nscore_mean"].notna() &
        (contr["ProtSyn_Nscore_mean"] > 0) &
        contr["Fraction_mean"].notna()
    )
    valid_max = (
        contr["Nscore_Cactus_pair"].notna() &
        contr["ProtSyn_Nscore_max"].notna() &
        (contr["ProtSyn_Nscore_max"] > 0) &
        contr["Fraction_max"].notna()
    )

    cactus_better_mean = (contr.loc[valid_mean, "Fraction_mean"] > 1).sum()
    cactus_better_max = (contr.loc[valid_max, "Fraction_max"] > 1).sum()

    print("\n2. 'Cactus picked better' counts (pair-level, contradictory only):")
    print(f"   Valid pairs (mean fraction defined): {valid_mean.sum()}")
    print(f"   Cactus better (Fraction_mean > 1):  {cactus_better_mean}")
    print(f"   Valid pairs (max fraction defined):  {valid_max.sum()}")
    print(f"   Cactus better (Fraction_max > 1):   {cactus_better_max}")

    ps_max_zero = (contr["ProtSyn_Nscore_max"].notna() & (contr["ProtSyn_Nscore_max"] == 0)).sum()
    ps_mean_zero = (contr["ProtSyn_Nscore_mean"].notna() & (contr["ProtSyn_Nscore_mean"] == 0)).sum()

    print("\n3. ProtSyn aggregated N-score == 0 (pair-level, contradictory only):")
    print(f"   ProtSyn_Nscore_max == 0:  {ps_max_zero}")
    print(f"   ProtSyn_Nscore_mean == 0: {ps_mean_zero}")


# =============================================================================
# Main Workflow
# =============================================================================

def main():
    """
    End-to-end FN analysis pipeline (N-score focused).

    Outputs written:
        - FN_WITH_IDS_OUT (pair-level input with FN_ID)
        - FN_CLASSIFIED_OUT (gene-level classified + N-scores + relative columns)
        - OUT_PATH_CACTUS_BETTER (gene-level subset)
        - OUT_PATH_ANOMALOUS (gene-level subset)
        - PLOT_DISTRIBUTIONS (gene-level 2x2 plot)
        - FN_PAIR_LEVEL_OUT (pair-level classified/aggregated)
        - PLOT_PAIR_DISTRIBUTIONS_PAIR (pair-level 2x2 plot)
    """
    nbh_dict, geomean_dict, map_df, protsyn_df, fn_df = load_data()
    prot2gene, gene2prot, gene_to_species, relationship_dict = build_mappings(map_df, geomean_dict)

    protsyn_df = prepare_protsyn_data(protsyn_df, gene_to_species)
    protsyn_index = build_protsyn_index(protsyn_df)

    fn_df = add_fn_ids(fn_df)
    gene_df = create_gene_level_rows(fn_df, protsyn_index, gene2prot)

    gene_df = compute_nscores(gene_df, nbh_dict, relationship_dict)

    # NEW: theoretical max + relative scores; prints affected counts
    gene_df = add_theoretical_max_and_relative_scores(gene_df, nbh_dict, window_max=20)
    # NEW: Save the edge-affected subset for reporting / appendix
    edge_df = save_edge_affected_gene_rows(gene_df, window_max=20)

    # Save gene-level classification AFTER relative columns are added
    gene_df.to_csv(FN_CLASSIFIED_OUT, sep="\t", index=False)
    print(f"\n✓ Saved gene-level classification (with relative score columns): {FN_CLASSIFIED_OUT}")

    better_df = extract_cactus_better_cases(gene_df)
    anomalous_df = extract_anomalous_cases(gene_df)

    plot_distributions(gene_df)
    print_summary_statistics(gene_df)

    pair_df = build_pair_level_table(fn_df, gene_df)
    plot_pair_distributions(pair_df)
    print_pair_level_summary(pair_df)

    print("\n" + "=" * 60)
    print("✓ FN N-score Analysis Complete!")
    print("=" * 60)
    print("\nOutputs:")
    print(f"  1. Pair-level with IDs:  {FN_WITH_IDS_OUT}")
    print(f"  2. Gene-level table:     {FN_CLASSIFIED_OUT}")
    print(f"  3. Cactus better cases:  {OUT_PATH_CACTUS_BETTER}")
    print(f"  4. Anomalous cases:      {OUT_PATH_ANOMALOUS}")
    print(f"  5. Gene-level plots:     {PLOT_DISTRIBUTIONS}")
    print(f"  6. Pair-level table:     {FN_PAIR_LEVEL_OUT}")
    print(f"  7. Pair-level plots:     {PLOT_PAIR_DISTRIBUTIONS_PAIR}")
    print(f"  8. Edge-affected rows:   {EDGE_AFFECTED_GENE_ROWS_OUT}")



if __name__ == "__main__":
    main()
