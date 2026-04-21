import pandas as pd
import re
import os

# Split false positives into species-specific files.


# ============================================================
# CORE FUNCTIONS
# ============================================================

def species_code(species_name: str) -> str:
    """
    Convert species name to 3-letter code.
    
    Args:
        species_name: Full species name (e.g., "Canis lupus")
    
    Returns:
        3-letter code (e.g., "can") or truncated name for unknown species
    
    Examples:
        >>> species_code("Canis lupus")
        'can'
        >>> species_code("Mus musculus")
        'mus'
        >>> species_code("")
        ''
    """
    if not isinstance(species_name, str) or species_name.strip() == "":
        return ""
    
    key = " ".join(species_name.split()[:2]).lower()  # genus species
    code_map = {
        "canis lupus": "can",
        "macaca fascicularis": "mac",
        "rattus norvegicus": "rat",
        "mus musculus": "mus",
    }
    return code_map.get(key, re.sub(r"[^a-z]+", "_", key)[:8])


def make_pair(qc: str, tc: str) -> str:
    """
    Create pair name from two species codes.
    
    Args:
        qc: Query species code
        tc: Target species code
    
    Returns:
        Sorted pair name (e.g., "can_mus")
    
    Examples:
        >>> make_pair("can", "mus")
        'can_mus'
        >>> make_pair("mus", "can")
        'can_mus'
        >>> make_pair("can", "can")
        'can_can'
    """
    if not qc or not tc:
        return ""
    if qc == tc:
        return f"{qc}_{tc}"
    return "_".join(sorted([qc, tc]))


def build_protein_to_species_map(map_df: pd.DataFrame) -> dict:
    """
    Build mapping from protein ID to species name.
    
    Args:
        map_df: DataFrame with columns 'ProteinID' and 'Species'
    
    Returns:
        Dictionary mapping protein ID -> species name
    
    Raises:
        ValueError: If 'Species' column not found
    """
    if "Species" not in map_df.columns:
        raise ValueError("The map must contain a 'Species' column.")
    
    return dict(zip(
        map_df["ProteinID"].astype(str),
        map_df["Species"].astype(str)
    ))


def build_protein_to_gene_map(map_df: pd.DataFrame, gene_col: str = "GeneID") -> dict:
    """
    Build mapping from protein ID to gene ID.
    
    Args:
        map_df: DataFrame with columns 'ProteinID' and gene_col
        gene_col: Name of gene ID column (default: "GeneID")
    
    Returns:
        Dictionary mapping protein ID -> gene ID, or None if gene_col not found
    
    Warnings:
        Prints warning if proteins map to multiple genes
    """
    if gene_col not in map_df.columns:
        print(f"NOTE: No '{gene_col}' column found; will not map genes.")
        return None
    
    # Check for duplicates
    multi = (map_df.dropna(subset=["ProteinID", gene_col])
                  .groupby("ProteinID")[gene_col].nunique())
    n_multi = int((multi > 1).sum())
    if n_multi:
        print(f"WARNING: {n_multi} ProteinID(s) map to multiple GeneIDs. Using first occurrence.")
    
    # Build mapping
    return (map_df.dropna(subset=["ProteinID"])
                  .drop_duplicates(subset=["ProteinID"], keep="first")
                  .set_index("ProteinID")[gene_col]
                  .to_dict())


def add_species_and_genes(fp_df: pd.DataFrame, 
                         prot2species: dict, 
                         prot2gene: dict = None) -> pd.DataFrame:
    """
    Add species, gene, and pair columns to false positives dataframe.
    
    Args:
        fp_df: DataFrame with Query_Protein and Target_Protein columns
        prot2species: Mapping protein ID -> species name
        prot2gene: Optional mapping protein ID -> gene ID
    
    Returns:
        New DataFrame with added columns:
            - Query_Species, Target_Species
            - Query_Gene, Target_Gene (if prot2gene provided)
            - Q_code, T_code (species codes)
            - Pair (species pair name)
    """
    fp_df = fp_df.copy()
    
    # Add species
    fp_df["Query_Species"] = fp_df["Query_Protein"].map(prot2species)
    fp_df["Target_Species"] = fp_df["Target_Protein"].map(prot2species)
    
    # Add genes if mapping provided
    if prot2gene is not None:
        fp_df["Query_Gene"] = fp_df["Query_Protein"].map(prot2gene)
        fp_df["Target_Gene"] = fp_df["Target_Protein"].map(prot2gene)
    
    # Add species codes
    fp_df["Q_code"] = fp_df["Query_Species"].apply(species_code)
    fp_df["T_code"] = fp_df["Target_Species"].apply(species_code)
    
    # Add pair names
    fp_df["Pair"] = [
        make_pair(qc, tc) 
        for qc, tc in zip(fp_df["Q_code"], fp_df["T_code"])
    ]
    
    return fp_df


def check_missing_species(fp_df: pd.DataFrame) -> dict:
    """
    Check for missing species mappings.
    
    Args:
        fp_df: DataFrame with Query_Species and Target_Species columns
    
    Returns:
        Dictionary with counts of missing query and target species
    """
    n_missing_q = fp_df["Query_Species"].isna().sum()
    n_missing_t = fp_df["Target_Species"].isna().sum()
    
    return {
        "missing_query": int(n_missing_q),
        "missing_target": int(n_missing_t)
    }


def split_by_species_pair(fp_df: pd.DataFrame, 
                          out_dir: str, 
                          prot2gene: dict = None) -> dict:
    """
    Split false positives by species pair and save to files.
    
    Args:
        fp_df: DataFrame with Pair column and protein/species columns
        out_dir: Output directory path
        prot2gene: Optional mapping to include gene columns
    
    Returns:
        Dictionary with:
            - 'saved_files': List of saved file paths
            - 'unmapped_file': Path to unmapped file (if any)
            - 'pair_counts': Dictionary of pair name -> count
    """
    os.makedirs(out_dir, exist_ok=True)
    
    result = {
        'saved_files': [],
        'unmapped_file': None,
        'pair_counts': {}
    }
    
    # Save unmapped pairs
    unmapped = fp_df[fp_df["Pair"] == ""].copy()
    if len(unmapped) > 0:
        unmapped_path = os.path.join(out_dir, "false_positives_nscore_UNMAPPED.tsv")
        unmapped.to_csv(unmapped_path, sep="\t", index=False)
        print(f"Saved unmapped FP rows to: {unmapped_path}")
        result['unmapped_file'] = unmapped_path
    
    # Save species-specific pairs
    for pair_name, sub in fp_df[fp_df["Pair"] != ""].groupby("Pair"):
        out_path = os.path.join(out_dir, f"false_positives_nscore_{pair_name}.tsv")
        
        # Select columns to save
        cols = ["Query_Protein", "Target_Protein", "Query_Species", "Target_Species"]
        if prot2gene is not None:
            cols += ["Query_Gene", "Target_Gene"]
        
        sub[cols].to_csv(out_path, sep="\t", index=False)
        print(f"Saved: {out_path}, n={len(sub)}")
        
        result['saved_files'].append(out_path)
        result['pair_counts'][pair_name] = len(sub)
    
    return result


# ============================================================
# MAIN PIPELINE
# ============================================================

def process_false_positives(fp_path: str, map_path: str, out_dir: str, gene_col: str = "GeneID"):
    """
    Main pipeline to process false positives.
    
    Args:
        fp_path: Path to false positives TSV file
        map_path: Path to protein-gene mapping TSV file
        out_dir: Output directory for split files
        gene_col: Name of gene ID column in mapping file
    
    Returns:
        Tuple of (processed_df, results_dict)
    """
    print("="*60)
    print("PROCESSING FALSE POSITIVES")
    print("="*60)
    
    # 1. Load data
    print(f"\n1. Loading data...")
    fp_df = pd.read_csv(fp_path, sep="\t", dtype=str)
    map_df = pd.read_csv(map_path, sep="\t", dtype=str)
    print(f"   Loaded {len(fp_df)} false positive pairs")
    print(f"   Loaded {len(map_df)} protein mappings")
    
    # 2. Build mappings
    print(f"\n2. Building mappings...")
    prot2species = build_protein_to_species_map(map_df)
    prot2gene = build_protein_to_gene_map(map_df, gene_col)
    print(f"   Protein->Species: {len(prot2species)} entries")
    if prot2gene:
        print(f"   Protein->Gene: {len(prot2gene)} entries")
    
    # 3. Add species and gene information
    print(f"\n3. Adding species and gene information...")
    fp_df = add_species_and_genes(fp_df, prot2species, prot2gene)
    
    # 4. Check for missing mappings
    print(f"\n4. Checking for missing mappings...")
    missing = check_missing_species(fp_df)
    print(f"   Missing Query species: {missing['missing_query']}")
    print(f"   Missing Target species: {missing['missing_target']}")
    
    # 5. Show distribution
    print(f"\n5. Species pair distribution:")
    pair_counts = fp_df["Pair"].value_counts(dropna=False)
    for pair, count in pair_counts.items():
        print(f"   {pair}: {count}")
    
    # 6. Split and save
    print(f"\n6. Splitting by species pair...")
    results = split_by_species_pair(fp_df, out_dir, prot2gene)
    
    print(f"\n{'='*60}")
    print(f"COMPLETE! Saved {len(results['saved_files'])} species-pair files")
    print(f"{'='*60}\n")
    
    return fp_df, results


# ============================================================
# SCRIPT ENTRY POINT
# ============================================================

def main():
    """Script entry point with hardcoded paths."""
    # Input paths
    fp_path = "results/evaluation/mammalia/vs_prot_syn/false_positives_nscore.tsv"
    map_path = "material/sex_experiment/gene_to_protein_map_FIXED.tsv"
    out_dir = "results/evaluation/mammalia/vs_prot_syn/disagreement_investigation/false_positives/"

    
    # Run pipeline
    process_false_positives(fp_path, map_path, out_dir, gene_col="GeneID")


if __name__ == "__main__":
    main()