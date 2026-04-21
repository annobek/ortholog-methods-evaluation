import pandas as pd
from pathlib import Path
import pickle
import re

# Evaluate cactus vs BUSCO+Compleasm (BUSComp reference); 
# Create CactBUSComp-set - union set of cactus and BUSComp 

# SPECIES PAIRS CONFIGURATION - Define expected query -> target directions
# Change according to species and their name codes

EXPECTED_SPECIES_PAIRS = [
    ("can", "mac"),  # Canis lupus familiaris -> Macaca fascicularis
    ("can", "rat"),  # Canis lupus familiaris -> Rattus norvegicus
    ("can", "mus"),  # Canis lupus familiaris -> Mus musculus
    ("mac", "rat"),  # Macaca fascicularis -> Rattus norvegicus
    ("mac", "mus"),  # Macaca fascicularis -> Mus musculus
    ("rat", "mus"),  # Rattus norvegicus -> Mus musculus
]

SPECIES_NAMES_CODE = {
    "canis lupus": "can",
    "macaca fascicularis": "mac",
    "mus musculus": "mus",
    "rattus norvegicus": "rat",
}

# Input:

# Reference set: BUSComp (BUSCO+Compleasm)
REF_PATH = "/storage/EasyVectorOmics/busco/results/mammalia/evaluation/BUSComp_set.tsv"

# Test set: Cactus files
TEST_PATH = [
    "results/reciprocal_pairs/can_mac_one2one_fix.tsv",
    "results/reciprocal_pairs/can_rat_one2one_fix.tsv",
    "results/reciprocal_pairs/can_mus_one2one_fix.tsv",
    "results/reciprocal_pairs/mac_rat_one2one_fix.tsv",
    "results/reciprocal_pairs/mac_mus_one2one_fix.tsv",
    "results/reciprocal_pairs/rat_mus_one2one_fix.tsv",
]

# Map protein-species (if BUSComp doesn't have species columns, otherwise not needed)
MAP_PATH = "material/sex_experiment/gene_to_protein_map_FIXED.tsv"

# Output:

OUT_METRICS = "results/evaluation/mammalia/vs_buscomp_set/cactus_BUSComp_metrics.txt"
OUT_METRICS_TSV = "results/evaluation/mammalia/vs_buscomp_set/cactus_BUSComp_metrics.tsv"
OUT_MISSED_BY_TEST = "results/evaluation/mammalia/vs_buscomp_set/all_false_negatives_cactus.tsv"
OUT_MISSED_BY_REF = "results/evaluation/mammalia/vs_buscomp_set/all_false_positives_cactus.tsv"
OUT_SPLIT_DIR_FN = "results/evaluation/mammalia/vs_buscomp_set"
OUT_SPLIT_DIR_FP = "results/evaluation/mammalia/vs_buscomp_set"
OUT_UNMAPPED_FN = "results/evaluation/mammalia/vs_buscomp_set/false_negatives_UNMAPPED.tsv"
OUT_UNMAPPED_FP = "results/evaluation/mammalia/vs_buscomp_set/false_positives_UNMAPPED.tsv"
OUT_CACTBUSCOMP_SET_TSV = "results/evaluation/mammalia/vs_buscomp_set/CactBUSComp_set.tsv"
OUT_CACTBUSCOMP_SET_PKL = "results/evaluation/mammalia/vs_buscomp_set/CactBUSComp_set.pkl"
# ----------------------------


def load_pairs(files, query_col="Query_Protein", target_col="Target_Protein", sep="\t"):
    """
    Load multiple Cactus CSVs and extract only the columns needed for the protein pairs.
    Automatically ignores other columns present in the file.
    """
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep=sep, dtype=str)
        if query_col not in df.columns or target_col not in df.columns:
            raise ValueError(f"Missing required columns in {f}: expected '{query_col}' and '{target_col}'")
        dfs.append(df[[query_col, target_col]])
    return pd.concat(dfs, ignore_index=True)


def load_buscomp_pairs(file_path, sep="\t"):
    """
    Load BUSComp file which may have different column names.
    Handles both old format (Protein1/Protein2) and new format (Query_Protein/Target_Protein).
    """
    df = pd.read_csv(file_path, sep=sep, dtype=str)
    
    # Check for Query_Protein/Target_Protein columns (new format)
    if "Query_Protein" in df.columns and "Target_Protein" in df.columns:
        return df
    # Check for Protein1/Protein2 columns (old format)
    elif "Protein1" in df.columns and "Protein2" in df.columns:
        df = df.rename(columns={"Protein1": "Query_Protein", "Protein2": "Target_Protein"})
        return df
    else:
        raise ValueError(f"BUSComp file must contain either (Query_Protein, Target_Protein) or (Protein1, Protein2) columns. Found: {df.columns.tolist()}")


# ----------------------------
# Robust pair handling
# ----------------------------

def build_pair_index(df, qcol, tcol):
    """
    Build:
      - keys: set of undirected pair keys (tuple(sorted((q, t))))
      - directed: dict mapping undirected key -> (q, t) as it appears in df (first occurrence)

    This lets us:
      - compute metrics using keys (undirected)
      - write FN/FP using directed orientation from the appropriate source
    """
    keys = set()
    directed = {}
    flipped = 0

    for q, t in zip(df[qcol], df[tcol]):
        if pd.isna(q) or pd.isna(t) or q == t:
            continue

        key = tuple(sorted((q, t)))
        keys.add(key)

        # If the same undirected key appears with opposite orientation, warn.
        if key in directed and directed[key] != (q, t):
            flipped += 1

        # Keep first occurrence so output is stable.
        directed.setdefault(key, (q, t))

    if flipped:
        print(f"WARNING: {flipped} pair(s) appear with BOTH orientations in columns {qcol}/{tcol}. Keeping first occurrence.")
    return keys, directed


def compare_sets(ref_set, test_set):
    """Return metrics dictionary comparing reference vs test pairs."""
    tp = len(ref_set & test_set)
    fn = len(ref_set - test_set)
    fp = len(test_set - ref_set)
    recall = tp / (tp + fn) if tp + fn else 0
    precision = tp / (tp + fp) if tp + fp else 0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0
    # Handle empty sets
    union_size = len(ref_set | test_set)
    jaccard = len(ref_set & test_set) / union_size if union_size else 0
    return dict(tp=tp, fn=fn, fp=fp, recall=recall, precision=precision, f1=f1, jaccard=jaccard)


def save_metrics_as_tsv(metrics_dict, ref_count, test_count, output_path):
    """
    Save evaluation metrics as a TSV file with one row of values.
    
    Args:
        metrics_dict: Dictionary with tp, fn, fp, recall, precision, f1, jaccard
        ref_count: Number of pairs in reference set
        test_count: Number of pairs in test set
        output_path: Path to save the TSV file
    """
    # Create DataFrame with one row
    metrics_df = pd.DataFrame({
        "Num_Pairs_Reference": [ref_count],
        "Num_Pairs_Test": [test_count],
        "Num_TP": [metrics_dict['tp']],
        "Num_FN": [metrics_dict['fn']],
        "Num_FP": [metrics_dict['fp']],
        "Precision": [metrics_dict['precision']],
        "Recall": [metrics_dict['recall']],
        "F1_score": [metrics_dict['f1']],
        "Jaccard_Similarity": [metrics_dict['jaccard']],
        "%_of_ref_recovered": [metrics_dict['recall'] * 100]  # Convert to percentage
    })
    
    metrics_df.to_csv(output_path, sep="\t", index=False)
    print(f"Saved metrics TSV to: {output_path}")


# ----------------------------
# Species / gene mapping helpers
# ----------------------------

def species_code(species_name: str) -> str:
    """Map species name (robust by first two words) to a short code."""
    if not isinstance(species_name, str) or species_name.strip() == "":
        return ""
    key = " ".join(species_name.split()[:2]).lower()
    code_map = SPECIES_NAMES_CODE
    return code_map.get(key, re.sub(r"[^a-z]+", "_", key)[:8])


def build_expected_directions(species_pairs):
    """
    Convert a list of (query_code, target_code) tuples into a lookup dict.
    
    Args:
        species_pairs: List of tuples like [("can", "mac"), ("can", "rat"), ...]
    
    Returns:
        dict: Maps frozenset({query, target}) -> (query, target)
              This allows lookup regardless of order while preserving direction
    """
    expected_dir = {}
    for query_code, target_code in species_pairs:
        if not query_code or not target_code:
            raise ValueError(f"Invalid species pair: ({query_code}, {target_code}). Both codes must be non-empty.")
        if query_code == target_code:
            raise ValueError(f"Invalid species pair: ({query_code}, {target_code}). Query and target must be different.")
        
        key = frozenset((query_code, target_code))
        if key in expected_dir:
            raise ValueError(f"Duplicate species pair detected: {query_code}, {target_code}")
        
        expected_dir[key] = (query_code, target_code)
    
    return expected_dir


def add_map_columns(df, prot2species, prot2gene=None):
    """Add Query/Target Species (+ code) and Gene columns to df."""
    df = df.copy()
    df["Query_Species"] = df["Query_Protein"].map(prot2species)
    df["Target_Species"] = df["Target_Protein"].map(prot2species)
    df["Q_code"] = df["Query_Species"].apply(species_code)
    df["T_code"] = df["Target_Species"].apply(species_code)

    if prot2gene is not None:
        df["Query_Gene"] = df["Query_Protein"].map(prot2gene)
        df["Target_Gene"] = df["Target_Protein"].map(prot2gene)

    return df


def canonicalize_direction_per_pair(df, expected_dir):
    """
    Canonicalize direction using expected pair direction PER PAIR (not global).
    Keeps original columns for auditing:
      Orig_Query_Protein, Orig_Target_Protein, Orig_Query_Species, Orig_Target_Species, ...
    Then swaps Query/Target when row is reversed vs expected_dir for that species pair.

    Adds:
      Pair = "<left>_<right>" according to expected direction
      Swapped = True/False
    """
    df = df.copy()

    # Preserve original for auditing
    df["Orig_Query_Protein"] = df["Query_Protein"]
    df["Orig_Target_Protein"] = df["Target_Protein"]
    df["Orig_Query_Species"] = df["Query_Species"]
    df["Orig_Target_Species"] = df["Target_Species"]
    df["Orig_Q_code"] = df["Q_code"]
    df["Orig_T_code"] = df["T_code"]
    if "Query_Gene" in df.columns and "Target_Gene" in df.columns:
        df["Orig_Query_Gene"] = df["Query_Gene"]
        df["Orig_Target_Gene"] = df["Target_Gene"]

    swapped = []

    for i, row in df.iterrows():
        qc, tc = row["Q_code"], row["T_code"]

        # If mapping missing, we can't canonicalize reliably.
        if not qc or not tc:
            df.at[i, "Pair"] = ""
            swapped.append(False)
            continue

        key = frozenset((qc, tc))
        if key not in expected_dir:
            # Not one of the expected species pairs
            df.at[i, "Pair"] = ""
            swapped.append(False)
            continue

        exp_q, exp_t = expected_dir[key]
        df.at[i, "Pair"] = f"{exp_q}_{exp_t}"

        # If current is reversed vs expected, swap.
        need_swap = (qc == exp_t and tc == exp_q)
        swapped.append(bool(need_swap))

        if need_swap:
            # swap proteins
            qprot, tprot = row["Query_Protein"], row["Target_Protein"]
            df.at[i, "Query_Protein"] = tprot
            df.at[i, "Target_Protein"] = qprot

            # swap species
            qsp, tsp = row["Query_Species"], row["Target_Species"]
            df.at[i, "Query_Species"] = tsp
            df.at[i, "Target_Species"] = qsp

            # swap codes
            df.at[i, "Q_code"] = exp_q
            df.at[i, "T_code"] = exp_t

            # swap genes (if present)
            if "Query_Gene" in df.columns and "Target_Gene" in df.columns:
                qg, tg = row["Query_Gene"], row["Target_Gene"]
                df.at[i, "Query_Gene"] = tg
                df.at[i, "Target_Gene"] = qg

    df["Swapped"] = swapped
    return df


def write_split_outputs(df, out_dir, unmapped_path, prefix):
    """
    Write:
      - per species-pair files (prefix_<pair>.tsv)
      - unmapped rows

    Output columns are kept minimal:
      Query_Protein, Target_Protein, Query_Species, Target_Species, Query_Gene, Target_Gene
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    unmapped = df[df["Pair"] == ""].copy()
    if len(unmapped) > 0:
        unmapped.to_csv(unmapped_path, sep="\t", index=False)
        print("Saved unmapped rows to:", unmapped_path)

    mappable = df[df["Pair"] != ""].copy()
    for pair_name, sub in mappable.groupby("Pair"):
        out_path = str(Path(out_dir) / f"{prefix}_{pair_name}.tsv")

        cols = [
            "Query_Protein", "Target_Protein",
            "Query_Species", "Target_Species",
            "Query_Gene", "Target_Gene"
        ]

        # If Gene columns don't exist (e.g. map lacks GeneID), fall back gracefully
        cols_existing = [c for c in cols if c in sub.columns]

        sub[cols_existing].to_csv(out_path, sep="\t", index=False)
        print("Saved:", out_path, "n=", len(sub))


def validate_canonical_direction(df, expected_dir, label):
    """
    Ensures every row with a known Pair is oriented exactly as expected.
    Raises AssertionError if any violations exist.
    """
    bad_rows = []

    for i, row in df.iterrows():
        if row.get("Pair", "") == "":
            continue  # unmapped / unknown pair
        qc, tc = row["Q_code"], row["T_code"]
        exp_q, exp_t = expected_dir[frozenset((qc, tc))]
        if (qc, tc) != (exp_q, exp_t):
            bad_rows.append((i, row["Query_Protein"], row["Target_Protein"], qc, tc, exp_q, exp_t))

    if bad_rows:
        print(f"\n[{label}] ERROR: {len(bad_rows)} row(s) not in canonical direction.")
        print("First 10:")
        for r in bad_rows[:10]:
            print("  idx=", r[0], "Q=", r[1], "T=", r[2], "got=", (r[3], r[4]), "expected=", (r[5], r[6]))
        raise AssertionError(f"{label}: canonical direction validation failed")

    print(f"[{label}] OK: all mappable rows are in canonical direction.")


def create_cactbuscomp_set_with_direction(ref_set, ref_dir, test_set, test_dir, expected_dir, prot2species):
    """
    Create CactBUSComp-set = union of BUSComp and Cactus pair sets with canonical direction.
    
    For each pair in the union:
    - If pair exists in ref (BUSComp), use ref's direction
    - Otherwise, use test's (Cactus) direction
    - Then canonicalize according to expected_dir
    
    Returns:
      - cactbuscomp_keys: set of undirected keys
      - cactbuscomp_directed: dict mapping undirected key -> canonical (query, target)
    """
    cactbuscomp_keys = ref_set | test_set
    cactbuscomp_directed = {}
    
    for key in cactbuscomp_keys:
        # Prefer ref direction if available, otherwise use test
        if key in ref_dir:
            pair = ref_dir[key]
        else:
            pair = test_dir[key]
        
        # Now canonicalize this pair according to expected directions
        q_prot, t_prot = pair
        
        # Get species for these proteins
        q_species = prot2species.get(q_prot, "")
        t_species = prot2species.get(t_prot, "")
        
        # Get species codes
        q_code = species_code(q_species)
        t_code = species_code(t_species)
        
        # Check if we need to swap
        if q_code and t_code:
            pair_key = frozenset((q_code, t_code))
            if pair_key in expected_dir:
                exp_q, exp_t = expected_dir[pair_key]
                # If current is reversed vs expected, swap
                if (q_code == exp_t and t_code == exp_q):
                    pair = (t_prot, q_prot)  # swap
        
        cactbuscomp_directed[key] = pair
    
    return cactbuscomp_keys, cactbuscomp_directed


def main():
    # --- Load data ---
    print("Loading BUSComp reference set...")
    ref_df = load_buscomp_pairs(REF_PATH)
    
    print("Loading Cactus test set...")
    test_df = load_pairs(TEST_PATH)

    # --- Convert to undirected keys, and also keep directed orientation per source ---
    ref_set, ref_dir = build_pair_index(ref_df, "Query_Protein", "Target_Protein")
    test_set, test_dir = build_pair_index(test_df, "Query_Protein", "Target_Protein")

    print(f"\nPairs in ref_set (BUSComp): {len(ref_set)}")
    print(f"Pairs in test_set (Cactus): {len(test_set)}")

    # --- Identify missed pairs (by undirected key) ---
    missed_by_test = ref_set - test_set   # BUSComp has, Cactus misses (FN)
    missed_by_ref = test_set - ref_set    # Cactus has, BUSComp misses (FP)

    # --- Save FN/FP with preserved direction ---
    # FN should be written in REF direction (BUSComp direction)
    fn_rows = [ref_dir[k] for k in missed_by_test]
    # FP should be written in TEST direction (Cactus direction)
    fp_rows = [test_dir[k] for k in missed_by_ref]

    Path(OUT_MISSED_BY_TEST).parent.mkdir(parents=True, exist_ok=True)
    fn_df = pd.DataFrame(fn_rows, columns=["Query_Protein", "Target_Protein"])
    fp_df = pd.DataFrame(fp_rows, columns=["Query_Protein", "Target_Protein"])

    fn_df.to_csv(OUT_MISSED_BY_TEST, sep="\t", index=False)
    fp_df.to_csv(OUT_MISSED_BY_REF, sep="\t", index=False)

    # --- Compute metrics ---
    m = compare_sets(ref_set, test_set)
    text = (
        f"Pairs in ref_set (BUSComp): {len(ref_set)}\n"
        f"Pairs in test_set (Cactus): {len(test_set)}\n\n"
        f"TP: {m['tp']}, FN: {m['fn']}, FP: {m['fp']}\n"
        f"Recall: {m['recall']:.3f}\n"
        f"Precision: {m['precision']:.3f}\n"
        f"F1 Score: {m['f1']:.3f}\n"
        f"Jaccard Similarity: {m['jaccard']:.3f}\n"
        f"Recovered {m['recall']:.1%} of reference 1:1 orthologs "
        f"({m['tp']} out of {m['tp'] + m['fn']})\n"
    )

    # --- Save metrics ---
    with open(OUT_METRICS, "w") as f:
        f.write(text)
    print("\n" + text)

    # Save metrics as TSV
    save_metrics_as_tsv(m, len(ref_set), len(test_set), OUT_METRICS_TSV)

    # ------------------------------------------------------------
    # Annotate FN/FP with species (+ gene) and split by pair
    # ------------------------------------------------------------

    # Check if BUSComp already has species columns
    buscomp_has_species = "Query_Species" in ref_df.columns and "Target_Species" in ref_df.columns
    
    if buscomp_has_species:
        print("\nBUSComp file already contains species information. Using it for mapping.")
        # Extract protein-to-species mapping from BUSComp
        prot2species_q = dict(zip(ref_df["Query_Protein"].astype(str), ref_df["Query_Species"].astype(str)))
        prot2species_t = dict(zip(ref_df["Target_Protein"].astype(str), ref_df["Target_Species"].astype(str)))
        prot2species = {**prot2species_q, **prot2species_t}
        
        # Check for gene columns
        if "Query_Gene" in ref_df.columns and "Target_Gene" in ref_df.columns:
            prot2gene_q = dict(zip(ref_df["Query_Protein"].astype(str), ref_df["Query_Gene"].astype(str)))
            prot2gene_t = dict(zip(ref_df["Target_Protein"].astype(str), ref_df["Target_Gene"].astype(str)))
            prot2gene = {**prot2gene_q, **prot2gene_t}
        else:
            prot2gene = None
    else:
        print("\nBUSComp file doesn't contain species information. Loading from mapping file.")
        map_df = pd.read_csv(MAP_PATH, sep="\t", dtype=str)
        if "ProteinID" not in map_df.columns or "Species" not in map_df.columns:
            raise ValueError("Map file must contain columns: ProteinID, Species")

        prot2species = dict(zip(map_df["ProteinID"].astype(str), map_df["Species"].astype(str)))

        prot2gene = None
        if "GeneID" in map_df.columns:
            prot2gene = (
                map_df.dropna(subset=["ProteinID"])
                      .drop_duplicates(subset=["ProteinID"], keep="first")
                      .set_index("ProteinID")["GeneID"]
                      .to_dict()
            )

    # Build expected directions from the configured species pairs
    expected_dir = build_expected_directions(EXPECTED_SPECIES_PAIRS)
    
    # Print expected directions once and verify
    print("\nExpected directions (from EXPECTED_SPECIES_PAIRS):")
    for q, t in EXPECTED_SPECIES_PAIRS:
        print(f"  {q} -> {t}")

    # Add mapping columns
    fn_df = add_map_columns(fn_df, prot2species, prot2gene)
    fp_df = add_map_columns(fp_df, prot2species, prot2gene)

    # Canonicalize direction PER species pair using expected directions
    fn_df = canonicalize_direction_per_pair(fn_df, expected_dir)
    fp_df = canonicalize_direction_per_pair(fp_df, expected_dir)

    # assert every mappable row matches the expected direction
    validate_canonical_direction(fn_df, expected_dir, "FN")
    validate_canonical_direction(fp_df, expected_dir, "FP")

    # Split and write per-pair files for FN and FP
    write_split_outputs(fn_df, OUT_SPLIT_DIR_FN, OUT_UNMAPPED_FN, prefix="false_negatives")
    write_split_outputs(fp_df, OUT_SPLIT_DIR_FP, OUT_UNMAPPED_FP, prefix="false_positives")

    # ------------------------------------------------------------
    # Create and save CactBUSComp-set with canonical direction
    # ------------------------------------------------------------
    
    print("\nCreating CactBUSComp set (union of BUSComp and Cactus)...")
    cactbuscomp_keys, cactbuscomp_directed = create_cactbuscomp_set_with_direction(
        ref_set, ref_dir, test_set, test_dir, expected_dir, prot2species
    )
    
    print(f"CactBUSComp set size: {len(cactbuscomp_keys)}")
    
    # Save as pickle (undirected keys for set operations)
    with open(OUT_CACTBUSCOMP_SET_PKL, "wb") as f:
        pickle.dump(cactbuscomp_keys, f)
    print(f"Saved CactBUSComp set (undirected keys) to: {OUT_CACTBUSCOMP_SET_PKL}")
    
    # Save as single TSV with canonical direction and full annotations
    cactbuscomp_pairs = [cactbuscomp_directed[k] for k in cactbuscomp_keys]
    cactbuscomp_df = pd.DataFrame(cactbuscomp_pairs, columns=["Query_Protein", "Target_Protein"])
    
    # Add species mapping to CactBUSComp set
    cactbuscomp_df = add_map_columns(cactbuscomp_df, prot2species, prot2gene)
    
    # Canonicalize CactBUSComp set direction
    cactbuscomp_df = canonicalize_direction_per_pair(cactbuscomp_df, expected_dir)
    
    # Validate CactBUSComp set direction
    validate_canonical_direction(cactbuscomp_df, expected_dir, "CactBUSComp")
    
    # Save as single comprehensive file with all annotations
    cols = [
        "Query_Protein", "Target_Protein",
        "Query_Species", "Target_Species",
    ]
    if "Query_Gene" in cactbuscomp_df.columns and "Target_Gene" in cactbuscomp_df.columns:
        cols += ["Query_Gene", "Target_Gene"]
    
    # Add the Pair column to help identify species pairs
    if "Pair" in cactbuscomp_df.columns:
        cols.append("Pair")
    
    cactbuscomp_df[cols].to_csv(OUT_CACTBUSCOMP_SET_TSV, sep="\t", index=False)
    print(f"Saved CactBUSComp set to: {OUT_CACTBUSCOMP_SET_TSV}")
    print(f"  Columns: {', '.join(cols)}")
    print(f"  Total pairs: {len(cactbuscomp_df)}")
    
    # Show summary by species pair
    if "Pair" in cactbuscomp_df.columns:
        print("\nCactBUSComp pairs by species combination:")
        pair_counts = cactbuscomp_df[cactbuscomp_df["Pair"] != ""].groupby("Pair").size().sort_index()
        for pair, count in pair_counts.items():
            print(f"  {pair}: {count} pairs")
        unmapped_count = len(cactbuscomp_df[cactbuscomp_df["Pair"] == ""])
        if unmapped_count > 0:
            print(f"  Unmapped: {unmapped_count} pairs")


if __name__ == "__main__":
    main()