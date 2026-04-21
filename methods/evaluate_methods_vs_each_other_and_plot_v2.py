import pandas as pd
import itertools
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================================
# CONFIGURATION
# ============================================================================

# Methods to compare against each other (calculate Jaccard)
METHODS_CALCULATE = {
    "Prot-syn (n-score)": "material/sex_experiment/one_to_one_nscore.tsv",
    "Prot-syn (sum)": "material/sex_experiment/one_to_one_sum.tsv",
    "Tree-Majority": "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/majority/complete_tree_conserved_orthologs_2.tsv",
    "Tree-Max-Score": "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/standard/complete_tree_conserved_orthologs_2.tsv",
    "Tree-Whitelist": "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/whitelist/complete_tree_conserved_orthologs_2.tsv",
    "PTP": "/storage/EasyVectorOmics/phylotreepruner/results/mammalia/ptp_one_to_one_pairs_new.tsv",
    "OrthoFinder": "results/evaluation/mammalia/vs_orthofinder/orthofinder_1to1.tsv"
}

# Methods with pre-calculated Jaccard (extract from TSV files)
# Structure: {method_name: {"vs_method_name": "path/to/metrics.tsv"}}
METHODS_EXTRACT = {
    "Cactus": {
      "BUSComp-set": "results/evaluation/mammalia/vs_buscomp_set/cactus_BUSComp_metrics.tsv"
    },
    "Prot-syn (n-score)": {
      "Cactus": "results/evaluation/mammalia/vs_prot_syn/nscore_approach/cactus_prot_syn_nscore_metrics.tsv",
      "CactBUSComp-set": "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/nscore_approach/CactBUSComp_prot_syn_nscore_metrics.tsv",
      "BUSComp-set": "results/evaluation/mammalia/vs_buscomp_set/vs_prot_syn/nscore_approach/BUSComp_prot_syn_nscore_metrics.tsv"
    },
    "Prot-syn (sum)": {
      "Cactus": "results/evaluation/mammalia/vs_prot_syn/sum_approach/cactus_prot_syn_sum_metrics.tsv",
      "CactBUSComp-set": "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/sum_approach/CactBUSComp_prot_syn_sum_metrics.tsv",
      "BUSComp-set": "results/evaluation/mammalia/vs_buscomp_set/vs_prot_syn/sum_approach/BUSComp_prot_syn_sum_metrics.tsv"
    },
    "Tree-Majority": {
      "Cactus": "results/evaluation/mammalia/vs_tree/majority/cactus_tree_majority_metrics.tsv",
      "CactBUSComp-set": "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/majority/CactBUSComp_tree_majority_metrics.tsv",
      "BUSComp-set": "results/evaluation/mammalia/vs_buscomp_set/vs_tree/majority/BUSComp_tree_majority_metrics.tsv"
    },
    "Tree-Max-Score": {
      "Cactus": "results/evaluation/mammalia/vs_tree/max_score/cactus_tree_max_score_metrics.tsv",
      "CactBUSComp-set": "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/max_score/CactBUSComp_tree_max_score_metrics.tsv",
      "BUSComp-set": "results/evaluation/mammalia/vs_buscomp_set/vs_tree/max_score/BUSComp_tree_max_score_metrics.tsv"
    },
    "Tree-Whitelist": {
      "Cactus": "results/evaluation/mammalia/vs_tree/whitelist/cactus_tree_whitelist_metrics.tsv",
      "CactBUSComp-set": "results/evaluation/mammalia/all_vs_cactbuscomp/vs_tree/whitelist/CactBUSComp_tree_whitelist_metrics.tsv",
      "BUSComp-set": "results/evaluation/mammalia/vs_buscomp_set/vs_tree/whitelist/BUSComp_tree_whitelist_metrics.tsv"
    },
    "PTP": {
      "Cactus": "results/evaluation/mammalia/vs_ptp/cactus_ptp_metrics.tsv",
      "CactBUSComp-set": "results/evaluation/mammalia/all_vs_cactbuscomp/vs_ptp/CactBUSComp_ptp_metrics.tsv",
      "BUSComp-set": "results/evaluation/mammalia/vs_buscomp_set/vs_ptp/BUSComp_ptp_metrics.tsv"
    },
    "OrthoFinder": {
      "Cactus": "results/evaluation/mammalia/vs_orthofinder/cactus_orthofinder_metrics.tsv",
      "CactBUSComp-set": "results/evaluation/mammalia/all_vs_cactbuscomp/vs_orthofinder/CactBUSComp_orthofinder_metrics.tsv",
      "BUSComp-set": "results/evaluation/mammalia/vs_buscomp_set/vs_orthofinder/BUSComp_orthofinder_metrics.tsv"
    },
    "BUSCO": {
        "Compleasm": "/storage/EasyVectorOmics/busco/results/mammalia/evaluation/busco_compleasm_metrics.tsv"
    }
}

# Output settings
OUT_METRICS = "results/evaluation/mammalia/all_vs_all/combined_jaccard_matrix_v2_1.txt"
OUT_PLOT = "results/evaluation/mammalia/all_vs_all/combined_jaccard_heatmap_v2_1.png"
PLOT_TITLE = "Pairwise Jaccard Similarity Between Orthology Methods"

# ============================================================================
# FUNCTIONS
# ============================================================================

def normalize_id(x):
    """Make IDs comparable across files."""
    if pd.isna(x):
        return None
    x = str(x).strip()
    x = x.replace("protein", "")
    x = x.replace(" ", "_")
    return x


def load_ortholog_pairs(path, sep="\t"):
    """
    Load ortholog pairs from any supported file format.
    Returns standardized columns: ['Species1', 'Protein1', 'Species2', 'Protein2'].
    """
    df = pd.read_csv(path, sep=sep, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    cols = {c.lower(): c for c in df.columns}

    # Detect columns
    species_cols = [cols[c] for c in cols if "species" in c]
    prot_cols = [cols[c] for c in cols if "protein" in c]
    gene_cols = [cols[c] for c in cols if "gene" in c]

    # Prefer Protein columns
    if len(prot_cols) >= 2:
        use_prot = prot_cols
    elif len(gene_cols) >= 2:
        print(f"  No Protein columns found, using Gene columns instead")
        use_prot = gene_cols
    else:
        raise ValueError(f"Could not detect suitable ID columns. Found: {df.columns.tolist()}")

    # Create output
    df_out = df[[species_cols[0], use_prot[0], species_cols[1], use_prot[1]]].copy()
    df_out.columns = ["Species1", "Protein1", "Species2", "Protein2"]

    # Normalize
    df_out["Protein1"] = df_out["Protein1"].astype(str).str.strip()
    df_out["Protein2"] = df_out["Protein2"].astype(str).str.strip()

    # Drop invalid / self pairs
    df_out = df_out[df_out["Protein1"] != df_out["Protein2"]].drop_duplicates()

    return df_out


def load_pair_set(path):
    """Load ortholog pairs and return as a set of sorted tuples."""
    df = load_ortholog_pairs(path)
    return {tuple(sorted((p1, p2))) for p1, p2 in zip(df["Protein1"], df["Protein2"])}


def extract_jaccard_from_tsv(path):
    """Extract Jaccard similarity value from a metrics TSV file."""
    try:
        df = pd.read_csv(path, sep='\t')
        if 'Jaccard_Similarity' in df.columns:
            return df['Jaccard_Similarity'].iloc[0]
        else:
            print(f"  Warning: 'Jaccard_Similarity' column not found in {path}")
            return np.nan
    except Exception as e:
        print(f"  Warning: Could not read {path}: {e}")
        return np.nan


def calculate_jaccard(set_a, set_b):
    """Calculate Jaccard similarity between two sets."""
    inter = len(set_a & set_b)
    union = len(set_a | set_b)
    return inter / union if union > 0 else 0


def build_jaccard_matrix(methods_calc, methods_extract):
    """
    Build a complete Jaccard similarity matrix combining:
    - Calculated Jaccard for methods_calc (all vs all)
    - Extracted Jaccard from TSV files for methods_extract
    """
    # Step 1: Load pair sets for methods that need calculation
    print("Loading ortholog pairs for calculation...")
    pair_sets = {}
    for name, path in methods_calc.items():
        try:
            pair_sets[name] = load_pair_set(path)
            print(f"  ✓ {name}: {len(pair_sets[name]):,} pairs")
        except Exception as e:
            print(f"  ✗ {name}: Failed to load - {e}")
    
    # Step 2: Get all method names, including reference methods from METHODS_EXTRACT
    all_methods = list(methods_calc.keys())
    
    # Add methods from METHODS_EXTRACT keys (if not already present)
    for method in methods_extract.keys():
        if method not in all_methods:
            all_methods.append(method)
    
    # Add reference methods that appear in METHODS_EXTRACT comparisons
    reference_methods = set()
    for method, comparisons in methods_extract.items():
        for ref_method in comparisons.keys():
            reference_methods.add(ref_method)
    
    for ref_method in sorted(reference_methods):
        if ref_method not in all_methods:
            all_methods.append(ref_method)
    
    print(f"\nAll methods in matrix ({len(all_methods)}): {', '.join(all_methods)}")
    n = len(all_methods)
    
    # Step 3: Initialize matrix
    matrix = np.full((n, n), np.nan)
    np.fill_diagonal(matrix, 1.0)
    
    # Step 4: Calculate Jaccard for methods_calc (all vs all)
    print("\nCalculating pairwise Jaccard similarities...")
    calc_names = list(methods_calc.keys())
    for i, j in itertools.combinations(range(len(calc_names)), 2):
        name_i, name_j = calc_names[i], calc_names[j]
        if name_i in pair_sets and name_j in pair_sets:
            jacc = calculate_jaccard(pair_sets[name_i], pair_sets[name_j])
            idx_i = all_methods.index(name_i)
            idx_j = all_methods.index(name_j)
            matrix[idx_i, idx_j] = jacc
            matrix[idx_j, idx_i] = jacc
            print(f"  {name_i} vs {name_j}: {jacc:.3f}")
    
    # Step 5: Extract Jaccard from TSV files
    print("\nExtracting Jaccard from TSV files...")
    for method_a, comparisons in methods_extract.items():
        idx_a = all_methods.index(method_a)
        for method_b, tsv_path in comparisons.items():
            if method_b in all_methods:
                idx_b = all_methods.index(method_b)
                jacc = extract_jaccard_from_tsv(tsv_path)
                if not np.isnan(jacc):
                    matrix[idx_a, idx_b] = jacc
                    matrix[idx_b, idx_a] = jacc
                    print(f"  ✓ {method_a} vs {method_b}: {jacc:.3f}")
                else:
                    print(f"  ✗ {method_a} vs {method_b}: Failed to extract")
            else:
                print(f"  ✗ {method_a} vs {method_b}: {method_b} not in methods list")
    
    return matrix, all_methods


def save_matrix(matrix, method_names, output_file):
    """Save the Jaccard matrix to a text file."""
    df_jaccard = pd.DataFrame(matrix, index=method_names, columns=method_names)
    
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        f.write("=== Jaccard Similarity Matrix ===\n\n")
        f.write(df_jaccard.round(3).to_string())
        f.write("\n")
    
    print(f"\n✓ Saved matrix to {output_file}")


def plot_jaccard_heatmap(matrix, method_names, title, output_file, show_lower=False):
    """
    Create and save a Jaccard similarity heatmap.
    
    Parameters:
    -----------
    matrix : numpy array
        Jaccard similarity matrix
    method_names : list
        Names of methods (for axis labels)
    title : str
        Plot title
    output_file : str
        Path to save the figure
    show_lower : bool
        If False, hide lower triangle (since matrix is symmetric)
    """
    n = len(method_names)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create heatmap
    im = ax.imshow(matrix, cmap="coolwarm", vmin=0, vmax=1)
    
    # Set ticks and labels
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(method_names, rotation=45, ha="right")
    ax.set_yticklabels(method_names)
    
    # Optionally hide lower triangle
    if not show_lower:
        for i in range(n):
            for j in range(n):
                if i > j:
                    ax.add_patch(plt.Rectangle(
                        (j-0.5, i-0.5), 1, 1, 
                        color='white', ec='lightgray', linewidth=0.5
                    ))
    
    # Annotate with numbers
    for i in range(n):
        for j in range(n):
            if show_lower or i <= j:
                value = matrix[i, j]
                if not np.isnan(value):
                    ax.text(j, i, f"{value:.2f}", 
                           ha="center", va="center", color="black", fontsize=9)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label="Jaccard Similarity")
    
    # Set title
    ax.set_title(title, fontsize=14, pad=20)
    
    # Save
    plt.tight_layout()
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved heatmap to {output_file}")
    
    plt.show()
    plt.close()


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    print("="*70)
    print("Building Combined Jaccard Similarity Matrix")
    print("="*70)
    
    # Build the matrix
    matrix, method_names = build_jaccard_matrix(METHODS_CALCULATE, METHODS_EXTRACT)
    
    # Save matrix
    save_matrix(matrix, method_names, OUT_METRICS)
    
    # Create visualization
    print("\nCreating heatmap...")
    plot_jaccard_heatmap(
        matrix, 
        method_names, 
        PLOT_TITLE, 
        OUT_PLOT,
        show_lower=False  # Set to True if you want to show both triangles
    )
    
    print("\n" + "="*70)
    print("Done!")
    print("="*70)
