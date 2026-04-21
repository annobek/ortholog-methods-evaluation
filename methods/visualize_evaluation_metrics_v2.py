import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


# ============================================================================
# CONFIGURATION - Specify exact paths to your TSV files
# ============================================================================

# Define methods (rows in heatmaps)
METHODS = [
    "Prot-syn (n-score)",
    "Prot-syn (sum)",
    "Tree-Majority",
    "Tree-Max-Score",
    "Tree-Whitelist",
    "PTP",
    "OrthoFinder"
]

# Define references (columns in heatmaps)
REFERENCES = ["Cactus", "CactBUSComp-set", "BUSComp-set"]

# Specify the path to each TSV file
# Structure: {method: {reference: "path/to/file.tsv"}}
TSV_PATHS = {
    "Prot-syn (n-score)": {
      "Cactus": "results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/cactus_prot_syn_nscore_metrics.tsv",
      "CactBUSComp-set": "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/nscore_approach/test_no_transitivity/CactBUSComp_prot_syn_nscore_metrics.tsv",
      "BUSComp-set": "results/evaluation/mammalia/vs_buscomp_set/vs_prot_syn/nscore_approach/test_no_transitivity/BUSComp_prot_syn_nscore_metrics.tsv"
    },
    "Prot-syn (sum)": {
      "Cactus": "results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/cactus_prot_syn_sum_metrics.tsv",
      "CactBUSComp-set": "results/evaluation/mammalia/all_vs_cactbuscomp/vs_prot_syn/sum_approach/test_no_transitivity/CactBUSComp_prot_syn_sum_metrics.tsv",
      "BUSComp-set": "results/evaluation/mammalia/vs_buscomp_set/vs_prot_syn/sum_approach/test_no_transitivity/BUSComp_prot_syn_sum_metrics.tsv"
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
    }
}

# Output settings
OUTPUT_FILE = "results/evaluation/mammalia/metrics_methods_vs_reference_sets_v2.png"
PLOT_TITLE = "Agreement between methods of 1:1 ortholog evaluation vs reference\n(Cactus, CactBUSComp-set and BUSComp-set)"

# ============================================================================
# MAIN CODE 
# ============================================================================

def read_metrics_from_tsv(filepath):
    """Read metrics from a TSV file."""
    try:
        df = pd.read_csv(filepath, sep='\t')
        metrics = {
            'Precision': df['Precision'].iloc[0],
            'Recall': df['Recall'].iloc[0],
            'F1_score': df['F1_score'].iloc[0],
            'Jaccard_Similarity': df['Jaccard_Similarity'].iloc[0]
        }
        return metrics
    except Exception as e:
        print(f"Warning: Could not read {filepath}: {e}")
        return None


def load_all_metrics(methods, references, tsv_paths):
    """Load metrics from specified paths."""
    # Initialize arrays for each metric
    precision = np.full((len(methods), len(references)), np.nan)
    recall = np.full((len(methods), len(references)), np.nan)
    f1_score = np.full((len(methods), len(references)), np.nan)
    jaccard = np.full((len(methods), len(references)), np.nan)
    
    # Read each file
    for i, method in enumerate(methods):
        for j, reference in enumerate(references):
            if method in tsv_paths and reference in tsv_paths[method]:
                filepath = Path(tsv_paths[method][reference])
                
                if filepath.exists():
                    metrics = read_metrics_from_tsv(filepath)
                    if metrics:
                        precision[i, j] = metrics['Precision']
                        recall[i, j] = metrics['Recall']
                        f1_score[i, j] = metrics['F1_score']
                        jaccard[i, j] = metrics['Jaccard_Similarity']
                        print(f" Loaded: {method} vs {reference}")
                else:
                    print(f" File not found: {filepath}")
            else:
                print(f" No path specified for: {method} vs {reference}")
    
    return {
        'Precision': precision,
        'Recall': recall,
        'F1-score': f1_score,
        'Jaccard': jaccard
    }


def create_heatmap_grid(metrics_data, methods, references, title, output_file):
    """Create 2x2 grid of heatmaps."""
    fig, axs = plt.subplots(2, 2, figsize=(13, 13))
    axs = axs.flatten()
    
    metric_names = ['Recall', 'Precision', 'F1-score', 'Jaccard']
    
    for ax, metric_name in zip(axs, metric_names):
        data = metrics_data[metric_name]
        
        # Create heatmap
        im = ax.imshow(data, cmap='viridis', vmin=0, vmax=1)
        
        # Set ticks and labels
        ax.set_xticks(range(len(references)))
        ax.set_xticklabels(references, rotation=45, ha="right")
        ax.set_yticks(range(len(methods)))
        ax.set_yticklabels(methods)
        
        # Add title
        ax.set_title(metric_name, fontsize=14)
        
        # Add values to cells
        for i in range(len(methods)):
            for j in range(len(references)):
                value = data[i, j]
                if not np.isnan(value):
                    ax.text(j, i, f"{value:.3f}",
                           ha="center", va="center", color="black", fontsize=10)
    
    # Add colorbar
    #fig.colorbar(im, ax=axs, orientation='vertical', 
    #            fraction=0.046, pad=0.04, label='Score')
    
    # Set main title
    fig.suptitle(title, fontsize=16, y=0.98)
    #plt.tight_layout(rect=[0, 0, 0.96, 0.96])
    plt.tight_layout()

    # Save figure
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    print(f"\n Saved figure to {output_path}")
    
    plt.show()
    plt.close()


# ============================================================================
# RUN
# ============================================================================

if __name__ == "__main__":
    print("Loading metrics from TSV files...")
    print("="*60)
    metrics = load_all_metrics(METHODS, REFERENCES, TSV_PATHS)
    
    print("\n" + "="*60)
    print("Creating heatmap visualization...")
    create_heatmap_grid(metrics, METHODS, REFERENCES, PLOT_TITLE, OUTPUT_FILE)
    
    print("\nDone!")