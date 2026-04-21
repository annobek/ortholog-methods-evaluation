#!/usr/bin/env python3
"""
Ready-to-Run Synteny Plot WITH LEGEND
Creates both the plot AND a text file showing what each number means
"""

from synteny_with_legend import plot_synteny_with_legend
import pandas as pd

# =============================================================================
# EDIT THIS SECTION WITH YOUR DATA
# =============================================================================

# Your synteny data files

# Forward direction

synteny_files = [
    # Dog vs others
    'results/genes_in_blocks/mammalia/region_gene_tables/can_mac_synteny_tables_regions_fix.tsv',
    'results/genes_in_blocks/mammalia/region_gene_tables/can_rat_synteny_tables_regions_fix.tsv',
    'results/genes_in_blocks/mammalia/region_gene_tables/can_mus_synteny_tables_regions_fix.tsv',
    
    # Macaque vs others
    'results/genes_in_blocks/mammalia/region_gene_tables/mac_rat_synteny_tables_regions_fix.tsv',
    'results/genes_in_blocks/mammalia/region_gene_tables/mac_mus_synteny_tables_regions_fix.tsv',
    
    # Rat vs mouse
    'results/genes_in_blocks/mammalia/region_gene_tables/rat_mus_synteny_tables_regions_fix.tsv',
]

# Species pairs (must match files above)
species_pairs = [
    ('canis_lupus_familiaris', 'macaca_fascicularis'),
    ('canis_lupus_familiaris', 'rattus_norvegicus'),
    ('canis_lupus_familiaris', 'mus_musculus'),
    ('macaca_fascicularis', 'rattus_norvegicus'),
    ('macaca_fascicularis', 'mus_musculus'),
    ('rattus_norvegicus', 'mus_musculus'),
]

# Output files
# forward direction
output_plot = 'results/linear_synteny_plots/linear_synteny_plot_forward.png'         # The plot image
output_legend = 'results/linear_synteny_plots/label_scaffold_legend_forward.tsv'       # The legend (number -> name)


# Filtering settings
TOP_SCAFFOLDS = 10
MIN_REGION_LENGTH = 2000000

# Visual settings  
SHOW_LABELS = True
COLOR_SCAFFOLDS = True
LABEL_FONTSIZE = 8

# =============================================================================
# RUN THE PLOT GENERATION
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("SYNTENY PLOT WITH LEGEND GENERATOR")
    print("=" * 70)
    print(f"\nWill create:")
    print(f"  1. {output_plot} - The synteny plot")
    print(f"  2. {output_legend} - Legend showing scaffold names")
    print(f"\nProcessing {len(synteny_files)} file(s)...")
    print("-" * 70)
    
    try:
        fig, ax, scaffold_map = plot_synteny_with_legend(
            data_files=synteny_files,
            species_pairs=species_pairs,
            output=output_plot,
            top_scaffolds=TOP_SCAFFOLDS,
            min_length=MIN_REGION_LENGTH,
            show_scaffold_labels=SHOW_LABELS,
            color_scaffolds=COLOR_SCAFFOLDS,
            label_fontsize=LABEL_FONTSIZE,
            show_legend=True,
            legend_file=output_legend
        )
        
        print("\n" + "=" * 70)
        print("SUCCESS!")
        print("=" * 70)
        print(f"\nGenerated files:")
        print(f"  ✓ {output_plot}")
        print(f"  ✓ {output_legend} (TSV format)")
        print(f"\nLegend format: Tab-separated with columns:")
        print(f"  - Label: Scaffold number (1, 2, 3, ...)")
        print(f"  - Species: Species name")
        print(f"  - Scaffold_Name: Actual chromosome/scaffold ID")
        print(f"  - Length_bp: Scaffold length in base pairs")
        print(f"\nNote: Labels are numbered within each species independently")
        print(f"      (e.g., both Dog and Mouse can have a 'Scaffold #1')")
        print("=" * 70)
        
        # Also print a quick preview of the legend
        print(f"\nQuick preview of {output_legend}:")
        print("-" * 70)
        legend_preview = pd.read_csv(output_legend, sep='\t')
        print(legend_preview.head(15).to_string(index=False))
        if len(legend_preview) > 15:
            print(f"\n... ({len(legend_preview) - 15} more rows)")
        
    except FileNotFoundError as e:
        print(f"\nERROR: Could not find file: {e}")
        print("Please check your file paths in the synteny_files list")
        
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        print("\nTroubleshooting:")
        print("1. Check file paths are correct")
        print("2. Verify TSV format with correct columns")
        print("3. Try increasing min_length or decreasing top_scaffolds")
