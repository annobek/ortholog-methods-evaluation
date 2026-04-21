"""
Neighborhood Length Analysis Script
====================================

PURPOSE:
--------
This script analyzes the genomic length of gene neighborhoods in synteny analysis,
comparing three categories:
1. ALL neighborhoods (baseline from all genes)
2. FALSE NEGATIVES (pairs found by Cactus but missed/contradicted by ProtSyn)
3. FALSE POSITIVES (pairs found by ProtSyn but not by Cactus)

LOGIC:
------
For each gene, we:
1. Look up its neighborhood (list of neighboring gene IDs from neighborhoods.pkl)
2. Get genomic coordinates for all genes in the neighborhood from GTF files
3. Calculate neighborhood span: max(gene_end) - min(gene_start) in kilobases

For pair-level data (FN and FP):
- Calculate neighborhood length for BOTH genes in the orthologous pair
- Take the MEAN of the two lengths as the pair's neighborhood length

INPUTS:
-------
--fn                : False negatives gene-level table (TSV)
                      Columns: FN_ID, Pair, Focus_Gene, Cactus_Ortholog_Gene, 
                               FN_Gene_Position, etc.
--fp                : False positives table (TSV)
                      Columns: Gene1, Gene2_prot_syn, Pair, Species1, Species2
--neighborhoods     : Pickle file with gene neighborhoods
                      Format: {gene_id: [neighbor1, neighbor2, ...]}
--gtf-can/rat/mus/mac : GTF annotation files for each species
                        Must contain: gene entries with db_xref "GeneID:XXXXX"

OUTPUTS:
--------
--output-plot       : Raincloud plot (PNG) showing length distributions per species pair
--output-summary    : Summary statistics (TSV) - mean, median, std per category/pair
--output-detailed   : Full dataset (TSV) with all calculated lengths

USAGE EXAMPLE:
--------------
python neighborhood_analysis.py \
    --fn false_negatives_gene_level.tsv \
    --fp false_positives.tsv \
    --neighborhoods neighborhood.pkl \
    --gtf-can canis.gtf \
    --gtf-rat rattus.gtf \
    --gtf-mus mus.gtf \
    --gtf-mac macaca.gtf \
    --output-plot results/neighborhood_plot.png \
    --output-summary results/summary.tsv
"""

import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import re
import argparse
import sys

# ============================================================================
# STEP 1: Parse GTF and build gene coordinate dictionary
# ============================================================================

def extract_gene_id(attributes_str):
    """
    Extract numeric GeneID from GTF attributes column.
    
    The function searches for the pattern db_xref "GeneID:XXXXX" in the 
    attributes string and extracts the numeric ID. Handles cases where 
    multiple db_xref entries exist.
    
    Parameters:
    -----------
    attributes_str : str
        The attributes column from a GTF file (9th column)
        Example: 'gene_id "P2rx5"; db_xref "GeneID:113995"; db_xref "RGD:620256"'
    
    Returns:
    --------
    str or None
        The numeric GeneID as a string (e.g., "113995"), or None if not found
    
    Example:
    --------
    >>> attr = 'gene_id "P2rx5"; db_xref "GeneID:113995"; db_xref "RGD:620256"'
    >>> extract_gene_id(attr)
    '113995'
    """
    # Search for all occurrences of db_xref "GeneID:digits"
    # The pattern captures only the numeric part after GeneID:
    matches = re.findall(r'db_xref\s+"GeneID:(\d+)"', attributes_str)
    
    if matches:
        # Return the first GeneID found
        # Note: There should typically be only one GeneID per gene entry
        return matches[0]
    
    return None


def build_gene_coordinate_dict(gtf_files):
    """
    Build a lookup dictionary mapping gene IDs to their genomic coordinates.
    
    Parses GTF files for multiple species and creates a unified dictionary
    with gene locations. Only processes entries where feature == 'gene'.
    
    Parameters:
    -----------
    gtf_files : dict
        Dictionary mapping species codes to GTF file paths
        Format: {'species_code': 'path/to/file.gtf'}
        Example: {'can': '/data/canis.gtf', 'rat': '/data/rattus.gtf'}
    
    Returns:
    --------
    dict
        Nested dictionary with structure:
        {
            'gene_id': {
                'species': 'species_code',
                'chrom': 'chromosome_name',
                'start': int (genomic start position),
                'end': int (genomic end position)
            }
        }
        Example: {'113995': {'species': 'rat', 'chrom': 'NC_086028.1', 
                             'start': 58276329, 'end': 58291108}}
    
    Side Effects:
    -------------
    Prints progress information to stdout including:
    - Number of gene entries found per species
    - Number of GeneIDs extracted per species
    - Sample GeneIDs from each species
    - Warning if any genes lack GeneID annotations
    """
    gene_coords = {}
    genes_without_id = 0
    
    for species_code, gtf_file in gtf_files.items():
        print(f"Processing GTF for {species_code}...")
        
        # Read GTF file
        # GTF format: 9 tab-separated columns, lines starting with # are comments
        df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None,
                         names=['chrom', 'source', 'feature', 'start', 'end', 
                                'score', 'strand', 'frame', 'attributes'])
        
        # Filter only gene entries (not transcripts, exons, etc.)
        df_genes = df[df['feature'] == 'gene'].copy()
        
        print(f"  Found {len(df_genes)} gene entries")
        
        species_count = 0
        for _, row in df_genes.iterrows():
            gene_id = extract_gene_id(row['attributes'])
            
            if gene_id:
                # Store gene coordinates with species information
                gene_coords[gene_id] = {
                    'species': species_code,
                    'chrom': row['chrom'],
                    'start': row['start'],
                    'end': row['end']
                }
                species_count += 1
            else:
                genes_without_id += 1
        
        print(f"  Extracted {species_count} numeric GeneIDs")
        
        # Debug: print sample gene IDs to verify format
        species_genes = [k for k, v in gene_coords.items() if v['species'] == species_code]
        if len(species_genes) > 0:
            sample_ids = species_genes[:5]
            print(f"  Sample GeneIDs from {species_code}: {sample_ids}")
    
    if genes_without_id > 0:
        print(f"\nWarning: {genes_without_id} gene entries had no GeneID in db_xref")
    
    print(f"\nTotal genes in dictionary: {len(gene_coords)}")
    return gene_coords


# ============================================================================
# STEP 2: Calculate neighborhood length
# ============================================================================

def calc_neighborhood_length(gene_id, neighborhoods, gene_coords):
    """
    Calculate the genomic span of a gene's neighborhood.
    
    Takes a gene and its neighborhood (list of neighboring genes), looks up
    the genomic coordinates of all genes in the neighborhood, and calculates
    the total genomic distance they span.
    
    Parameters:
    -----------
    gene_id : str
        The focal gene ID (as string, e.g., "113995")
    neighborhoods : dict
        Dictionary mapping gene IDs to lists of neighbor gene IDs
        Format: {gene_id: [neighbor1_id, neighbor2_id, ...]}
        Example: {'113995': ['448785', '475235', '84024']}
    gene_coords : dict
        Dictionary from build_gene_coordinate_dict() containing genomic coordinates
    
    Returns:
    --------
    float or None
        Neighborhood length in kilobases (kb), or None if:
        - Gene not found in neighborhoods dict
        - No coordinates available for neighborhood genes
        
    Algorithm:
    ----------
    1. Look up neighborhood genes for focal gene
    2. Get genomic coordinates for all genes in neighborhood
    3. Include focal gene's coordinates if available
    4. Calculate span: (max end position - min start position) / 1000
    
    Example:
    --------
    If a neighborhood contains 3 genes at positions:
    Gene A: chr1:1000-2000
    Gene B: chr1:5000-6000  
    Gene C: chr1:8000-9000
    Length = (9000 - 1000) / 1000 = 8.0 kb
    """
    # Check if gene has a neighborhood
    if gene_id not in neighborhoods:
        return None
    
    neighbor_genes = neighborhoods[gene_id]
    
    # Collect coordinates for all genes in the neighborhood
    coords = []
    
    # Get coordinates for each neighbor
    for ngene in neighbor_genes:
        if ngene in gene_coords:
            coords.append((gene_coords[ngene]['start'], gene_coords[ngene]['end']))
    
    # Also include the focal gene itself
    if gene_id in gene_coords:
        coords.append((gene_coords[gene_id]['start'], gene_coords[gene_id]['end']))
    
    # Return None if no coordinates found
    if not coords:
        return None
    
    # Calculate genomic span in kilobases
    min_start = min(c[0] for c in coords)
    max_end = max(c[1] for c in coords)
    length_kb = (max_end - min_start) / 1000
    
    return length_kb


# ============================================================================
# STEP 3: Process all data categories
# ============================================================================

def process_all_neighborhoods(neighborhoods, gene_coords):
    """
    Calculate neighborhood lengths for ALL genes in the neighborhoods dictionary.
    
    This provides the baseline distribution of neighborhood lengths across
    all genes, regardless of whether they're in FN or FP sets.
    
    Parameters:
    -----------
    neighborhoods : dict
        Dictionary mapping gene IDs to neighbor lists
        Format: {gene_id: [neighbor1, neighbor2, ...]}
    gene_coords : dict
        Dictionary from build_gene_coordinate_dict()
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with columns:
        - gene_id: str, the gene identifier
        - length_kb: float, neighborhood length in kilobases
        - species: str, species code (can/rat/mus/mac)
        - category: str, always 'all' for this function
        
    Side Effects:
    -------------
    Prints diagnostic information:
    - Number of genes in each input dictionary
    - Overlap between neighborhoods and GTF coordinates
    - Sample gene IDs to help diagnose mismatches
    - Warning if no overlap found (indicates ID format mismatch)
    """
    print("\n=== Processing ALL neighborhoods ===")
    
    # Debug: check overlap between neighborhoods and gene_coords
    neighborhood_genes = set(neighborhoods.keys())
    coord_genes = set(gene_coords.keys())
    overlap = neighborhood_genes & coord_genes
    
    print(f"  Genes in neighborhoods dict: {len(neighborhood_genes)}")
    print(f"  Genes in GTF coords dict: {len(coord_genes)}")
    print(f"  Overlap: {len(overlap)}")
    
    # Sample some IDs to see format
    print(f"  Sample neighborhood gene IDs: {list(neighborhood_genes)[:5]}")
    print(f"  Sample GTF gene IDs: {list(coord_genes)[:5]}")
    
    if len(overlap) == 0:
        print("  ERROR: No overlap between neighborhood and GTF gene IDs!")
        print("  This suggests a gene ID format mismatch.")
        return pd.DataFrame()
    
    all_lengths = []
    
    # Calculate length for each gene's neighborhood
    for gene_id in neighborhoods.keys():
        length = calc_neighborhood_length(gene_id, neighborhoods, gene_coords)
        
        if length is not None:
            # Determine species from gene_coords
            species = gene_coords.get(gene_id, {}).get('species', 'unknown')
            
            all_lengths.append({
                'gene_id': gene_id,
                'length_kb': length,
                'species': species,
                'category': 'all'
            })
    
    print(f"  Calculated {len(all_lengths)} neighborhood lengths")
    return pd.DataFrame(all_lengths)


def process_false_negatives(fn_gene_level, neighborhoods, gene_coords):
    """
    Calculate mean neighborhood lengths for false negative gene pairs.
    
    False negatives are pairs where Cactus found a relationship but ProtSyn
    either found nothing or a contradictory result. For each pair, we calculate
    the neighborhood length for BOTH genes and take their mean.
    
    Parameters:
    -----------
    fn_gene_level : pd.DataFrame
        Gene-level false negatives table with columns:
        - FN_ID: int, unique identifier for each FN pair
        - Pair: str, species pair code (e.g., 'can_rat')
        - Focus_Gene: int/str, gene ID in first species
        - Cactus_Ortholog_Gene: int/str, ortholog gene ID in second species
        - FN_Gene_Position: str, either 'Gene1' or 'Gene2'
        
        Each FN_ID should appear in exactly 2 rows (one for each gene in the pair)
    
    neighborhoods : dict
        Gene neighborhoods dictionary
    gene_coords : dict
        Gene coordinates dictionary
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with columns:
        - fn_id: int, false negative pair identifier
        - pair: str, species pair code
        - focal_gene: str, gene ID from first species
        - partner_gene: str, gene ID from second species
        - focal_length_kb: float, neighborhood length of focal gene
        - partner_length_kb: float, neighborhood length of partner gene
        - mean_length_kb: float, mean of the two lengths
        - category: str, always 'false_negative'
    
    Side Effects:
    -------------
    Prints:
    - Sample gene IDs to verify format
    - Number of successfully calculated pair lengths
    - Warnings if no pairs calculated with diagnostic info
    """
    print("\n=== Processing FALSE NEGATIVES ===")
    
    # Debug: check gene ID formats in the input data
    sample_focal = fn_gene_level['Focus_Gene'].head(5).tolist()
    sample_cactus = fn_gene_level['Cactus_Ortholog_Gene'].head(5).tolist()
    print(f"  Sample Focus_Gene IDs: {sample_focal}")
    print(f"  Sample Cactus_Ortholog_Gene IDs: {sample_cactus}")
    
    fn_lengths = []
    
    # Group by FN_ID to process pairs
    for fn_id, group in fn_gene_level.groupby('FN_ID'):
        # Each FN_ID should have exactly 2 rows (Gene1 and Gene2 of the pair)
        if len(group) != 2:
            continue
        
        pair = group['Pair'].iloc[0]
        
        # Get both genes from the pair
        # Gene1 row contains focal gene in Focus_Gene and partner in Cactus_Ortholog_Gene
        gene1_row = group[group['FN_Gene_Position'] == 'Gene1'].iloc[0]
        
        focal_gene = str(gene1_row['Focus_Gene'])
        partner_gene = str(gene1_row['Cactus_Ortholog_Gene'])
        
        # Calculate neighborhood lengths for both genes
        focal_length = calc_neighborhood_length(focal_gene, neighborhoods, gene_coords)
        partner_length = calc_neighborhood_length(partner_gene, neighborhoods, gene_coords)
        
        # Only include if both lengths calculated successfully
        if focal_length is not None and partner_length is not None:
            mean_length = (focal_length + partner_length) / 2
            
            fn_lengths.append({
                'fn_id': fn_id,
                'pair': pair,
                'focal_gene': focal_gene,
                'partner_gene': partner_gene,
                'focal_length_kb': focal_length,
                'partner_length_kb': partner_length,
                'mean_length_kb': mean_length,
                'category': 'false_negative'
            })
    
    print(f"  Calculated {len(fn_lengths)} FN pair lengths")
    
    # Diagnostic output if no pairs were processed
    if len(fn_lengths) == 0:
        print("  WARNING: No FN lengths calculated!")
        print("  Checking if genes exist in neighborhoods and coords...")
        focal_in_neighborhoods = sum(1 for g in fn_gene_level['Focus_Gene'].astype(str).unique() 
                                     if g in neighborhoods)
        focal_in_coords = sum(1 for g in fn_gene_level['Focus_Gene'].astype(str).unique() 
                              if g in gene_coords)
        print(f"    Focal genes in neighborhoods: {focal_in_neighborhoods}")
        print(f"    Focal genes in coords: {focal_in_coords}")
    
    return pd.DataFrame(fn_lengths)


def process_false_positives(fp_df, neighborhoods, gene_coords):
    """
    Calculate mean neighborhood lengths for false positive gene pairs.
    
    False positives are pairs where ProtSyn found a relationship but Cactus
    found nothing. For each pair, we calculate the neighborhood length for
    BOTH genes and take their mean.
    
    Parameters:
    -----------
    fp_df : pd.DataFrame
        False positives table with columns:
        - Gene1: int/str, gene ID in first species
        - Gene2_prot_syn: int/str, gene ID in second species (found by ProtSyn)
        - Pair: str, species pair code (e.g., 'can_rat')
        - Species1: str, full species name for Gene1
        - Species2_prot_syn: str, full species name for Gene2
    
    neighborhoods : dict
        Gene neighborhoods dictionary
    gene_coords : dict
        Gene coordinates dictionary
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with columns:
        - gene1: str, gene ID from first species
        - gene2: str, gene ID from second species
        - pair: str, species pair code
        - gene1_length_kb: float, neighborhood length of gene1
        - gene2_length_kb: float, neighborhood length of gene2
        - mean_length_kb: float, mean of the two lengths
        - category: str, always 'false_positive'
    
    Side Effects:
    -------------
    Prints:
    - Sample gene IDs to verify format
    - Number of successfully calculated pair lengths
    - Warnings if no pairs calculated with diagnostic info
    """
    print("\n=== Processing FALSE POSITIVES ===")
    
    # Debug: check gene ID formats
    sample_gene1 = fp_df['Gene1'].head(5).tolist()
    sample_gene2 = fp_df['Gene2_prot_syn'].head(5).tolist()
    print(f"  Sample Gene1 IDs: {sample_gene1}")
    print(f"  Sample Gene2_prot_syn IDs: {sample_gene2}")
    
    fp_lengths = []
    
    for idx, row in fp_df.iterrows():
        gene1 = str(row['Gene1'])
        gene2 = str(row['Gene2_prot_syn'])
        pair = row['Pair']
        
        # Calculate neighborhood lengths for both genes
        gene1_length = calc_neighborhood_length(gene1, neighborhoods, gene_coords)
        gene2_length = calc_neighborhood_length(gene2, neighborhoods, gene_coords)
        
        # Only include if both lengths calculated successfully
        if gene1_length is not None and gene2_length is not None:
            mean_length = (gene1_length + gene2_length) / 2
            
            fp_lengths.append({
                'gene1': gene1,
                'gene2': gene2,
                'pair': pair,
                'gene1_length_kb': gene1_length,
                'gene2_length_kb': gene2_length,
                'mean_length_kb': mean_length,
                'category': 'false_positive'
            })
    
    print(f"  Calculated {len(fp_lengths)} FP pair lengths")
    
    # Diagnostic output if no pairs were processed
    if len(fp_lengths) == 0:
        print("  WARNING: No FP lengths calculated!")
        print("  Checking if genes exist in neighborhoods and coords...")
        gene1_in_neighborhoods = sum(1 for g in fp_df['Gene1'].astype(str).unique() 
                                     if g in neighborhoods)
        gene1_in_coords = sum(1 for g in fp_df['Gene1'].astype(str).unique() 
                              if g in gene_coords)
        print(f"    Gene1 in neighborhoods: {gene1_in_neighborhoods}")
        print(f"    Gene1 in coords: {gene1_in_coords}")
    
    return pd.DataFrame(fp_lengths)


# ============================================================================
# STEP 4: Create visualizations
# ============================================================================

def plot_raincloud(data_df, output_file='neighborhood_lengths_raincloud.png'):
    """
    Create a raincloud plot showing neighborhood length distributions.
    
    A raincloud plot combines a violin plot (distribution shape), box plot
    (quartiles), and individual data points to show the full distribution
    of neighborhood lengths across three categories for each species pair.
    
    Parameters:
    -----------
    data_df : pd.DataFrame
        Combined data with columns:
        - length_kb: float, neighborhood length in kilobases
        - category: str, one of ['all', 'false_negative', 'false_positive']
        - pair: str, species pair code (e.g., 'can_rat')
    
    output_file : str, optional
        Path where the plot PNG will be saved
        Default: 'neighborhood_lengths_raincloud.png'
    
    Returns:
    --------
    None
    
    Outputs:
    --------
    - Saves a high-resolution PNG plot to output_file
    - Creates one subplot per species pair, stacked vertically
    - Uses color coding: gray (all), coral (FN), skyblue (FP)
    
    Side Effects:
    -------------
    - Prints message when plot is saved
    - Prints warning if ptitprince library not available (uses violin plot fallback)
    - Closes the matplotlib figure after saving
    
    Notes:
    ------
    Requires ptitprince library for true raincloud plots. If not available,
    falls back to violin + box plot combination.
    Install with: pip install ptitprince
    """
    if len(data_df) == 0:
        print("ERROR: No data to plot!")
        return
    
    # Check if ptitprince library is available for raincloud plots
    try:
        import ptitprince as pt
        use_raincloud = True
    except ImportError:
        print("ptitprince not available, using violin plots instead")
        print("Install with: pip install ptitprince")
        use_raincloud = False
    
    # Get unique species pairs and sort for consistent ordering
    pairs = sorted(data_df['pair'].unique())
    n_pairs = len(pairs)
    
    # Create figure with one subplot per species pair
    fig, axes = plt.subplots(n_pairs, 1, figsize=(12, 4*n_pairs), sharex=True)
    if n_pairs == 1:
        axes = [axes]
    
    # Define color palette for categories
    # Create a list of colors in the same order as category_order
    category_order = ['all', 'false_negative', 'false_positive']
    color_list = ['lightgray', 'coral', 'skyblue']
    
    # Create plot for each species pair
    for idx, pair in enumerate(pairs):
        ax = axes[idx]
        pair_data = data_df[data_df['pair'] == pair]
        
        if use_raincloud:
            # Full raincloud plot with distribution, box, and points
            pt.RainCloud(
                data=pair_data,
                x='category',
                y='length_kb',
                order=category_order,
                palette=color_list,  # Use list instead of dict
                ax=ax,
                orient='h',
                width_viol=0.6,
                width_box=0.2
            )
        else:
            # Fallback: violin + box plot
            # Create a proper palette dict for seaborn
            palette_dict = dict(zip(category_order, color_list))
            
            sns.violinplot(
                data=pair_data,
                y='category',
                x='length_kb',
                order=category_order,
                palette=palette_dict,
                ax=ax,
                orient='h'
            )
            # Overlay box plot for quartiles
            sns.boxplot(
                data=pair_data,
                y='category',
                x='length_kb',
                order=category_order,
                ax=ax,
                width=0.3,
                boxprops={'facecolor': 'none'},
                orient='h'
            )
        
        # Customize subplot
        ax.set_title(f'Species Pair: {pair}', fontsize=14, fontweight='bold')
        ax.set_xlabel('Neighborhood Length (kb)', fontsize=12)
        ax.set_ylabel('Category', fontsize=12)
        ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to {output_file}")
    plt.close()


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

def parse_args():
    """
    Parse command-line arguments.
    
    Returns:
    --------
    argparse.Namespace
        Object containing all parsed arguments as attributes
    """
    parser = argparse.ArgumentParser(
        description='Analyze neighborhood lengths for false negatives and false positives',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python neighborhood_analysis.py \\
    --fn false_negatives_gene_level.tsv \\
    --fp false_positives.tsv \\
    --neighborhoods neighborhood.pkl \\
    --gtf-can canis.gtf \\
    --gtf-rat rattus.gtf \\
    --gtf-mus mus.gtf \\
    --gtf-mac macaca.gtf \\
    --output-plot results/neighborhood_plot.png \\
    --output-summary results/neighborhood_summary.tsv
        """
    )
    
    # Input files
    parser.add_argument('--fn', required=True,
                        help='Path to false negatives gene-level file (TSV)')
    parser.add_argument('--fp', required=True,
                        help='Path to false positives file (TSV)')
    parser.add_argument('--neighborhoods', required=True,
                        help='Path to neighborhoods pickle file')
    
    # GTF files for each species
    parser.add_argument('--gtf-can', required=True,
                        help='Path to Canis lupus familiaris GTF')
    parser.add_argument('--gtf-rat', required=True,
                        help='Path to Rattus norvegicus GTF')
    parser.add_argument('--gtf-mus', required=True,
                        help='Path to Mus musculus GTF')
    parser.add_argument('--gtf-mac', required=True,
                        help='Path to Macaca fascicularis GTF')
    
    # Output files
    parser.add_argument('--output-plot', default='neighborhood_lengths_raincloud.png',
                        help='Output path for raincloud plot (default: neighborhood_lengths_raincloud.png)')
    parser.add_argument('--output-summary', default='neighborhood_length_summary.tsv',
                        help='Output path for summary statistics TSV (default: neighborhood_length_summary.tsv)')
    parser.add_argument('--output-detailed', default='neighborhood_lengths_detailed.tsv',
                        help='Output path for detailed results TSV (default: neighborhood_lengths_detailed.tsv)')
    
    return parser.parse_args()


# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def main():
    """
    Main analysis workflow.
    
    Orchestrates the complete analysis pipeline:
    1. Parse command-line arguments
    2. Load input data (FN, FP, neighborhoods, GTFs)
    3. Build gene coordinate dictionary from GTFs
    4. Process all three categories (all, FN, FP)
    5. Combine data for visualization
    6. Generate plots and save results
    
    Returns:
    --------
    None
    
    Outputs:
    --------
    - Raincloud plot (PNG)
    - Summary statistics (TSV)
    - Detailed data (TSV)
    
    Exits:
    ------
    Returns early if no data can be calculated (with error messages)
    """
    args = parse_args()
    
    # Define GTF files from command-line arguments
    gtf_files = {
        'can': args.gtf_can,
        'rat': args.gtf_rat,
        'mus': args.gtf_mus,
        'mac': args.gtf_mac
    }
    
    # ========================================================================
    # LOAD DATA
    # ========================================================================
    print("="*80)
    print("NEIGHBORHOOD LENGTH ANALYSIS")
    print("="*80)
    print("\nLoading data...")
    print(f"  FN file: {args.fn}")
    print(f"  FP file: {args.fp}")
    print(f"  Neighborhoods: {args.neighborhoods}")
    
    neighborhoods = pickle.load(open(args.neighborhoods, 'rb'))
    print(f"  Loaded {len(neighborhoods)} neighborhoods")
    
    fn_gene_level = pd.read_csv(args.fn, sep='\t')
    print(f"  Loaded {len(fn_gene_level)} FN gene-level records")
    
    fp_df = pd.read_csv(args.fp, sep='\t')
    print(f"  Loaded {len(fp_df)} FP records")
    
    # ========================================================================
    # BUILD GENE COORDINATE DICTIONARY
    # ========================================================================
    gene_coords = build_gene_coordinate_dict(gtf_files)
    
    # ========================================================================
    # PROCESS ALL CATEGORIES
    # ========================================================================
    all_df = process_all_neighborhoods(neighborhoods, gene_coords)
    fn_df = process_false_negatives(fn_gene_level, neighborhoods, gene_coords)
    fp_df_processed = process_false_positives(fp_df, neighborhoods, gene_coords)
    
    # Check if we have any data at all
    if len(all_df) == 0 and len(fn_df) == 0 and len(fp_df_processed) == 0:
        print("\n" + "="*80)
        print("ERROR: No neighborhood lengths calculated!")
        print("="*80)
        print("This suggests a gene ID format mismatch.")
        print("Please check the debug output above.")
        sys.exit(1)
    
    # ========================================================================
    # PREPARE DATA FOR PLOTTING
    # ========================================================================
    print("\n=== Preparing data for visualization ===")
    
    # For "all": use individual gene neighborhood lengths
    if len(all_df) > 0:
        all_plot_df = all_df[['length_kb', 'species']].copy()
        all_plot_df['category'] = 'all'
    else:
        all_plot_df = pd.DataFrame(columns=['length_kb', 'species', 'category'])
    
    # Determine which species pairs we have data for
    pairs_in_data = set()
    if len(fn_df) > 0:
        pairs_in_data.update(fn_df['pair'].unique())
    if len(fp_df_processed) > 0:
        pairs_in_data.update(fp_df_processed['pair'].unique())
    
    if len(pairs_in_data) == 0:
        print("ERROR: No species pairs found in FN or FP data!")
        sys.exit(1)
    
    print(f"  Found species pairs: {sorted(pairs_in_data)}")
    
    # Replicate "all" data for each species pair
    # (since "all" is baseline for comparison with each pair)
    all_expanded = []
    for pair in pairs_in_data:
        temp = all_plot_df.copy()
        temp['pair'] = pair
        all_expanded.append(temp)
    
    if len(all_expanded) > 0:
        all_plot_df = pd.concat(all_expanded, ignore_index=True)
    else:
        all_plot_df = pd.DataFrame(columns=['length_kb', 'category', 'pair'])
    
    # For FN and FP: use mean lengths (already calculated)
    if len(fn_df) > 0:
        fn_plot_df = fn_df[['mean_length_kb', 'pair']].copy()
        fn_plot_df.rename(columns={'mean_length_kb': 'length_kb'}, inplace=True)
        fn_plot_df['category'] = 'false_negative'
    else:
        fn_plot_df = pd.DataFrame(columns=['length_kb', 'pair', 'category'])
    
    if len(fp_df_processed) > 0:
        fp_plot_df = fp_df_processed[['mean_length_kb', 'pair']].copy()
        fp_plot_df.rename(columns={'mean_length_kb': 'length_kb'}, inplace=True)
        fp_plot_df['category'] = 'false_positive'
    else:
        fp_plot_df = pd.DataFrame(columns=['length_kb', 'pair', 'category'])
    
    # Combine all data for plotting
    plot_data = pd.concat([
        all_plot_df[['length_kb', 'category', 'pair']],
        fn_plot_df,
        fp_plot_df
    ], ignore_index=True)
    
    print(f"\n  Total data points for plotting: {len(plot_data)}")
    print(f"    - All: {len(all_plot_df)}")
    print(f"    - FN: {len(fn_plot_df)}")
    print(f"    - FP: {len(fp_plot_df)}")
    
    if len(plot_data) == 0:
        print("\nERROR: No data to plot!")
        sys.exit(1)
    
    # ========================================================================
    # CREATE VISUALIZATION
    # ========================================================================
    print("\n=== Creating visualization ===")
    plot_raincloud(plot_data, output_file=args.output_plot)
    
    # ========================================================================
    # SAVE RESULTS
    # ========================================================================
    print("\n=== Saving results ===")
    
    # Summary statistics: mean, std, quartiles for each pair/category
    summary = plot_data.groupby(['pair', 'category'])['length_kb'].describe()
    summary.to_csv(args.output_summary, sep="\t")
    print(f"  Summary statistics saved to {args.output_summary}")
    
    # Detailed data: all individual length measurements
    plot_data.to_csv(args.output_detailed, sep="\t", index=False)
    print(f"  Detailed data saved to {args.output_detailed}")
    
    # ========================================================================
    # PRINT SUMMARY
    # ========================================================================
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print(summary)
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)

if __name__ == '__main__':
    main()