#!/usr/bin/env python3
"""
Analyze chromosomal distances for conflicting gene pairs between Cactus and ProtSyn.

Usage:
    python script.py --fn_table <path> --neighborhood_dict <path> \
                     --grass_scores <path> --output_tsv <path> --output_plot <path> \
                     --gene_gtf <file1> <file2> ... \
                     --cactus_genes <file1> <file2> ... \
                     --cactus_regions <file1> <file2> ...
"""

import sys
import pickle
import re
import pandas as pd
import matplotlib.pyplot as plt
from math import sqrt


def load_pickle(filepath):
    """Load a pickle file."""
    with open(filepath, 'rb') as f:
        return pickle.load(f)


def parse_gtf_attributes(attr_string):
    """
    Parse GTF attributes string and return a dictionary.
    Handles both quoted and unquoted values.
    """
    attributes = {}
    # Match key-value pairs: key "value" or key value
    pattern = r'(\w+)\s+"([^"]+)"|(\w+)\s+([^;]+)'
    
    for match in re.finditer(pattern, attr_string):
        if match.group(1):  # quoted value
            key = match.group(1)
            value = match.group(2)
        else:  # unquoted value
            key = match.group(3)
            value = match.group(4).strip()
        attributes[key] = value
    
    return attributes


def extract_gene_id(attributes_str):
    """
    Extract unique numeric gene ID from GTF attributes string.
    Searches for GeneID directly in the attributes string.
    Returns None if no numeric GeneID found.
    """
    # Search directly in the full attributes string for GeneID
    match = re.search(r'GeneID:(\d+)', attributes_str)
    if match:
        return match.group(1)
    return None


def load_gene_coordinates_from_gtf(gtf_files):
    """
    Load gene coordinates from GTF/GFF files.
    Returns dict: gene_id -> (chrom, start, end)
    Only processes lines with feature type 'gene'.
    Uses GeneID from db_xref as unique identifier.
    If a gene appears multiple times, keeps the last occurrence.
    """
    coords = {}
    
    for gtf_file in gtf_files:
        print(f"  Reading {gtf_file}...")
        line_count = 0
        gene_count = 0
        
        with open(gtf_file, 'r') as f:
            for line in f:
                # Skip comments
                if line.startswith('#'):
                    continue
                
                line_count += 1
                parts = line.strip().split('\t')
                
                # GTF/GFF should have 9 columns
                if len(parts) != 9:
                    continue
                
                chrom = parts[0]
                feature_type = parts[2]
                start = parts[3]
                end = parts[4]
                attributes_str = parts[8]
                
                # Only process gene features
                if feature_type != 'gene':
                    continue
                
                gene_count += 1
                
                # Extract gene ID directly from attributes string
                gene_id = extract_gene_id(attributes_str)
                
                if gene_id is None:
                    continue
                
                # Store coordinates (overwrites if duplicate)
                try:
                    coords[gene_id] = (chrom, int(start), int(end))
                except ValueError:
                    continue
        
        print(f"  Processed {line_count} lines, found {gene_count} gene features")
    
    return coords

def load_cactus_gene_to_region(gene_files):
    """
    Load multiple Cactus gene->region files.
    Returns dict: (species, gene_id) -> region_id
    """
    gene_to_region = {}
    for filepath in gene_files:
        df = pd.read_csv(filepath, sep='\t')
        for _, row in df.iterrows():
            species = row['Species']
            gene_id = str(row['Gene-ID'])  # Convert to string to match GTF
            region_id = row['Region_ID']
            gene_to_region[(species, gene_id)] = region_id
    return gene_to_region


def load_cactus_regions(region_files):
    """
    Load multiple Cactus synteny region files.
    Returns dict: region_id -> (query_start, query_end, target_start, target_end)
    """
    regions = {}
    for filepath in region_files:
        df = pd.read_csv(filepath, sep='\t')
        for _, row in df.iterrows():
            region_id = row['Region_ID']
            q_start = int(row['Query_Start'])
            q_end = int(row['Query_End'])
            t_start = int(row['Target_Start'])
            t_end = int(row.get('Target_End', t_start))  # adjust if needed
            regions[region_id] = (q_start, q_end, t_start, t_end)
    return regions


def compute_cactus_length(focus_species, focus_gene, gene_to_region, regions):
    """
    Compute geometric mean of Cactus region lengths.
    Returns None if gene not found in Cactus.
    """
    key = (focus_species, focus_gene)
    if key not in gene_to_region:
        return None
    
    region_id = gene_to_region[key]
    if region_id not in regions:
        return None
    
    q_start, q_end, t_start, t_end = regions[region_id]
    len_q = q_end - q_start + 1
    len_t = t_end - t_start + 1
    
    return sqrt(len_q * len_t)


def compute_neighborhood_max_distance(focus_gene, neighborhood, grass_scores, gene_coords):
    """
    Compute max distance from focus gene to furthest scored neighbor on same chromosome.
    Returns None if no valid neighbors found.
    """
    if focus_gene not in neighborhood:
        return None
    
    if focus_gene not in grass_scores:
        return None
    
    if focus_gene not in gene_coords:
        return None
    
    neighbors = neighborhood[focus_gene]
    scores_dict = dict(grass_scores[focus_gene])  # convert list of tuples to dict
    focus_chrom, focus_start, focus_end = gene_coords[focus_gene]
    focus_pos = (focus_start + focus_end) // 2
    
    distances = []
    for neighbor in neighbors:
        if neighbor not in scores_dict:
            continue
        if neighbor not in gene_coords:
            continue
        
        neighbor_chrom, neighbor_start, neighbor_end = gene_coords[neighbor]
        
        # Only consider neighbors on same chromosome
        if neighbor_chrom != focus_chrom:
            continue
        
        neighbor_pos = (neighbor_start + neighbor_end) // 2
        distance = abs(neighbor_pos - focus_pos)
        distances.append(distance)
    
    if not distances:
        return None
    
    return max(distances)


def parse_arguments(args):
    """Parse command line arguments."""
    params = {
        'fn_table': None,
        'gene_gtf': [],
        'neighborhood_dict': None,
        'grass_scores': None,
        'output_tsv': None,
        'output_plot': None,
        'cactus_genes': [],
        'cactus_regions': []
    }
    
    i = 0
    while i < len(args):
        arg = args[i]
        
        if arg == '--fn_table':
            params['fn_table'] = args[i + 1]
            i += 2
        elif arg == '--gene_gtf':
            i += 1
            while i < len(args) and not args[i].startswith('--'):
                params['gene_gtf'].append(args[i])
                i += 1
        elif arg == '--neighborhood_dict':
            params['neighborhood_dict'] = args[i + 1]
            i += 2
        elif arg == '--grass_scores':
            params['grass_scores'] = args[i + 1]
            i += 2
        elif arg == '--output_tsv':
            params['output_tsv'] = args[i + 1]
            i += 2
        elif arg == '--output_plot':
            params['output_plot'] = args[i + 1]
            i += 2
        elif arg == '--cactus_genes':
            i += 1
            while i < len(args) and not args[i].startswith('--'):
                params['cactus_genes'].append(args[i])
                i += 1
        elif arg == '--cactus_regions':
            i += 1
            while i < len(args) and not args[i].startswith('--'):
                params['cactus_regions'].append(args[i])
                i += 1
        else:
            print(f"Unknown argument: {arg}")
            print(__doc__)
            sys.exit(1)
    
    # Validate required parameters
    required = ['fn_table', 'neighborhood_dict', 'grass_scores', 
                'output_tsv', 'output_plot']
    missing = [p for p in required if params[p] is None]
    
    if missing:
        print(f"Error: Missing required parameters: {', '.join(missing)}")
        print(__doc__)
        sys.exit(1)
    
    if not params['gene_gtf']:
        print("Error: No GTF files provided (--gene_gtf)")
        sys.exit(1)
    
    if not params['cactus_genes']:
        print("Error: No Cactus gene files provided (--cactus_genes)")
        sys.exit(1)
    
    if not params['cactus_regions']:
        print("Error: No Cactus region files provided (--cactus_regions)")
        sys.exit(1)
    
    return params


def main():
    # Parse arguments
    params = parse_arguments(sys.argv[1:])
    
    # ==========================
    # LOAD DATA
    # ==========================
    print("Loading gene-level FN table...")
    fn_df = pd.read_csv(params['fn_table'], sep='\t')
    print(f"  Loaded {len(fn_df)} rows")
    
    print(f"\nLoading gene coordinates from {len(params['gene_gtf'])} GTF files...")
    gene_coords = load_gene_coordinates_from_gtf(params['gene_gtf'])
    print(f"Loaded {len(gene_coords)} unique gene coordinates")
    print(f"Sample gene IDs from GTF: {list(gene_coords.keys())[:5]}")
    
    print("\nLoading neighborhood and GRASS scores...")
    neighborhood = load_pickle(params['neighborhood_dict'])
    grass_scores = load_pickle(params['grass_scores'])
    print(f"Neighborhood has {len(neighborhood)} genes")
    print(f"GRASS scores has {len(grass_scores)} genes")
    print(f"Sample gene IDs from neighborhood: {list(neighborhood.keys())[:5]}")
    
    print(f"\nLoading {len(params['cactus_genes'])} Cactus gene->region mappings...")
    gene_to_region = load_cactus_gene_to_region(params['cactus_genes'])
    print(f"Loaded {len(gene_to_region)} gene-region mappings")
    print(f"Sample keys from Cactus genes: {list(gene_to_region.keys())[:5]}")
    
    print(f"\nLoading {len(params['cactus_regions'])} Cactus synteny regions...")
    regions = load_cactus_regions(params['cactus_regions'])
    print(f"Loaded {len(regions)} synteny regions")
    
    # ==========================
    # FILTER CONTRADICTORY ROWS
    # ==========================
    print("\nFiltering contradictory pairs...")
    contradictory = fn_df[
        (fn_df['FN_Class'] == 'contradictory') &
        (fn_df['ProtSyn_Ortholog_Gene'].notna()) &
        (fn_df['ProtSyn_Ortholog_Gene'] != '')
    ].copy()
    
    print(f"Found {len(contradictory)} contradictory pairs")
    if len(contradictory) > 0:
        print(f"Sample Focus_Species: {contradictory['Focus_Species'].iloc[0]}")
        print(f"Sample Focus_Gene: {contradictory['Focus_Gene'].iloc[0]}")
    
    # ==========================
    # COMPUTE DISTANCES
    # ==========================
    print("\nComputing distances...")
    results = []
    missing_cactus = 0
    missing_neighborhood = 0

    # Diagnostic counters
    not_in_neighborhood = 0
    not_in_grass = 0
    not_in_coords = 0
    no_valid_neighbors = 0
    
    for idx, row in contradictory.iterrows():
        fn_id = row['FN_ID']
        pair = row['Pair']
        focus_species = row['Focus_Species']
        focus_gene = str(row['Focus_Gene'])  # Convert to string
        protsyn_ortholog = row['ProtSyn_Ortholog_Gene']
        
        # Compute Cactus length
        cactus_len = compute_cactus_length(focus_species, focus_gene, gene_to_region, regions)
        
        # Detailed neighborhood diagnostics
        if focus_gene not in neighborhood:
            not_in_neighborhood += 1
            neighborhood_max_dist = None
        elif focus_gene not in grass_scores:
            not_in_grass += 1
            neighborhood_max_dist = None
        elif focus_gene not in gene_coords:
            not_in_coords += 1
            neighborhood_max_dist = None
        else:
            neighborhood_max_dist = compute_neighborhood_max_distance(
                focus_gene, neighborhood, grass_scores, gene_coords
            )
            if neighborhood_max_dist is None:
                no_valid_neighbors += 1


        if cactus_len is None:
            missing_cactus += 1
        if neighborhood_max_dist is None:
            missing_neighborhood += 1
        
        # Only keep if both values are available
        if cactus_len is not None and neighborhood_max_dist is not None:
            region_id = gene_to_region.get((focus_species, focus_gene), '')
            results.append({
                'FN_ID': fn_id,
                'Pair': pair,
                'Focus_Species': focus_species,
                'Focus_Gene': focus_gene,
                'ProtSyn_Ortholog_Gene': protsyn_ortholog,
                'Region_ID': region_id,
                'cactus_len': cactus_len,
                'neighborhood_max_dist': neighborhood_max_dist
            })
    
    print(f"\nComputed distances for {len(results)} pairs")
    print(f"Missing Cactus data: {missing_cactus} ({missing_cactus/len(contradictory)*100:.1f}%)")
    print(f"Missing neighborhood data: {missing_neighborhood} ({missing_neighborhood/len(contradictory)*100:.1f}%)")
    print(f"\nNeighborhood missing breakdown:")
    print(f"  Not in neighborhood.pkl: {not_in_neighborhood}")
    print(f"  Not in grass_scores.pkl: {not_in_grass}")
    print(f"  Not in GTF coords: {not_in_coords}")
    print(f"  No valid neighbors found: {no_valid_neighbors}")
    
    # ==========================
    # SAVE RESULTS
    # ==========================
    if len(results) > 0:
        results_df = pd.DataFrame(results)
        results_df.to_csv(params['output_tsv'], sep='\t', index=False)
        print(f"\nSaved results to {params['output_tsv']}")
    else:
        print("\nWARNING: No results to save!")
        return
    
    # ==========================
    # PLOT
    # ==========================
    if len(results) > 0:
        plt.figure(figsize=(10, 8))
        plt.scatter(results_df['cactus_len'], results_df['neighborhood_max_dist'], 
                   alpha=0.6, s=50)
        plt.xlabel('Cactus Region Length (geometric mean)', fontsize=12)
        plt.ylabel('Neighborhood Max Distance', fontsize=12)
        plt.title('Chromosomal Distances for Conflicting Pairs', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(params['output_plot'], dpi=300)
        print(f"Saved plot to {params['output_plot']}")
        plt.close()
    else:
        print("No data to plot")


if __name__ == '__main__':
    main()