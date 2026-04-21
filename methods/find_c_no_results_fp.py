import pandas as pd

# Among FP pairs finds those that have no result in Cactus, and saves them for further investigation

'''
Check Gene1: if found anywhere in cactus → ignore the whole FP pair
If Gene1 not found → check Gene2: if found anywhere → ignore
Only if BOTH genes are not found anywhere → C_no_result
'''


def load_false_positives(fp_file):
    """Load false positives table."""
    return pd.read_csv(fp_file, sep='\t')

def load_cactus_orthologs(cactus_file):
    """Load cactus ortholog table."""
    return pd.read_csv(cactus_file, sep='\t')

def check_gene_in_cactus(gene, query_genes_set, target_genes_set):
    """
    Check if gene is in cactus table (query or target column).
    
    Args:
        gene: Gene ID to check
        query_genes_set: Set of genes from Query_Gene column
        target_genes_set: Set of genes from Target_Gene column
    
    Returns:
        'query' if found in Query_Gene column
        'target' if found in Target_Gene column
        'none' if not found anywhere
    """
    if gene in query_genes_set:
        return 'query'
    elif gene in target_genes_set:
        return 'target'
    else:
        return 'none'

def should_be_no_result(gene, query_genes_set, target_genes_set):
    """
    Apply the decision logic for one gene.
    
    Logic:
    - If gene in Query column -> False (cactus found it, ignore)
    - If gene in Target column -> False (contradictory, ignore)
    - If gene nowhere -> True (potential C_no_result)
    
    Returns:
        True if this gene indicates C_no_result, False otherwise
    """
    location = check_gene_in_cactus(gene, query_genes_set, target_genes_set)
    return location == 'none'

def find_no_results_bidirectional(fp_df, cactus_files_dict):
    """
    Find genes from FP table that have no results in cactus (bidirectional check).
    
    For each FP row with Gene1 and Gene2:
    1. Check Gene1 in cactus (Query and Target columns)
       - If in Query -> ignore entire row
       - If in Target -> ignore entire row
       - If nowhere -> potential C_no_result, continue to step 2
    2. Check Gene2 in cactus (Query and Target columns)
       - If in Query -> ignore entire row
       - If in Target -> ignore entire row
       - If nowhere -> C_no_result (both genes not found)
    
    Args:
        fp_df: False positives dataframe with all pairs
        cactus_files_dict: Dict mapping species pair to cactus dataframe
    
    Returns:
        Dataframe with C_no_result cases
    """
    no_results_list = []
    
    # Pre-compute gene sets for each species pair for faster lookup
    cactus_gene_sets = {}
    for species_pair, cactus_df in cactus_files_dict.items():
        query_genes = set(cactus_df['Query_Gene'].dropna())
        target_genes = set(cactus_df['Target_Gene'].dropna())
        cactus_gene_sets[species_pair] = (query_genes, target_genes)
        print(f"Prepared gene sets for {species_pair}: {len(query_genes)} query, {len(target_genes)} target")
    
    for idx, row in fp_df.iterrows():
        gene1 = row['Gene1']
        gene2 = row['Gene2']
        species1 = row['Species1']
        species2 = row['Species2']
        
        # Create species pair key
        species_pair = get_species_pair_key(species1, species2)
        
        # Get corresponding cactus gene sets
        if species_pair not in cactus_gene_sets:
            print(f"Warning: No cactus file for species pair {species_pair}")
            continue
        
        query_genes_set, target_genes_set = cactus_gene_sets[species_pair]
        
        # Check both genes
        gene1_is_no_result = should_be_no_result(gene1, query_genes_set, target_genes_set)
        gene2_is_no_result = should_be_no_result(gene2, query_genes_set, target_genes_set)
        
        # Only add to C_no_result if BOTH genes are not found anywhere
        if gene1_is_no_result and gene2_is_no_result:
            no_results_list.append({
                'Protein1': row['Protein1'],
                'Gene1': row['Gene1'],
                'Species1': row['Species1'],
                'Protein2_prot_syn': row['Protein2'],
                'Gene2_prot_syn': row['Gene2'],
                'Species2_prot_syn': row['Species2'],
                'Pair': row['Pair'],
                'Swapped': row['Swapped']
            })
    
    return pd.DataFrame(no_results_list)

def get_species_pair_key(species1, species2):
    """
    Normalize species names to create consistent key.
    """
    s1 = species1.lower().replace(' ', '_')
    s2 = species2.lower().replace(' ', '_')
    return f"{s1}-{s2}"

def load_cactus_files(cactus_file_paths):
    """
    Load multiple cactus files into a dictionary.
    """
    cactus_dict = {}
    for species_pair, file_path in cactus_file_paths.items():
        cactus_dict[species_pair] = load_cactus_orthologs(file_path)
        print(f"Loaded cactus file for {species_pair}: {len(cactus_dict[species_pair])} rows")
    
    return cactus_dict

def process_all_false_positives(fp_file, cactus_file_paths, output_file):
    """
    Process all false positives against multiple cactus files.
    """
    # Load data
    fp_df = load_false_positives(fp_file)
    print(f"Loaded {len(fp_df)} false positive pairs")
    
    cactus_dict = load_cactus_files(cactus_file_paths)
    
    # Find no results
    no_results = find_no_results_bidirectional(fp_df, cactus_dict)
    
    # Save results
    no_results.to_csv(output_file, sep='\t', index=False)
    print(f"\nFound {len(no_results)} C_no_result cases")
    print(f"Results saved to {output_file}")
    
    return no_results

# Example usage:
if __name__ == "__main__":
    fp_file = "results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/all_false_positives_sum.tsv"
    
    cactus_files = {
        "canis_lupus_familiaris-macaca_fascicularis": "results/reciprocal_pairs/can_mac_one2one_fix.tsv",
        "canis_lupus_familiaris-rattus_norvegicus": "results/reciprocal_pairs/can_rat_one2one_fix.tsv",
        "canis_lupus_familiaris-mus_musculus": "results/reciprocal_pairs/can_mus_one2one_fix.tsv",
        "macaca_fascicularis-rattus_norvegicus": "results/reciprocal_pairs/mac_rat_one2one_fix.tsv",
        "macaca_fascicularis-mus_musculus": "results/reciprocal_pairs/mac_mus_one2one_fix.tsv",
        "rattus_norvegicus-mus_musculus": "results/reciprocal_pairs/rat_mus_one2one_fix.tsv",
    }
    
    output_file = "results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_positives/all_C_no_results.tsv"
    
    process_all_false_positives(fp_file, cactus_files, output_file)