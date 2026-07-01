import pandas as pd
import sys

def load_data(fn_file, protsyn_file):
    """Load the false negative and prot-syn ortholog tables."""
    fn_df = pd.read_csv(fn_file, sep='\t')
    protsyn_df = pd.read_csv(protsyn_file, sep='\t')
    
    # Normalize species names to handle potential formatting differences
    # Remove underscores and convert to lowercase for matching
    protsyn_df['Species1_norm'] = protsyn_df['Species1'].str.replace('_', ' ').str.lower()
    protsyn_df['Species2_norm'] = protsyn_df['Species2'].str.replace('_', ' ').str.lower()
    
    return fn_df, protsyn_df

def find_protsyn_ortholog(query_gene, query_species, target_species, protsyn_df):
    """
    Search for the query gene in prot-syn table for the specific species pair.
    
    Returns:
        tuple: (found, protein, gene) where found is boolean, 
               protein and gene are the prot-syn assignments (None if not found)
    """
    # Normalize species names for matching
    query_species_norm = query_species.lower()
    target_species_norm = target_species.lower()
    
    # Search where query gene is in Gene1 position
    match1 = protsyn_df[
        (protsyn_df['Gene1'] == query_gene) & 
        (protsyn_df['Species1_norm'] == query_species_norm) &
        (protsyn_df['Species2_norm'] == target_species_norm)
    ]
    
    if not match1.empty:
        row = match1.iloc[0]
        return True, row['Protein2'], row['Gene2']
    
    # Search where query gene is in Gene2 position
    match2 = protsyn_df[
        (protsyn_df['Gene2'] == query_gene) & 
        (protsyn_df['Species2_norm'] == query_species_norm) &
        (protsyn_df['Species1_norm'] == target_species_norm)
    ]
    
    if not match2.empty:
        row = match2.iloc[0]
        return True, row['Protein1'], row['Gene1']
    
    # No match found
    return False, None, None

def classify_false_negatives(fn_df, protsyn_df):
    """
    Classify each false negative as PS-no-results or PS-contradictory.
    """
    results = []
    
    for idx, row in fn_df.iterrows():
        query_gene = row['Query_Gene']
        query_protein = row['Query_Protein']
        query_species = row['Query_Species']
        target_species = row['Target_Species']
        cactus_target_protein = row['Target_Protein']
        cactus_target_gene = row['Target_Gene']
        
        # Search for prot-syn ortholog
        found, protsyn_protein, protsyn_gene = find_protsyn_ortholog(
            query_gene, query_species, target_species, protsyn_df
        )
        
        if not found:
            label = 'PS-no-results'
            protsyn_protein = None
            protsyn_gene = None
        else:
            # Check if prot-syn found a different gene than Cactus
            if protsyn_gene != cactus_target_gene:
                label = 'PS-contradictory'
            else:
                # This shouldn't happen (would mean it's not a false negative)
                label = 'PS-agreement'  # Edge case
        
        results.append({
            'Query_Protein': query_protein,
            'Query_Gene': query_gene,
            'Query_Species': query_species,
            'Target_Species': target_species,
            'Label': label,
            'Cactus_Target_Protein': cactus_target_protein,
            'Cactus_Target_Gene': cactus_target_gene,
            'ProtSyn_Target_Protein': protsyn_protein,
            'ProtSyn_Target_Gene': protsyn_gene
        })
        
        # Progress indicator
        if (idx + 1) % 100 == 0:
            print(f"Processed {idx + 1}/{len(fn_df)} rows...", file=sys.stderr)
    
    return pd.DataFrame(results)

def main():
    # File paths - modify these as needed
    fn_file = 'results/evaluation/mammalia/vs_prot_syn/nscore_approach/false_negatives_can_mus.tsv'
    protsyn_file = 'material/sex_experiment/one_to_one_nscore.tsv'
    output_file = 'results/evaluation/mammalia/vs_prot_syn/nscore_approach/disagreement_investigation/false_negatives/classified_false_negatives_can_mus.tsv'
    
    print("Loading data...", file=sys.stderr)
    fn_df, protsyn_df = load_data(fn_file, protsyn_file)
    
    print(f"Loaded {len(fn_df)} false negatives", file=sys.stderr)
    print(f"Loaded {len(protsyn_df)} prot-syn ortholog pairs", file=sys.stderr)
    
    print("\nClassifying false negatives...", file=sys.stderr)
    classified_df = classify_false_negatives(fn_df, protsyn_df)
    
    # Save results
    classified_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nResults saved to {output_file}", file=sys.stderr)
    
    # Print summary statistics
    print("\n=== Classification Summary ===", file=sys.stderr)
    print(classified_df['Label'].value_counts(), file=sys.stderr)
    
    # Show a few examples
    print("\n=== First 5 rows ===", file=sys.stderr)
    print(classified_df.head().to_string(), file=sys.stderr)

if __name__ == "__main__":
    main()
