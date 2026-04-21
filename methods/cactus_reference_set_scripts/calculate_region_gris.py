import pandas as pd
import numpy as np
import os

"""
Script: Region GRIS Score Calculation

Purpose: Calculates GRIS (Genomic Region Identity Score) for synteny blocks by matching gene pairs 
between species using protein similarity scores from BLAST/DIAMOND results.

Input:
- All-vs-all BLAST/DIAMOND results with geometric mean of protein similarity scores (.txt)
- Protein-to-gene mapping file (.tsv)
- Tables with syntenic genes for each species pair (*_synteny_tables_genes_fix.tsv)

Output:
- Tables with region GRIS scores for each species pair (*_region_gris_{method}.tsv)
  Columns: Region_ID, Query_Genes, Target_Genes, Num_Pairs, reg-GRIS

Prints:
- Number of valid DIAMOND rows loaded
- Example DIAMOND rows
- Number of protein mappings loaded
- Total gene pairs mapped
- Example mapped gene pairs
- For each species pair:
  - Number of regions with/without gene pairs
  - Unique regions with at least one gene
  - Regions with GRIS data
  - Sum verification (regions with + without pairs)
  - Output file path
"""

# From the gene tables build the species pairs based on their GRIS score (aka geometrical mean)
# Calculate region GRIS score for each region


# Input files: 
# - all-vs-all BLAST/DIAMOND results with geometrical mean of similarity score of proteins (.txt)
# - Protein-gene map (.tsv)
# - Tables with syntenic genes for each species pair (.tsv)
# Output:
# - Tables with region GRIS scores for each species pair (.tsv)
# - Statistical information is printed as stdout


# Full names of the species (for search for genes in syntenic gene tables)
full_names = {
    "can": "Canis lupus familiaris",
    "mac": "Macaca fascicularis",
    "rat": "Rattus norvegicus",
    "mus": "Mus musculus", 
}


def extract_from_blast(blast_file):
    """
    Extract protein pairs and their similarity scores from BLAST/DIAMOND results.
    
    Input:
    - blast_file: str, path to BLAST/DIAMOND output file (CSV format)
    
    Output:
    - DataFrame with columns: qseqid (str), sseqid (str), geometric_mean (float)
    
    Prints:
    - Number of valid DIAMOND rows after cleanup
    - First 3 example rows
    """
    colnames = ["qseqid","sseqid","qstart","qend","qlen","sstart","send","slen","pident","evalue","bitscore","overlap","geometric_mean"]
    #blast_df = pd.read_csv(blast_file, sep=",", usecols=["qseqid", "sseqid", "geometric_mean"], engine="python")
    blast_df = pd.read_csv(blast_file, sep=",", names=colnames, engine="python")
    blast_df = blast_df[["qseqid","sseqid","geometric_mean"]]
    blast_df.columns = blast_df.columns.str.strip()  # remove stray spaces
    blast_df = blast_df.dropna(subset=["qseqid", "sseqid", "geometric_mean"])
    blast_df = blast_df[blast_df["qseqid"] != blast_df["sseqid"]]
    blast_df["geometric_mean"] = pd.to_numeric(blast_df["geometric_mean"], errors="coerce")
    blast_df = blast_df.dropna(subset=["geometric_mean"])  

    print(f"Loaded {len(blast_df):,} valid DIAMOND rows after cleanup")
    print("Example:", blast_df.head(3))

    return blast_df[["qseqid","sseqid","geometric_mean"]]


def extract_from_map(map_file):
    """
    Extract protein-to-gene ID mapping.
    
    Input:
    - map_file: str, path to TSV file with columns GeneID and ProteinID (order may vary)
    
    Output:
    - dict, mapping ProteinID (str) -> GeneID (str)
    
    Prints: Nothing
    """
    map_df = pd.read_csv(map_file, sep="\t")

    # Address different column order
    if map_df.columns[0].lower().startswith("protein"):
        map_df = map_df[[map_df.columns[1], map_df.columns[0]]]  # swap columns


    # For mammals and plants
    map_df = map_df[["GeneID", "ProteinID"]]
    map_information = map_df[["GeneID", "ProteinID"]]

    # For drosophila
    #map_df = map_df[["ProteinID", "GeneID"]]
    #map_information = map_df[["ProteinID", "GeneID"]]

    # Create dictionary for faster lookup: ProteinID -> GeneID
    protein_to_gene_map = dict(zip(map_information['ProteinID'], map_information['GeneID']))

    return protein_to_gene_map


def map_proteins_to_genes(geomean_df, protein_to_gene):
    """
    Map protein IDs to gene IDs using the protein-to-gene dictionary.
    
    Input:
    - geomean_df: DataFrame with columns qseqid, sseqid, geometric_mean (from extract_from_blast)
    - protein_to_gene: dict, ProteinID -> GeneID mapping (from extract_from_map)
    
    Output:
    - DataFrame with columns: Query (str), Target (str), geometric_mean (float)
    
    Prints:
    - Total number of gene pairs mapped
    """
    # Map both qseqid and sseqid to gene IDs
    geomean_df['Query'] = geomean_df['qseqid'].map(protein_to_gene)
    geomean_df['Target'] = geomean_df['sseqid'].map(protein_to_gene)
    
    # Keep only rows where both proteins mapped to genes
    gene_pairs = geomean_df.dropna(subset=['Query', 'Target'])
    
    # Keep only relevant columns
    gene_pairs = gene_pairs[['Query', 'Target', 'geometric_mean']]
    
    print(f"Total gene pairs mapped: {len(gene_pairs)}\n")
    
    return gene_pairs


def calculate_region_gris(genes_table_file, gene_pairs_df, species1, species2, score_method='mean'):
    """
    Calculate GRIS scores for all synteny regions.
    
    Input:
    - genes_table_file: str, path to synteny genes table TSV with columns Region_ID, Species, Gene-ID
    - gene_pairs_df: DataFrame with columns Query, Target, geometric_mean (from map_proteins_to_genes)
    - species1: str, first species code (e.g., "can")
    - species2: str, second species code (e.g., "mac")
    - score_method: str, aggregation method - 'sum', 'mean', 'normalized', or 'geomean' (default 'mean')
    
    Output:
    - DataFrame with columns: Region_ID (int), Query_Genes (int), Target_Genes (int), Num_Pairs (int), reg_GRIS (float)
    
    Prints:
    - Number of regions with gene pairs
    - Number of regions with no pairs
    """
    # Assign full species names to their code
    species1_full = full_names[species1]
    species2_full = full_names[species2]

    genes_df = pd.read_csv(genes_table_file, sep="\t", dtype=str) #
          
    # Normalize columns to avoid whitespace-caused misses
    for col in ["Region_ID", "Species", "Gene-ID"]:
        if col in genes_df.columns:
            genes_df[col] = genes_df[col].astype(str).str.strip()
        else:
            raise KeyError(f"Missing column {col} in {genes_table_file}")

    # Keep only these two species
    genes_df = genes_df[genes_df["Species"].isin([species1_full, species2_full])].copy()

    # Split into per-species tables of (Region_ID, Gene-ID)
    sp1 = genes_df.loc[genes_df["Species"] == species1_full, ["Region_ID", "Gene-ID"]].rename(columns={"Gene-ID": "Query"})
    sp2 = genes_df.loc[genes_df["Species"] == species2_full, ["Region_ID", "Gene-ID"]].rename(columns={"Gene-ID": "Target"})

    # normalize data types for merging
    gene_pairs_df["Query"]  = gene_pairs_df["Query"].astype(str)
    gene_pairs_df["Target"] = gene_pairs_df["Target"].astype(str)

    # Inner-join gene_pairs onto each side to find all Query/Target that exist in-region
    left  = gene_pairs_df.merge(sp1, on="Query",  how="inner")   # adds Region_ID for the query side -> Region_ID_x
    both  = left.merge(sp2, on="Target", how="inner", suffixes=("_q", "_t"))


    # Keep only pairs that fall in the SAME region (same numeric Region_ID)
    # Handle numeric/string safely
    both["Region_ID_q"] = both["Region_ID_q"].astype(str).str.strip()
    both["Region_ID_t"] = both["Region_ID_t"].astype(str).str.strip()
    in_same_region = both[both["Region_ID_q"] == both["Region_ID_t"]].copy()
    in_same_region.rename(columns={"Region_ID_q": "Region_ID"}, inplace=True)

    # If the same gene pair appears multiple times in a region, collapse it
    in_same_region = in_same_region.drop_duplicates(subset=["Region_ID", "Query", "Target"])

    # number of pairs per region + summary score
    if score_method == "sum":
        agg = in_same_region.groupby("Region_ID")["geometric_mean"].agg(Num_Pairs="count", reg_GRIS="sum").reset_index()
    elif score_method == "mean":
        agg = in_same_region.groupby("Region_ID")["geometric_mean"].agg(Num_Pairs="count", reg_GRIS="mean").reset_index()
    elif score_method == "normalized":
        # need |genes_q| and |genes_t| per region for denominator -> compute below
        base = in_same_region.groupby("Region_ID")["geometric_mean"].agg(Num_Pairs="count", sum_score="sum").reset_index()
    elif score_method == "geomean":
        # geometric mean per region
        grp = in_same_region.groupby("Region_ID")["geometric_mean"]
        agg = grp.apply(lambda s: pd.Series({"Num_Pairs": s.size, "reg_GRIS": float(np.exp(np.log(s[s > 0]).mean())) if (s > 0).any() else 0.0})).reset_index()
    else:
        # default to sum
        agg = in_same_region.groupby("Region_ID")["geometric_mean"].agg(Num_Pairs="count", reg_GRIS="sum").reset_index()

    # per-region unique gene counts on each side (for reporting)
    q_counts = sp1.groupby("Region_ID")["Query"].nunique().rename("Query_Genes")
    t_counts = sp2.groupby("Region_ID")["Target"].nunique().rename("Target_Genes")
    size_tbl = pd.concat([q_counts, t_counts], axis=1).reset_index().fillna(0)
    size_tbl["Query_Genes"]  = size_tbl["Query_Genes"].astype(int)
    size_tbl["Target_Genes"] = size_tbl["Target_Genes"].astype(int)

    # Combine
    results_df = size_tbl.merge(agg, on="Region_ID", how="left")

    # For regions with no pairs, fill Num_Pairs/reg_GRIS
    results_df["Num_Pairs"] = results_df["Num_Pairs"].fillna(0).astype(int)
    if score_method in ("sum", "mean", "geomean"):
        results_df["reg_GRIS"] = results_df["reg_GRIS"].fillna(0.0)
    elif score_method == "normalized":
        # compute normalized = sum_score / (|Q| * |T|)
        results_df = results_df.merge(base, on="Region_ID", how="left")
        results_df["sum_score"] = results_df["sum_score"].fillna(0.0)
        denom = (results_df["Query_Genes"] * results_df["Target_Genes"]).replace(0, np.nan)
        results_df["reg_GRIS"] = (results_df["sum_score"] / denom).fillna(0.0)
        results_df.drop(columns=["sum_score"], inplace=True)

    # Make Region_ID numeric-sortable but keep as string if needed
    results_df["Region_ID"] = pd.to_numeric(results_df["Region_ID"], errors="coerce").fillna(-1).astype(int)
    results_df = results_df.sort_values("Region_ID").reset_index(drop=True)

    # Quick report
    with_pairs = int((results_df["Num_Pairs"] > 0).sum())
    no_pairs   = int((results_df["Num_Pairs"] == 0).sum())
    print(f"\nRegions with gene pairs: {with_pairs}")
    print(f"Regions with no pairs:   {no_pairs}")

    return results_df


def main():
    """
    Main function to calculate GRIS scores for all species pairs.
    
    Prints:
    - Step-by-step progress messages
    - Statistics from each processing step
    - Output file paths for each species pair
    - Completion message
    """
    # File paths

    # Uncomment path needed 

    # For mammals
    blast_file = "/storage/EasyVectorOmics/synteny_algorithm/material/mammalian/geometric_mean.txt"
    map_file = "/storage/EasyVectorOmics/synteny_algorithm/material/mammalian/gene_to_protein_map.tsv"

    # For plants
    #blast_file = "/storage/EasyVectorOmics/synteny_algorithm/material/cardamine/geometric_mean.txt"
    #map_file = "/storage/EasyVectorOmics/synteny_algorithm/material/cardamine/gene_to_protein_map.tsv"

    # For drosophila
    #blast_file = "/storage/EasyVectorOmics/synteny_algorithm/material/drosophila/geometric_mean.txt"
    #map_file = "/storage/EasyVectorOmics/synteny_algorithm/material/drosophila/gene_to_protein_map.tsv"

    # Score method
    score_method = 'mean'
    
    # Species pairs to process
    # Uncomment for needed input

    # For mammals
    species_pairs = [
        ('can', 'mac', 'results/genes_in_blocks/mammalia/region_gene_tables/can_mac_synteny_tables_genes_fix.tsv'),
        ('mac', 'can', 'results/genes_in_blocks/mammalia/region_gene_tables/mac_can_synteny_tables_genes_fix.tsv'),
        ('can', 'rat', 'results/genes_in_blocks/mammalia/region_gene_tables/can_rat_synteny_tables_genes_fix.tsv'),
        ('rat', 'can', 'results/genes_in_blocks/mammalia/region_gene_tables/rat_can_synteny_tables_genes_fix.tsv'),
        ('can', 'mus', 'results/genes_in_blocks/mammalia/region_gene_tables/can_mus_synteny_tables_genes_fix.tsv'),
        ('mus', 'can', 'results/genes_in_blocks/mammalia/region_gene_tables/mus_can_synteny_tables_genes_fix.tsv'),
        ('mac', 'rat', 'results/genes_in_blocks/mammalia/region_gene_tables/mac_rat_synteny_tables_genes_fix.tsv'),
        ('rat', 'mac', 'results/genes_in_blocks/mammalia/region_gene_tables/rat_mac_synteny_tables_genes_fix.tsv'),
        ('mac', 'mus', 'results/genes_in_blocks/mammalia/region_gene_tables/mac_mus_synteny_tables_genes_fix.tsv'),
        ('mus', 'mac', 'results/genes_in_blocks/mammalia/region_gene_tables/mus_mac_synteny_tables_genes_fix.tsv'),
        ('rat', 'mus', 'results/genes_in_blocks/mammalia/region_gene_tables/rat_mus_synteny_tables_genes_fix.tsv'),
        ('mus', 'rat', 'results/genes_in_blocks/mammalia/region_gene_tables/mus_rat_synteny_tables_genes_fix.tsv')
    ]

    # For plants
    '''
    species_pairs = [
        ('ar', 'card', 'results/genes_in_blocks/brassicaceae/region_gene_tables/ar_card_synteny_tables_genes.tsv'),
        ('card', 'ar', 'results/genes_in_blocks/brassicaceae/region_gene_tables/card_ar_synteny_tables_genes.tsv')
    ]
    '''
    
    # For drosophila
    '''
    species_pairs = [
        ('dmel', 'dsec', 'results/genes_in_blocks/drosophila/region_gene_tables/dmel_dsec_synteny_tables_genes.tsv'),
        ('dsec', 'dmel', 'results/genes_in_blocks/drosophila/region_gene_tables/dsec_dmel_synteny_tables_genes.tsv'),
        ('dmel', 'dsim', 'results/genes_in_blocks/drosophila/region_gene_tables/dmel_dsim_synteny_tables_genes.tsv'),
        ('dsim', 'dmel', 'results/genes_in_blocks/drosophila/region_gene_tables/dsim_dmel_synteny_tables_genes.tsv'),
        ('dsec', 'dsim', 'results/genes_in_blocks/drosophila/region_gene_tables/dsec_dsim_synteny_tables_genes.tsv'),
        ('dsim', 'dsec', 'results/genes_in_blocks/drosophila/region_gene_tables/dsim_dsec_synteny_tables_genes.tsv')
    ]
    '''
    
    
    # main
    print("Step 1: Loading DIAMOND results...")
    geomean_df = extract_from_blast(blast_file)
    print(f"Loaded {len(geomean_df)} protein pairs")
    
    print("\nStep 2: Loading protein-to-gene mapping...")
    protein_to_gene = extract_from_map(map_file)
    print(f"Loaded {len(protein_to_gene)} protein mappings")
    
    print("\nStep 3: Mapping proteins to genes...")
    gene_pairs_df = map_proteins_to_genes(geomean_df, protein_to_gene)

    print("First five entries from mapped dataframe:")
    print(gene_pairs_df.head(5))
    
    print("\nStep 4: Calculating region GRIS scores...\n")
    
    # Create output directory
    output_dir = "results/region_gris"
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each species pair
    for sp1, sp2, genes_table in species_pairs:
        results_df = calculate_region_gris(
            genes_table, 
            gene_pairs_df, 
            sp1, 
            sp2, 
            score_method
        )
        
        # Save results
        output_file = f"{output_dir}/{sp1}_{sp2}_region_gris_{score_method}_fix.tsv"

        results_df.to_csv(output_file, sep="\t", index=False)
        print(f"Saved: {output_file}")

        genes_df = pd.read_csv(genes_table, sep="\t", dtype=str)
        print("Unique regions with at least one gene:", genes_df["Region_ID"].nunique())

        gris_df = pd.read_csv(output_file, sep="\t")
        print("Regions with GRIS data:", gris_df["Region_ID"].nunique())
        print("Sum of (regions with pairs + regions without pairs):", (gris_df["Num_Pairs"] > 0).sum() + (gris_df["Num_Pairs"] == 0).sum())

    
    print("All calculations complete!")
    
    

if __name__ == "__main__":
    main()