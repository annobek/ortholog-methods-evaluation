import pandas as pd
from pathlib import Path
import pickle
import re
import time

# -------------------------------------------------------------
# 1. Extract gene IDs from top 10 cases
# -------------------------------------------------------------

out_path_better = "results/evaluation/mammalia/vs_prot_syn/cases_cactus_picked_better_n_score.tsv"
neighborhood_dict_pkl = "/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/neighborhoods.pkl"
map_protein_genes_path = "/storage/EasyVectorOmics/synteny_algorithm/material/mammalian/gene_to_protein_map.tsv"


# Load all files and dicts
with open(neighborhood_dict_pkl, "rb") as f:
    nbh_dict = pickle.load(f)          # expected: {GeneID(str): [neighbor_gene1, neighbor_gene2, ...]}


map_df = pd.read_csv(map_protein_genes_path, sep="\t")



# Load the cases where Cactus picked better
cases_better = pd.read_csv(out_path_better, sep="\t")
top_cases = cases_better.head(10)

# Collect all genes involved (A_gene, B_gene, X_gene)
test_genes = set()
for _, row in top_cases.iterrows():
    test_genes.add(str(row["A_gene"]))
    test_genes.add(str(row["B_gene"]))
    test_genes.add(str(row["X_gene"]))

print(f"Genes from top 10 cases: {len(test_genes)}")
print(f"Gene IDs: {sorted(test_genes)}")

# Also add their neighborhoods for context
print("\nAdding neighborhood genes...")
all_genes_with_neighborhoods = test_genes.copy()
for gene in test_genes:
    neighbors = nbh_dict.get(gene, [])
    all_genes_with_neighborhoods.update(str(n) for n in neighbors)

print(f"Total genes including neighborhoods: {len(all_genes_with_neighborhoods)}")

# Save gene list
gene_list_file = "results/evaluation/mammalia/vs_prot_syn/test_genes_top10.txt"
with open(gene_list_file, "w") as f:
    for gene in sorted(all_genes_with_neighborhoods):
        f.write(f"{gene}\n")
print(f"Saved gene list to: {gene_list_file}")

# -------------------------------------------------------------
# 2. Filter GTF files
# -------------------------------------------------------------

def filter_gtf_optimized(input_gtf, output_gtf, gene_set):
    """
    Optimized GTF filtering for large files.
    Uses compiled regex for fast pattern matching.
    """
    print(f"\n  Filtering: {Path(input_gtf).name}")
    print(f"  Searching for {len(gene_set):,} genes...")
    
    start_time = time.time()
    
    # Create single regex pattern: GeneID:(gene1|gene2|gene3|...)
    # This is MUCH faster than checking each gene individually
    gene_pattern = '|'.join(re.escape(str(g)) for g in gene_set)
    regex = re.compile(f'GeneID:(?:{gene_pattern})')
    
    kept_lines = 0
    total_lines = 0
    
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            total_lines += 1
            
            # Progress every 1M lines
            if total_lines % 1000000 == 0:
                elapsed = time.time() - start_time
                rate = total_lines / elapsed
                print(f"    {total_lines:,} lines ({kept_lines:,} kept) | {rate:,.0f} lines/sec")
            
            # Keep comment lines
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            # Fast regex search
            if regex.search(line):
                outfile.write(line)
                kept_lines += 1
    
    elapsed = time.time() - start_time
    print(f"  ✓ Done in {elapsed:.1f}s: {total_lines:,} total → {kept_lines:,} kept")
    return kept_lines


# Define your GTF file paths
gtf_files = {
    "can": "/storage/EasyVectorOmics/synteny_algorithm/material/mammalian/gtf_Canis.gtf",
    "mac": "/storage/EasyVectorOmics/synteny_algorithm/material/mammalian/gtf_Macaca.gtf",
    "rat": "/storage/EasyVectorOmics/synteny_algorithm/material/mammalian/gtf_Rat.gtf",
    "mus": "/storage/EasyVectorOmics/synteny_algorithm/material/mammalian/gtf_Mus.gtf",
}


# Output directory
output_dir = Path("results/evaluation/mammalia/vs_prot_syn/filtered_gtfs_top10")
output_dir.mkdir(exist_ok=True)

# -------------------------------------------------------------
# Run the filtering with timing
# -------------------------------------------------------------

print("\n=== Filtering GTF files (OPTIMIZED) ===")
print(f"Gene set size: {len(all_genes_with_neighborhoods):,}")

total_start = time.time()

for species, gtf_path in gtf_files.items():
    output_path = output_dir / f"{species}_filtered.gtf"
    
    if not Path(gtf_path).exists():
        print(f"\n  SKIPPING {species}: File not found at {gtf_path}")
        continue
    
    file_size_mb = Path(gtf_path).stat().st_size / (1024**2)
    print(f"\n[{species}] File size: {file_size_mb:.1f} MB")
    
    kept = filter_gtf_optimized(gtf_path, output_path, all_genes_with_neighborhoods)

total_elapsed = time.time() - total_start
print(f"\n{'='*60}")
print(f" ALL FILES DONE in {total_elapsed/60:.1f} minutes")
print(f" Filtered GTFs saved to: {output_dir}")

# -------------------------------------------------------------
# 3. Verify the filtered files
# -------------------------------------------------------------

print("\n=== Verification ===")
for species, gtf_path in gtf_files.items():
    filtered_path = output_dir / f"{species}_filtered.gtf"
    
    # Count unique genes in filtered file
    genes_found = set()
    with open(filtered_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            # Extract GeneID
            if 'GeneID:' in line:
                start = line.find('GeneID:') + 7
                end = line.find('"', start)
                gene_id = line[start:end]
                genes_found.add(gene_id)
    
    print(f"{species}: {len(genes_found)} unique genes")

# Check which genes from our list are in each species
print("\n=== Genes by species ===")
species_gene_map = {
    "can": set(map_df[map_df["Species"].str.startswith("Canis")]["GeneID"].astype(str)),
    "mac": set(map_df[map_df["Species"].str.startswith("Macaca")]["GeneID"].astype(str)),
    "mus": set(map_df[map_df["Species"].str.startswith("Mus")]["GeneID"].astype(str)),
    "rat": set(map_df[map_df["Species"].str.startswith("Rattus")]["GeneID"].astype(str)),
}

for species, all_species_genes in species_gene_map.items():
    overlap = all_genes_with_neighborhoods & all_species_genes
    print(f"{species}: {len(overlap)} genes from our list")