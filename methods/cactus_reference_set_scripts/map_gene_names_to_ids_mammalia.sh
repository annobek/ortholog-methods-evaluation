#!/bin/bash
# Name of the job
#SBATCH --job-name=map_names_id
#SBATCH --output=slurm_outputs/map_names_id%j.out

# Apply gene name -> gene ID mapping to BED6 files

echo "Starting the script"

mkdir -p results/gene_only_annotation/mammalia/

# Gene-only 

map_names_id_genes() {
  local input_map=$1
  local input_bed6=$2
  local output_stable_bed6=$3

  awk 'BEGIN{OFS="\t"}
      FNR==NR {map[$1] = $2; next;}
      {if ($4 in map) $4 = map[$4]; print;}' \
      $input_map $input_bed6 \
      > $output_stable_bed6

  echo "Created mapped BED6: $output_stable_bed6"    
}

# Canis
map_names_id_genes "results/gene_only_annotation/mammalia/Canis_gene_to_id_fix.tsv" \
                   "results/gene_only_annotation/mammalia/Canis_genes_only_fix.bed6" \
                   "results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6"

# Macaca
map_names_id_genes "results/gene_only_annotation/mammalia/Macaca_gene_to_id_fix.tsv" \
                   "results/gene_only_annotation/mammalia/Macaca_genes_only_fix.bed6" \
                   "results/gene_only_annotation/mammalia/Macaca_genes_only_stableid_fix.bed6"

# Rat
map_names_id_genes "results/gene_only_annotation/mammalia/Rat_gene_to_id_fix.tsv" \
                   "results/gene_only_annotation/mammalia/Rat_genes_only_fix.bed6" \
                   "results/gene_only_annotation/mammalia/Rat_genes_only_stableid_fix.bed6"

# Mus
map_names_id_genes "results/gene_only_annotation/mammalia/Mus_gene_to_id_fix.tsv" \
                   "results/gene_only_annotation/mammalia/Mus_genes_only_fix.bed6" \
                   "results/gene_only_annotation/mammalia/Mus_genes_only_stableid_fix.bed6"


# CDS-only
#awk 'BEGIN{OFS="\t"}
#FNR==NR {map[$1] = $2; next;}
#{if ($4 in map) $4 = map[$4]; print;}' \
#results/ann_cds/Canis_cds_to_id.txt results/ann_cds/Canis_cds_only.bed6 \
#> results/ann_cds/Canis_cds_only_stableid.bed6


echo "Finished!"
