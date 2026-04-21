#!/bin/bash
# Name of the job
#SBATCH --job-name=gtf_to_bed
#SBATCH --output=slurm_outputs/gtf_to_bed%j.out

# Convert GTF files to BED for mammals with convert2bed from bedops

echo "Starting the script"

mkdir -p results/gene_only_annotation/mammalia/

convert_gtf_to_bed() {
  local input_filtered_gtf=$1
  local output_bed=$2
  local output_bed6=$3

  # Convert GTF to BED
  convert2bed --input=gtf < $input_filtered_gtf > $output_bed

  # Simplify to BED6 (chrom, start, end, ID, score=0, strand)
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,0,$6}' \
      $output_bed \
      > $output_bed6

  # Ensure no stray strand values are in the column 6 of the BED6 file (should print 0 for each species)
  awk '$6!="+" && $6!="-" && $6!="."' $output_bed6 | wc -l
}


# Mammals
#--------

# Canis
convert_gtf_to_bed "results/gene_only_annotation/mammalia/Canis_genes_only_fix.gtf" \
                   "results/gene_only_annotation/mammalia/Canis_genes_only_fix.bed" \
                   "results/gene_only_annotation/mammalia/Canis_genes_only_fix.bed6"

# Macaca
convert_gtf_to_bed "results/gene_only_annotation/mammalia/Macaca_genes_only_fix.gtf" \
                   "results/gene_only_annotation/mammalia/Macaca_genes_only_fix.bed" \
                   "results/gene_only_annotation/mammalia/Macaca_genes_only_fix.bed6"

# Rat
convert_gtf_to_bed "results/gene_only_annotation/mammalia/Rat_genes_only_fix.gtf" \
                    "results/gene_only_annotation/mammalia/Rat_genes_only_fix.bed" \
                    "results/gene_only_annotation/mammalia/Rat_genes_only_fix.bed6"

# Mus
convert_gtf_to_bed "results/gene_only_annotation/mammalia/Mus_genes_only_fix.gtf" \
                    "results/gene_only_annotation/mammalia/Mus_genes_only_fix.bed" \
                    "results/gene_only_annotation/mammalia/Mus_genes_only_fix.bed6"


echo "Finished script"
