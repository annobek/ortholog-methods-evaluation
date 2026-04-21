#!/bin/bash
# Name of the job
#SBATCH --job-name=plot_total_number_genes_mammalia
#SBATCH --output=slurm_outputs/plot_total_number_genes_mammalia%j.out

# Plot piechart of ratio of total number of genes and genes located in synteny blocks

echo "Starting the script"

python3 gene_relationship_classifier/methods/data_check/plot_total_number_genes_mammalia.py

echo "Finished script"

