#!/bin/bash
# Name of the job
#SBATCH --job-name=plot_total_number_genes
#SBATCH --output=slurm_outputs/plot_total_number_genes%j.out

# Plot piechart of ratio of total number of genes and genes located in synteny blocks (drosophila and plants)

echo "Starting the script"

python3 gene_relationship_classifier/methods/data_check/plot_total_number_genes.py

echo "Finished script"

