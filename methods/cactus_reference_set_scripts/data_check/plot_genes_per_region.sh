#!/bin/bash
# Name of the job
#SBATCH --job-name=plot_genes_per_region
#SBATCH --output=slurm_outputs/plot_genes_per_region%j.out

# Plot gene distribution in synteny regions

echo "Starting the script"

python3 gene_relationship_classifier/methods/data_check/plot_genes_per_region.py

echo "Finished script"
