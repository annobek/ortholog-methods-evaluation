#!/bin/bash
#SBATCH --job-name=calculate_region_gris
#SBATCH --output=slurm_outputs/calculate_region_gris%j.out

# Calculate region GRIS scores based on protein similarity

echo "Starting the script"

python3 gene_relationship_classifier/methods/calculate_region_gris.py

echo "Finished script"