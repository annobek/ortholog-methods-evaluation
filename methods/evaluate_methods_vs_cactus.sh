#!/bin/bash
#SBATCH --job-name=evaluate_methods_vs_cactus_v2_
#SBATCH --output=slurm_outputs/evaluate_methods_vs_cactus_v2_%j.out

# Evaluate algorithms against cactus set

echo "Starting the script"

python3 gene_relationship_classifier/methods/evaluate_methods_vs_cactus_v2.py

echo "Finished script"

