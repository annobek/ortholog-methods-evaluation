#!/bin/bash
#SBATCH --job-name=extract_one_to_one_orthofinder
#SBATCH --output=slurm_outputs/extract_one_to_one_orthofinder%j.out

# Extract 1:1 pairs from orthofinder

echo "Starting the script"

python3 gene_relationship_classifier/methods/extract_one_to_one_orthofinder.py

echo "Finished script"

