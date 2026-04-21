#!/bin/bash
#SBATCH --job-name=reciprocity_and_resolve_one_to_many
#SBATCH --output=slurm_outputs/reciprocity_and_resolve_one_to_many%j.out

# Extracts 1:1 pairs

echo "Starting the script"

python3 gene_relationship_classifier/methods/reciprocity_and_resolve_one_to_many.py

echo "Finished script"

