#!/bin/bash
#SBATCH --job-name=check_one_to_one_orthofinder
#SBATCH --output=slurm_outputs/check_one_to_one_orthofinder%j.out

# Verify 1:1 relationship in orthofinder results

echo "Starting the script"

python3 gene_relationship_classifier/methods/data_check/check_one_to_one_orthofinder.py

echo "Finished script"

