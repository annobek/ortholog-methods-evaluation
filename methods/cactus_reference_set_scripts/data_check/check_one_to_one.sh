#!/bin/bash
#SBATCH --job-name=check_one_to_one
#SBATCH --output=slurm_outputs/check_one_to_one%j.out

# Verify 1:1 relationship

echo "Starting the script"

python3 gene_relationship_classifier/methods/data_check/check_one_to_one.py

echo "Finished script"


