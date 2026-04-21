#!/bin/bash
#SBATCH --job-name=check_in_paralogs
#SBATCH --output=slurm_outputs/check_in_paralogs%j.out

# Verify in-paralogs

echo "Starting the script"

python3 gene_relationship_classifier/methods/data_check/check_in_paralogs.py

echo "Finished script"

