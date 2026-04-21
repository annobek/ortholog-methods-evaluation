#!/bin/bash
#SBATCH --job-name=neighborhood_score_cactus_fn_v4
#SBATCH --output=slurm_outputs/neighborhood_score_cactus_fn_v4_%j.out

# Investigate false negatives; calculate N-score; investigate contradictory results

echo "Starting the script"

python3 gene_relationship_classifier/methods/neighborhood_score_cactus_fn_v4.py

echo "Finished script"

