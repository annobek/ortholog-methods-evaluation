#!/bin/bash
# Name of the job
#SBATCH --job-name=cactus_support_score_fp
#SBATCH --output=slurm_outputs/cactus_support_score_fp%j.out
#SBATCH --error=slurm_outputs/cactus_support_score_fp%j.err


echo "Starting the script"

python3 gene_relationship_classifier/methods/cactus_support_score_fp.py

echo "Finished script"
