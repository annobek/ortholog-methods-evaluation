#!/bin/bash
# Name of the job
#SBATCH --job-name=run_linear_synteny
#SBATCH --output=slurm_outputs/run_linear_synteny%j.out
#SBATCH --output=slurm_outputs/run_linear_synteny%j.err

echo "Starting the script"

python3 -B gene_relationship_classifier/methods/run_linear_synteny.py

echo "Finished script"