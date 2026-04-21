#!/bin/bash
# Name of the job
#SBATCH --job-name=find_ps_no_results_and_contradictory_v2
#SBATCH --output=slurm_outputs/find_ps_no_results_and_contradictory_v2_%j.out

echo "Starting the script"

python3 gene_relationship_classifier/methods/find_ps_no_results_and_contradictory_v2.py

echo "Finished script"
