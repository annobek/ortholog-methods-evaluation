#!/bin/bash
# Name of the job
#SBATCH --job-name=find_c_no_results_fp
#SBATCH --output=slurm_outputs/find_c_no_results_fp%j.out

echo "Starting the script"

python3 gene_relationship_classifier/methods/find_c_no_results_fp.py

echo "Finished script"
