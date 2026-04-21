#!/bin/bash
# Name of the job
#SBATCH --job-name=plot_region_conservativeness
#SBATCH --output=slurm_outputs/plot_region_conservativeness%j.out

# Plot region conservativeness

echo "Starting the script"

python3 gene_relationship_classifier/methods/data_check/plot_region_conservativeness.py

echo "Finished script"
