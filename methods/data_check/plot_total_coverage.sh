#!/bin/bash
# Name of the job
#SBATCH --job-name=plot_total_coverage
#SBATCH --output=slurm_outputs/plot_total_coverage%j.out

# Plot total coverage

echo "Starting the script"

python3 gene_relationship_classifier/methods/data_check/plot_total_coverage.py

echo "Finished script"

