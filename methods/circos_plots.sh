#!/bin/bash
#SBATCH --job-name=circos_plots
#SBATCH --output=slurm_outputs/circos_plots%j.out
#SBATCH --mem=8G

# Create Circos plots for synteny visualization

echo "Starting the script"



circos -conf results/circos_plots/can_mac.conf  # For input and output information see .conf files 

echo "Finished script"