#!/bin/bash
#SBATCH --job-name=evaluate_methods_vs_each_other_and_plot_v3
#SBATCH --output=slurm_outputs/evaluate_methods_vs_each_other_and_plot_v3%j.out

# Evaluate algorithms against each other (Jaccard similarity only) and plot heatmap

echo "Starting the script"

python3 gene_relationship_classifier/methods/evaluate_methods_vs_each_other_and_plot_v3.py

echo "Finished script"

