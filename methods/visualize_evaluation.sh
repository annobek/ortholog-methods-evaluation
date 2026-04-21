#!/bin/bash
#SBATCH --job-name=visualize_evaluation_metrics_v2
#SBATCH --output=slurm_outputs/visualize_evaluation_metrics_v2%j.out

# Visualize evaluation metrics of methods vs differnet reference sets - Cactus, CactBUSComp and BUSComp

echo "Starting the script"

python3 gene_relationship_classifier/methods/visualize_evaluation_metrics_v2.py

echo "Finished script"

