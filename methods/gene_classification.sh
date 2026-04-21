#!/bin/bash
#SBATCH --job-name=gene_classification
#SBATCH --output=slurm_outputs/gene_classification%j.out

# analyze discarded genes after 1:1 step and classify 

echo "Starting the script"

python3 gene_relationship_classifier/methods/gene_classification.py

echo "Finished script"
