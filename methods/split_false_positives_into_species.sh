#!/bin/bash
# Name of the job
#SBATCH --job-name=split_false_positives_into_species
#SBATCH --output=slurm_outputs/split_false_positives_into_species%j.out

# split false positives file obtained after evaluation into species pairs


echo "Starting the script"

python3 gene_relationship_classifier/methods/split_false_positives_into_species.py

echo "Finished script"
