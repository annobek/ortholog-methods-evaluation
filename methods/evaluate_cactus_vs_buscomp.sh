#!/bin/bash
#SBATCH --job-name=evaluate_cactus_vs_buscomp
#SBATCH --output=slurm_outputs/evaluate_cactus_vs_buscomp%j.out

# Evaluate BUSCO+compleasm (BUSComp) against cactus and create union set CactBUSComp

echo "Starting the script"

python3 gene_relationship_classifier/methods/evaluate_cactus_vs_buscomp_v2.py

echo "Finished script"

