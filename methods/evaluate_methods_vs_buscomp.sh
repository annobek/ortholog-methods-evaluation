#!/bin/bash
#SBATCH --job-name=evaluate_methods_vs_buscomp
#SBATCH --output=slurm_outputs/evaluate_methods_vs_buscomp%j.out

# Evaluate BUSCO+compleasm (BUSComp-set) against methdos

echo "Starting the script"

python3 gene_relationship_classifier/methods/evaluate_methods_vs_buscomp_v2.py

echo "Finished script"

