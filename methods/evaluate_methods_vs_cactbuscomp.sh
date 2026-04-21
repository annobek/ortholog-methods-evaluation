#!/bin/bash
#SBATCH --job-name=evaluate_methods_vs_cactbuscomp_set
#SBATCH --output=slurm_outputs/evaluate_methods_vs_cactbuscomp_set%j.out

# Evaluate algorithms against CactBUSComp-set (cactus+BUSCO/compleasm)

echo "Starting the script"

python3 gene_relationship_classifier/methods/evaluate_methods_vs_cactbuscomp_v2.py

echo "Finished script"


