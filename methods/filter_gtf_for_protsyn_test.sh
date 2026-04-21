#!/bin/bash
# Name of the job
#SBATCH --job-name=filter_gtf_for_protsyn_test
#SBATCH --output=slurm_outputs/filter_gtf_for_protsyn_test%j.out


echo "Starting the script"

python3 gene_relationship_classifier/methods/filter_gtf_for_protsyn_test.py

echo "Finished script"
