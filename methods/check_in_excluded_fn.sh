#!/bin/bash
# Name of the job
#SBATCH --job-name=check_in_excluded_fn
#SBATCH --output=slurm_outputs/check_in_excluded_fn%j.out

echo "Starting the script"
python3 gene_relationship_classifier/methods/check_in_excluded_fn.py \
   --fn results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/false_negatives/false_negatives_classified_gene_level.tsv \
   --excluded-dog material/sex_experiment/excluded_genes_Canis_Lupus.tsv \
   --excluded-rat material/sex_experiment/excluded_genes_Rat.tsv \
   --excluded-mouse material/sex_experiment/excluded_genes_Mus.tsv \
   --out-dir results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/false_negatives/

echo "Finished script"

