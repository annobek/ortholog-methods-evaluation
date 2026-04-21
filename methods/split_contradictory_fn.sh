#!/bin/bash
# Name of the job
#SBATCH --job-name=split_contradictory_fn
#SBATCH --output=slurm_outputs/split_contradictory_fn%j.out
#SBATCH --output=slurm_outputs/split_contradictory_fn%j.err

# split FN contradictory (gene-level) into focus_gene-cactus_ortholog and focus_gene-ps_ortholog
# output: two TSV tables with split data + map of focus_genes to their corresponding cactus and ortholog partners

echo "Starting the script"

# prot-syn nscore
echo "Starting nscore approach."

python3 gene_relationship_classifier/methods/split_contradictory_fn.py \
  --fn results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/false_negatives/false_negatives_classified_gene_level.tsv \
  --out_prefix results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/false_negatives/contradictory

echo "Finished nscore approach."
echo "...................."

# prot-syn sum
echo "Starting sum approach."

python3 gene_relationship_classifier/methods/split_contradictory_fn.py \
  --fn results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_negatives/false_negatives_classified_gene_level.tsv \
  --out_prefix results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_negatives/contradictory

echo "Finished sum approach."

echo "Finished script"