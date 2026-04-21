#!/bin/bash
# Name of the job
#SBATCH --job-name=plot_cactus_score_no_result_fp
#SBATCH --output=slurm_outputs/plot_cactus_score_no_result_fp%j.out
#SBATCH --output=slurm_outputs/plot_cactus_score_no_result_fp%j.err

echo "Starting the script"

# nscore
echo "Plotting for nscore..."

python gene_relationship_classifier/methods/plot_cactus_score_no_result_fp.py \
    --input results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/false_positives/pairs_with_cactus_score.tsv \
    --output results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/false_positives/cactus_score_distribution_no_result.png

echo "Finished plotting for nscore."

echo "------------------------------------"

# sum
echo "Plotting for sum..."

python gene_relationship_classifier/methods/plot_cactus_score_no_result_fp.py \
    --input results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_positives/pairs_with_cactus_score.tsv \
    --output results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_positives/cactus_score_distribution_no_result.png

echo "Finished plotting for sum."


echo "Finished the script"
