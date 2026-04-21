#!/bin/bash
# Name of the job
#SBATCH --job-name=plot_cactus_score_contradictory_fn
#SBATCH --output=slurm_outputs/plot_cactus_score_contradictory_fn%j.out
#SBATCH --output=slurm_outputs/plot_cactus_score_contradictory_fn%j.err

# nscore
python3 gene_relationship_classifier/methods/plot_cactus_score_contradictory_fn.py \
  --focus_map \
    results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_negatives/contradictory_focus_candidate_map.tsv \
  --c_scores \
    results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_negatives/contradictory_c_pairs_with_cactus_score.tsv \
  --n_scores \
    results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_negatives/contradictory_n_pairs_with_cactus_score.tsv \
  --out_tsv \
    results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_negatives/contradictory_cactus_support_scores_merged.tsv \
  --out_png \
    results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_negatives/contradictory_cactus_support_scatter.png
