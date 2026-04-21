#!/bin/bash
# Name of the job
#SBATCH --job-name=neighborhood_lengths
#SBATCH --output=slurm_outputs/neighborhood_lengths%j.out


echo "Starting the script"

# prot-syn nscore
echo "Starting calculating and plotting for prot-syn nscore"

#python gene_relationship_classifier/methods/neighborhood_lengths.py \
#    --fn results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/false_negatives/false_negatives_classified_gene_level.tsv \
#    --fp results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/false_positives/all_C_no_results.tsv \
#    --neighborhoods /storage/EasyVectorOmics/synteny_algorithm/results/mammalian/neighborhoods.pkl \
#    --gtf-can material/sex_experiment/gtf_Canis.gtf \
#    --gtf-mac material/sex_experiment/gtf_Macaca.gtf \
#    --gtf-rat material/sex_experiment/gtf_Rat.gtf \
#    --gtf-mus material/sex_experiment/gtf_Mus.gtf \
#    --output-plot results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/neighborhood_lengths_raincloud.png \
#    --output-summary results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/neighborhood_summary.tsv \
#    --output-detailed results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/neighborhood_lengths_detailed.tsv

#echo "Finished calculating and plotting for prot-syn nscore"

#echo "------------------------------------------"

# prot-syn sum
echo "Starting calculating and plotting for prot-syn sum"

python gene_relationship_classifier/methods/neighborhood_lengths.py \
    --fn results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_negatives/false_negatives_classified_gene_level.tsv \
    --fp results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_positives/all_C_no_results.tsv \
    --neighborhoods /storage/EasyVectorOmics/synteny_algorithm/results/mammalian/neighborhoods.pkl \
    --gtf-can material/sex_experiment/gtf_Canis.gtf \
    --gtf-mac material/sex_experiment/gtf_Macaca.gtf \
    --gtf-rat material/sex_experiment/gtf_Rat.gtf \
    --gtf-mus material/sex_experiment/gtf_Mus.gtf \
    --output-plot results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/neighborhood_lengths_raincloud.png \
    --output-summary results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/neighborhood_summary.tsv \
    --output-detailed results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/neighborhood_lengths_detailed.tsv    


echo "Finished calculating and plotting for prot-syn sum"

echo "Finished the script"    