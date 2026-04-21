#!/bin/bash
# Name of the job
#SBATCH --job-name=max_chr_distance
#SBATCH --output=slurm_outputs/max_chr_distance%j.out


echo "Starting the script"

# prot-syn nscore
echo "Starting calculating and plotting for prot-syn nscore"

python gene_relationship_classifier/methods/max_chr_distance.py \
    --fn_table results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/false_negatives/false_negatives_classified_gene_level.tsv \
    --gene_gtf \
      material/sex_experiment/gtf_Canis.gtf \
      material/sex_experiment/gtf_Macaca.gtf \
      material/sex_experiment/gtf_Rat.gtf \
      material/sex_experiment/gtf_Mus.gtf \
    --neighborhood_dict /storage/EasyVectorOmics/synteny_algorithm/results/mammalian/neighborhoods.pkl \
    --grass_scores /storage/EasyVectorOmics/synteny_algorithm/results/mammalian/geometric_mean.pkl \
    --output_tsv results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/false_negatives/max_neighborhood_distances.tsv \
    --output_plot results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/disagreement_investigation/false_negatives/max_neighborhood_distances.png \
    --cactus_genes \
      results/genes_in_blocks/mammalia/region_gene_tables/can_mac_synteny_tables_genes_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/can_rat_synteny_tables_genes_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/can_mus_synteny_tables_genes_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/mac_rat_synteny_tables_genes_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/mac_mus_synteny_tables_genes_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/rat_mus_synteny_tables_genes_fix.tsv \
    --cactus_regions \
      results/genes_in_blocks/mammalia/region_gene_tables/can_mac_synteny_tables_regions_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/can_rat_synteny_tables_regions_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/can_mus_synteny_tables_regions_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/mac_rat_synteny_tables_regions_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/mac_mus_synteny_tables_regions_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/rat_mus_synteny_tables_regions_fix.tsv \

echo "Finished calculating and plotting for prot-syn nscore"

echo "------------------------------------------"

# prot-syn sum
echo "Starting calculating and plotting for prot-syn sum"
python gene_relationship_classifier/methods/max_chr_distance.py \
    --fn_table results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_negatives/false_negatives_classified_gene_level.tsv \
     --gene_gtf \
        material/sex_experiment/gtf_Canis.gtf \
        material/sex_experiment/gtf_Macaca.gtf \
        material/sex_experiment/gtf_Rat.gtf \
        material/sex_experiment/gtf_Mus.gtf \
    --neighborhood_dict /storage/EasyVectorOmics/synteny_algorithm/results/mammalian/neighborhoods.pkl \
    --grass_scores /storage/EasyVectorOmics/synteny_algorithm/results/mammalian/geometric_mean.pkl \
    --output_tsv results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_negatives/max_neighborhood_distances.tsv \
    --output_plot results/evaluation/mammalia/vs_prot_syn/sum_approach/test_no_transitivity/disagreement_investigation/false_negatives/max_neighborhood_distances.png \
    --cactus_genes \
      results/genes_in_blocks/mammalia/region_gene_tables/can_mac_synteny_tables_genes_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/can_rat_synteny_tables_genes_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/can_mus_synteny_tables_genes_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/mac_rat_synteny_tables_genes_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/mac_mus_synteny_tables_genes_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/rat_mus_synteny_tables_genes_fix.tsv \
    --cactus_regions \
      results/genes_in_blocks/mammalia/region_gene_tables/can_mac_synteny_tables_regions_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/can_rat_synteny_tables_regions_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/can_mus_synteny_tables_regions_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/mac_rat_synteny_tables_regions_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/mac_mus_synteny_tables_regions_fix.tsv \
      results/genes_in_blocks/mammalia/region_gene_tables/rat_mus_synteny_tables_regions_fix.tsv \

echo "Finished calculating and plotting for prot-syn sum"

echo "Finished the script"