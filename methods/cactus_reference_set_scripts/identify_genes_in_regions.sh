#!/bin/bash
# Name of the job
#SBATCH --job-name=find_genes_in_blocks
#SBATCH --output=slurm_outputs/find_genes_in_blocks%j.out

# Identify genes in synteny blocks with bedtools intersect

echo "Starting the script"

mkdir -p results/genes_in_blocks/mammalia/

find_genes() {

  local input_blocks_query=$1
  local input_blocks_target=$2
  local input_ann_query=$3
  local input_ann_target=$4
  local output_query=$5
  local output_target=$6
  local species1=$7
  local species2=$8

  echo ""
  echo "Processing $species1-$species2:"
  echo "---------------------"

  bedtools intersect -a $input_blocks_query -b $input_ann_query -wo > $output_query
  bedtools intersect -a $input_blocks_target -b $input_ann_target -wo > $output_target

  echo "Calculating statistics for $species1..."
  
  # Total synteny block coverage
  total_synteny_query=$(awk '{sum += ($3-$2)} END {print sum}' $input_blocks_query)

  # Synteny blocks containing no genes
  blocks_no_genes_query=$(awk 'FNR==NR {seen[$4]; next} !($4 in seen)' $output_query $input_blocks_query | wc -l)
  
  # Synteny blocks that contain at least one gene
  blocks_with_genes_query=$(cut -f1-4 $output_query | sort -u | wc -l)
  total_blocks_query=$(wc -l < $input_blocks_query)
  
  # Total genes from the annotation
  total_genes_query=$(wc -l < $input_ann_query)
  
  # Genes in syntenic regions
  genes_in_synteny_query=$(awk '{print $5"-"$6"-"$7"-"$8}' $output_query | sort -u | wc -l)
  
  # Gene coverage in syntenic regions
  gene_bases_in_synteny_query=$(awk '{sum += $NF} END {print sum+0}' $output_query)
  
  # Total gene bases
  total_gene_bases_query=$(awk '{sum += ($3-$2)} END {print sum}' $input_ann_query)
  
  # Calculate percentages
  pct_blocks_with_genes_query=$(echo "scale=2; $blocks_with_genes_query * 100 / $total_blocks_query" | bc)
  pct_genes_in_synteny_query=$(echo "scale=2; $genes_in_synteny_query * 100 / $total_genes_query" | bc)
  pct_gene_bases_syntenic_query=$(echo "scale=2; $gene_bases_in_synteny_query * 100 / $total_gene_bases_query" | bc)

  # Repeat for target species
  echo "Calculating statistics for $species2..."
  
  total_synteny_target=$(awk '{sum += ($3-$2)} END {print sum}' $input_blocks_target)
  blocks_no_genes_target=$(awk 'FNR==NR {seen[$4]; next} !($4 in seen)' $output_target $input_blocks_target | wc -l)
  blocks_with_genes_target=$(cut -f1-4 $output_target | sort -u | wc -l)
  total_blocks_target=$(wc -l < $input_blocks_target)
  total_genes_target=$(wc -l < $input_ann_target)
  genes_in_synteny_target=$(awk '{print $5"-"$6"-"$7"-"$8}' $output_target | sort -u | wc -l)
  gene_bases_in_synteny_target=$(awk '{sum += $NF} END {print sum+0}' $output_target)
  total_gene_bases_target=$(awk '{sum += ($3-$2)} END {print sum}' $input_ann_target)
  
  pct_blocks_with_genes_target=$(echo "scale=2; $blocks_with_genes_target * 100 / $total_blocks_target" | bc)
  pct_genes_in_synteny_target=$(echo "scale=2; $genes_in_synteny_target * 100 / $total_genes_target" | bc)
  pct_gene_bases_syntenic_target=$(echo "scale=2; $gene_bases_in_synteny_target * 100 / $total_gene_bases_target" | bc)

  # Print summary
  echo ""
  echo "SUMMARY for $species1-$species2:"
  echo "================================="
  echo ""
  echo "$species1 (Query):"
  echo "  Synteny blocks: $total_blocks_query"
  echo "  Blocks with no genes: $blocks_no_genes_query"
  echo "  Blocks with genes: $blocks_with_genes_query ($pct_blocks_with_genes_query%)"
  echo "  Total genes: $total_genes_query"
  echo "  Genes in synteny: $genes_in_synteny_query ($pct_genes_in_synteny_query%)"
  echo "  Gene bases in synteny: $(echo "scale=2; $gene_bases_in_synteny_query / 1000000" | bc) Mb ($pct_gene_bases_syntenic_query%)"
  echo ""
  echo "$species2 (Target):"
  echo "  Synteny blocks: $total_blocks_target"
  echo "  Blocks with no genes: $blocks_no_genes_target"
  echo "  Blocks with genes: $blocks_with_genes_target ($pct_blocks_with_genes_target%)"
  echo "  Total genes: $total_genes_target"
  echo "  Genes in synteny: $genes_in_synteny_target ($pct_genes_in_synteny_target%)"
  echo "  Gene bases in synteny: $(echo "scale=2; $gene_bases_in_synteny_target / 1000000" | bc) Mb ($pct_gene_bases_syntenic_target%)"
 
  echo "$species1-$species2 finished"
  echo "----------------------------------------------------------"

}


# Mammalia below


echo "Gene-Synteny Intersection Statistics" 
echo "====================================="

## Plants
##-------
#
## ar-card
#
#find_genes "results/halsynteny_output/brassicaceae/bed/ar_card_synteny_query_balanced.bed" \
#           "results/halsynteny_output/brassicaceae/bed/ar_card_synteny_target_balanced.bed" \
#           "results/gene_only_annotation/brassicaceae/Arabidopsis_genes_only.bed6" \
#           "results/gene_only_annotation/brassicaceae/Cardamine_genes_only.bed6" \
#           "results/genes_in_blocks/brassicaceae/ar_card_genes_in_blocks_query.bed" \
#           "results/genes_in_blocks/brassicaceae/ar_card_genes_in_blocks_target.bed" \
#           "ar" \
#           "card" 
#
## card-ar
#
#find_genes "results/halsynteny_output/brassicaceae/bed/card_ar_synteny_query_balanced.bed" \
#           "results/halsynteny_output/brassicaceae/bed/card_ar_synteny_target_balanced.bed" \
#           "results/gene_only_annotation/brassicaceae/Cardamine_genes_only.bed6" \
#           "results/gene_only_annotation/brassicaceae/Arabidopsis_genes_only.bed6" \
#           "results/genes_in_blocks/brassicaceae/card_ar_genes_in_blocks_query.bed" \
#           "results/genes_in_blocks/brassicaceae/card_ar_genes_in_blocks_target.bed" \
#           "card" \
#           "ar" 
#
#
## Drosophila
##-----------
#
## dmel-dsec
#
#find_genes "results/halsynteny_output/drosophila/bed/dmel_dsec_synteny_query_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dmel_dsec_synteny_target_new.bed" \
#           "results/gene_only_annotation/drosophila/dmel_genes_only.bed6" \
#           "results/gene_only_annotation/drosophila/dsec_genes_only.bed6" \
#           "results/genes_in_blocks/drosophila/dmel_dsec_genes_in_blocks_query.bed" \
#           "results/genes_in_blocks/drosophila/dmel_dsec_genes_in_blocks_target.bed" \
#           "dmel" \
#           "dsec" 
#
## dsec-dmel
#
#find_genes "results/halsynteny_output/drosophila/bed/dsec_dmel_synteny_query_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dsec_dmel_synteny_target_new.bed" \
#           "results/gene_only_annotation/drosophila/dsec_genes_only.bed6" \
#           "results/gene_only_annotation/drosophila/dmel_genes_only.bed6" \
#           "results/genes_in_blocks/drosophila/dsec_dmel_genes_in_blocks_query.bed" \
#           "results/genes_in_blocks/drosophila/dsec_dmel_genes_in_blocks_target.bed" \
#           "dsec" \
#           "dmel" 
#
#
## dmel-dsim
#
#find_genes "results/halsynteny_output/drosophila/bed/dmel_dsim_synteny_query_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dmel_dsim_synteny_target_new.bed" \
#           "results/gene_only_annotation/drosophila/dmel_genes_only.bed6" \
#           "results/gene_only_annotation/drosophila/dsim_genes_only.bed6" \
#           "results/genes_in_blocks/drosophila/dmel_dsim_genes_in_blocks_query.bed" \
#           "results/genes_in_blocks/drosophila/dmel_dsim_genes_in_blocks_target.bed" \
#           "dmel" \
#           "dsim" 
#
## dsim-dmel
#
#find_genes "results/halsynteny_output/drosophila/bed/dsim_dmel_synteny_query_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dsim_dmel_synteny_target_new.bed" \
#           "results/gene_only_annotation/drosophila/dsim_genes_only.bed6" \
#           "results/gene_only_annotation/drosophila/dmel_genes_only.bed6" \
#           "results/genes_in_blocks/drosophila/dsim_dmel_genes_in_blocks_query.bed" \
#           "results/genes_in_blocks/drosophila/dsim_dmel_genes_in_blocks_target.bed" \
#           "dsim" \
#           "dmel" 
#
## dsec-dsim
#
#find_genes "results/halsynteny_output/drosophila/bed/dsec_dsim_synteny_query_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dsec_dsim_synteny_target_new.bed" \
#           "results/gene_only_annotation/drosophila/dsec_genes_only.bed6" \
#           "results/gene_only_annotation/drosophila/dsim_genes_only.bed6" \
#           "results/genes_in_blocks/drosophila/dsec_dsim_genes_in_blocks_query.bed" \
#           "results/genes_in_blocks/drosophila/dsec_dsim_genes_in_blocks_target.bed" \
#           "dsec" \
#           "dsim" 
#
## dsim-dsec
#
#find_genes "results/halsynteny_output/drosophila/bed/dsim_dsec_synteny_query_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dsim_dsec_synteny_target_new.bed" \
#           "results/gene_only_annotation/drosophila/dsim_genes_only.bed6" \
#           "results/gene_only_annotation/drosophila/dsec_genes_only.bed6" \
#           "results/genes_in_blocks/drosophila/dsim_dsec_genes_in_blocks_query.bed" \
#           "results/genes_in_blocks/drosophila/dsim_dsec_genes_in_blocks_target.bed" \
#           "dsim" \
#           "dsec"
#
#
# Mammals
#--------

# can-mac

find_genes "results/halsynteny_output/mammalia/bed/can_mac_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/can_mac_synteny_blocks_50kb_target_fix.bed" \
           "results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6" \
           "results/gene_only_annotation/mammalia/Macaca_genes_only_stableid_fix.bed6" \
           "results/genes_in_blocks/mammalia/can_mac_genes_in_blocks_50kb_query_fix.bed" \
           "results/genes_in_blocks/mammalia/can_mac_genes_in_blocks_50kb_target_fix.bed" \
           "can" \
           "mac" 


# mac-can

find_genes "results/halsynteny_output/mammalia/bed/mac_can_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mac_can_synteny_blocks_50kb_target_fix.bed" \
           "results/gene_only_annotation/mammalia/Macaca_genes_only_stableid_fix.bed6" \
           "results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6" \
           "results/genes_in_blocks/mammalia/mac_can_genes_in_blocks_50kb_query_fix.bed" \
           "results/genes_in_blocks/mammalia/mac_can_genes_in_blocks_50kb_target_fix.bed" \
           "mac" \
           "can" 

# can-rat

find_genes "results/halsynteny_output/mammalia/bed/can_rat_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/can_rat_synteny_blocks_50kb_target_fix.bed" \
           "results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6" \
           "results/gene_only_annotation/mammalia/Rat_genes_only_stableid_fix.bed6" \
           "results/genes_in_blocks/mammalia/can_rat_genes_in_blocks_50kb_query_fix.bed" \
           "results/genes_in_blocks/mammalia/can_rat_genes_in_blocks_50kb_target_fix.bed" \
           "can" \
           "rat" 

# rat-can

find_genes "results/halsynteny_output/mammalia/bed/rat_can_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/rat_can_synteny_blocks_50kb_target_fix.bed" \
           "results/gene_only_annotation/mammalia/Rat_genes_only_stableid_fix.bed6" \
           "results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6" \
           "results/genes_in_blocks/mammalia/rat_can_genes_in_blocks_50kb_query_fix.bed" \
           "results/genes_in_blocks/mammalia/rat_can_genes_in_blocks_50kb_target_fix.bed" \
           "rat" \
           "can" 

# can-mus

find_genes "results/halsynteny_output/mammalia/bed/can_mus_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/can_mus_synteny_blocks_50kb_target_fix.bed" \
           "results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6" \
           "results/gene_only_annotation/mammalia/Mus_genes_only_stableid_fix.bed6" \
           "results/genes_in_blocks/mammalia/can_mus_genes_in_blocks_50kb_query_fix.bed" \
           "results/genes_in_blocks/mammalia/can_mus_genes_in_blocks_50kb_target_fix.bed" \
           "can" \
           "mus" 

# mus-can

find_genes "results/halsynteny_output/mammalia/bed/mus_can_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mus_can_synteny_blocks_50kb_target_fix.bed" \
           "results/gene_only_annotation/mammalia/Mus_genes_only_stableid_fix.bed6" \
           "results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6" \
           "results/genes_in_blocks/mammalia/mus_can_genes_in_blocks_50kb_query_fix.bed" \
           "results/genes_in_blocks/mammalia/mus_can_genes_in_blocks_50kb_target_fix.bed" \
           "mus" \
           "can" 


# mac-rat

find_genes "results/halsynteny_output/mammalia/bed/mac_rat_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mac_rat_synteny_blocks_50kb_target_fix.bed" \
           "results/gene_only_annotation/mammalia/Macaca_genes_only_stableid_fix.bed6" \
           "results/gene_only_annotation/mammalia/Rat_genes_only_stableid_fix.bed6" \
           "results/genes_in_blocks/mammalia/mac_rat_genes_in_blocks_50kb_query_fix.bed" \
           "results/genes_in_blocks/mammalia/mac_rat_genes_in_blocks_50kb_target_fix.bed" \
           "mac" \
           "rat" 

# rat-mac

find_genes "results/halsynteny_output/mammalia/bed/rat_mac_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/rat_mac_synteny_blocks_50kb_target_fix.bed" \
           "results/gene_only_annotation/mammalia/Rat_genes_only_stableid_fix.bed6" \
           "results/gene_only_annotation/mammalia/Macaca_genes_only_stableid_fix.bed6" \
           "results/genes_in_blocks/mammalia/rat_mac_genes_in_blocks_50kb_query_fix.bed" \
           "results/genes_in_blocks/mammalia/rat_mac_genes_in_blocks_50kb_target_fix.bed" \
           "rat" \
           "mac" 

# mac-mus

find_genes "results/halsynteny_output/mammalia/bed/mac_mus_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mac_mus_synteny_blocks_50kb_target_fix.bed" \
           "results/gene_only_annotation/mammalia/Macaca_genes_only_stableid_fix.bed6" \
           "results/gene_only_annotation/mammalia/Mus_genes_only_stableid_fix.bed6" \
           "results/genes_in_blocks/mammalia/mac_mus_genes_in_blocks_50kb_query_fix.bed" \
           "results/genes_in_blocks/mammalia/mac_mus_genes_in_blocks_50kb_target_fix.bed" \
           "mac" \
           "mus" 

# mus-mac

find_genes "results/halsynteny_output/mammalia/bed/mus_mac_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mus_mac_synteny_blocks_50kb_target_fix.bed" \
           "results/gene_only_annotation/mammalia/Mus_genes_only_stableid_fix.bed6" \
           "results/gene_only_annotation/mammalia/Macaca_genes_only_stableid_fix.bed6" \
           "results/genes_in_blocks/mammalia/mus_mac_genes_in_blocks_50kb_query_fix.bed" \
           "results/genes_in_blocks/mammalia/mus_mac_genes_in_blocks_50kb_target_fix.bed" \
           "mus" \
           "mac" 

# rat-mus

find_genes "results/halsynteny_output/mammalia/bed/rat_mus_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/rat_mus_synteny_blocks_50kb_target_fix.bed" \
           "results/gene_only_annotation/mammalia/Rat_genes_only_stableid_fix.bed6" \
           "results/gene_only_annotation/mammalia/Mus_genes_only_stableid_fix.bed6" \
           "results/genes_in_blocks/mammalia/rat_mus_genes_in_blocks_50kb_query_fix.bed" \
           "results/genes_in_blocks/mammalia/rat_mus_genes_in_blocks_50kb_target_fix.bed" \
           "rat" \
           "mus" 

# mus-rat

find_genes "results/halsynteny_output/mammalia/bed/mus_rat_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mus_rat_synteny_blocks_50kb_target_fix.bed" \
           "results/gene_only_annotation/mammalia/Mus_genes_only_stableid_fix.bed6" \
           "results/gene_only_annotation/mammalia/Rat_genes_only_stableid_fix.bed6" \
           "results/genes_in_blocks/mammalia/mus_rat_genes_in_blocks_50kb_query_fix.bed" \
           "results/genes_in_blocks/mammalia/mus_rat_genes_in_blocks_50kb_target_fix.bed" \
           "mus" \
           "rat" 

echo ""
echo "Finished the script"
