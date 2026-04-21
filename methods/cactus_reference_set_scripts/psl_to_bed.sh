#!/bin/bash
# Name of the job
#SBATCH --job-name=psl_to_bed
#SBATCH --output=slurm_outputs/psl_to_bed%j.out

# Convert PSL files (halSynteny results) to BED and add ID for each synteny block

echo "Starting the script"

mkdir -p results/halsynteny_output/mammalia/bed/

psl_to_bed() {
    local psl_input=$1
    local bed_output_all=$2
    local bed_output_query=$3
    local bed_output_target=$4
    
    echo "Processing $psl_input ..."

    awk -v bed_all="$bed_output_all" \
        -v bed_query="$bed_output_query" \
        -v bed_target="$bed_output_target" \
    'BEGIN { OFS="\t" }
     {
         block_id = "block_" NR
         # full table
         print block_id, $10, $12, $13, $14, $16, $17 > bed_all
         # query and target sides
         print $10, $12, $13, block_id > bed_query
         print $14, $16, $17, block_id > bed_target
     }' "$psl_input"

    echo "Processed and saved: $bed_output_query and $bed_output_target"

    # Check the number of blocks in query and target
    echo "Number of blocks in query"
    awk '{print $4}' "$bed_output_query" | sort | uniq -c | wc -l

    echo "Number of blocks in target"
    awk '{print $4}' "$bed_output_target" | sort | uniq -c | wc -l

    # Check if there are any differences (expected no output)
    echo "Differences between query and target file:"
    diff <(cut -f4 "$bed_output_query" | sort) <(cut -f4 "$bed_output_target" | sort)
}



# Mammalia below



##--------------------------------
## Brassicaceae
##--------------------------------
#
## ar-card
#
#psl_to_bed "results/halsynteny_output/brassicaceae/ar_card_synteny_balanced.psl" \
#           "results/halsynteny_output/brassicaceae/bed/ar_card_synteny_all_balanced.bed" \
#           "results/halsynteny_output/brassicaceae/bed/ar_card_synteny_query_balanced.bed" \
#           "results/halsynteny_output/brassicaceae/bed/ar_card_synteny_target_balanced.bed"
#
## card-ar
#
#psl_to_bed "results/halsynteny_output/brassicaceae/card_ar_synteny_balanced.psl" \
#           "results/halsynteny_output/brassicaceae/bed/card_ar_synteny_all_balanced.bed" \
#           "results/halsynteny_output/brassicaceae/bed/card_ar_synteny_query_balanced.bed" \
#           "results/halsynteny_output/brassicaceae/bed/card_ar_synteny_target_balanced.bed"
#
#
##-----------------------------
## Drosophila
##-----------------------------
#
## dmel-dsec
#psl_to_bed "results/halsynteny_output/drosophila/dmel_dsec_synteny_new.psl" \
#           "results/halsynteny_output/drosophila/bed/dmel_dsec_synteny_all_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dmel_dsec_synteny_query_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dmel_dsec_synteny_target_new.bed"
#
## dsec-dmel
#
#psl_to_bed "results/halsynteny_output/drosophila/dsec_dmel_synteny_new.psl" \
#           "results/halsynteny_output/drosophila/bed/dsec_dmel_synteny_all_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dsec_dmel_synteny_query_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dsec_dmel_synteny_target_new.bed"
#
## dmel-dsim
#
#psl_to_bed "results/halsynteny_output/drosophila/dmel_dsim_synteny_new.psl" \
#           "results/halsynteny_output/drosophila/bed/dmel_dsim_synteny_all_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dmel_dsim_synteny_query_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dmel_dsim_synteny_target_new.bed"
#
## dsim-dmel
#
#psl_to_bed "results/halsynteny_output/drosophila/dsim_dmel_synteny_new.psl" \
#           "results/halsynteny_output/drosophila/bed/dsim_dmel_synteny_all_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dsim_dmel_synteny_query_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dsim_dmel_synteny_target_new.bed"
#
## dsec-dsim
#
#psl_to_bed "results/halsynteny_output/drosophila/dsec_dsim_synteny_new.psl" \
#           "results/halsynteny_output/drosophila/bed/dsec_dsim_synteny_all_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dsec_dsim_synteny_query_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dsec_dsim_synteny_target_new.bed"
#
## dsim-dsec
#
#psl_to_bed "results/halsynteny_output/drosophila/dsim_dsec_synteny_new.psl" \
#           "results/halsynteny_output/drosophila/bed/dsim_dsec_synteny_all_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dsim_dsec_synteny_query_new.bed" \
#           "results/halsynteny_output/drosophila/bed/dsim_dsec_synteny_target_new.bed"
#
#
#---------------------------------------------
# Mammalia
#---------------------------------------------


# can-mac

psl_to_bed "results/halsynteny_output/mammalia/can_mac_synteny_50kb.psl" \
           "results/halsynteny_output/mammalia/bed/can_mac_synteny_blocks_50kb_all_fix.bed" \
           "results/halsynteny_output/mammalia/bed/can_mac_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/can_mac_synteny_blocks_50kb_target_fix.bed" 

# mac-can

psl_to_bed "results/halsynteny_output/mammalia/mac_can_synteny_50kb.psl" \
	       "results/halsynteny_output/mammalia/bed/mac_can_synteny_blocks_50kb_all_fix.bed" \
   	       "results/halsynteny_output/mammalia/bed/mac_can_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mac_can_synteny_blocks_50kb_target_fix.bed" 

# can-rat

psl_to_bed "results/halsynteny_output/mammalia/can_rat_synteny_50kb.psl" \
	       "results/halsynteny_output/mammalia/bed/can_rat_synteny_blocks_50kb_all_fix.bed" \
   	       "results/halsynteny_output/mammalia/bed/can_rat_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/can_rat_synteny_blocks_50kb_target_fix.bed" 

# rat-can

psl_to_bed "results/halsynteny_output/mammalia/rat_can_synteny_50kb.psl" \
	       "results/halsynteny_output/mammalia/bed/rat_can_synteny_blocks_50kb_all_fix.bed" \
	       "results/halsynteny_output/mammalia/bed/rat_can_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/rat_can_synteny_blocks_50kb_target_fix.bed" 


# can-mus

psl_to_bed "results/halsynteny_output/mammalia/can_mus_synteny_50kb.psl" \
	       "results/halsynteny_output/mammalia/bed/can_mus_synteny_blocks_50kb_all_fix.bed" \
  	       "results/halsynteny_output/mammalia/bed/can_mus_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/can_mus_synteny_blocks_50kb_target_fix.bed" 

# mus-can

psl_to_bed "results/halsynteny_output/mammalia/mus_can_synteny_50kb.psl" \
           "results/halsynteny_output/mammalia/bed/mus_can_synteny_blocks_50kb_all_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mus_can_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mus_can_synteny_blocks_50kb_target_fix.bed" 

# mac-rat

psl_to_bed "results/halsynteny_output/mammalia/mac_rat_synteny_50kb.psl" \
           "results/halsynteny_output/mammalia/bed/mac_rat_synteny_blocks_50kb_all_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mac_rat_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mac_rat_synteny_blocks_50kb_target_fix.bed" 

# rat-mac

psl_to_bed "results/halsynteny_output/mammalia/rat_mac_synteny_50kb.psl" \
           "results/halsynteny_output/mammalia/bed/rat_mac_synteny_blocks_50kb_all_fix.bed" \
           "results/halsynteny_output/mammalia/bed/rat_mac_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/rat_mac_synteny_blocks_50kb_target_fix.bed" 

# mac-mus

psl_to_bed "results/halsynteny_output/mammalia/mac_mus_synteny_50kb.psl" \
           "results/halsynteny_output/mammalia/bed/mac_mus_synteny_blocks_50kb_all_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mac_mus_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mac_mus_synteny_blocks_50kb_target_fix.bed" 

# mus-mac

psl_to_bed "results/halsynteny_output/mammalia/mus_mac_synteny_50kb.psl" \
           "results/halsynteny_output/mammalia/bed/mus_mac_synteny_blocks_50kb_all_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mus_mac_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mus_mac_synteny_blocks_50kb_target_fix.bed"

# rat-mus

psl_to_bed "results/halsynteny_output/mammalia/rat_mus_synteny_50kb.psl" \
           "results/halsynteny_output/mammalia/bed/rat_mus_synteny_blocks_50kb_all_fix.bed" \
           "results/halsynteny_output/mammalia/bed/rat_mus_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/rat_mus_synteny_blocks_50kb_target_fix.bed" 

# mus-rat

psl_to_bed "results/halsynteny_output/mammalia/mus_rat_synteny_50kb.psl" \
           "results/halsynteny_output/mammalia/bed/mus_rat_synteny_blocks_50kb_all_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mus_rat_synteny_blocks_50kb_query_fix.bed" \
           "results/halsynteny_output/mammalia/bed/mus_rat_synteny_blocks_50kb_target_fix.bed" 


cho "Finished script"
