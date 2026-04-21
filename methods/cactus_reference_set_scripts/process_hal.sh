#!/bin/bash
# Name of the job
#SBATCH --job-name=find_synteny
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
# output file name
#SBATCH --output=slurm_outputs/find_synteny%j.out
#SBATCH --error=slurm_outputs/find_synteny%j.err

# Parse Cactus results and Identify syntenic regions with halSynteny

# Execute your script
echo "Starting the script"

#mkdir -p results/halsynteny_output/brassicaceae
#mkdir -p results/halsynteny_output/drosophila
mkdir -p results/halsynteny_output/mammalia


run_halsynteny() {
  local query=$1              #species abbreviation 
  local target=$2
  local minBlockSize=$3       
  local maxAnchorDistance=$4
  local input_hal=$5
  local output_psl=$6
  local query_genome=$7       # species name stated in seqFile.txt
  local target_genome=$8

  echo "[$(date '+%F %T')] Starting $query-$target"

  halSynteny --minBlockSize $minBlockSize --maxAnchorDistance $maxAnchorDistance $input_hal $output_psl --queryGenome $query_genome --targetGenome $target_genome
  echo "Finished the pair"
}



# Mammalia below




# ------------------------------------------------------------
# Brassicaceae
#-------------------------------------------------------------
#
## ar-card (arabidopsis vs cardamine)
#
#run_halsynteny "ar" \
#               "card" \
#               4000 \
#               9000 \
#               "results/output_alignment_cardamine.hal" \
#               "results/halsynteny_output/brassicaceae/ar_card_synteny_balanced.psl" \
#               "arabidopsis" \
#               "cardamine"
#
## card-ar (cardamine vs arabidopsis)
#
#run_halsynteny "card" \
#               "ar" \
#               4000 \
#               9000 \
#               "results/output_alignment_cardamine.hal" \
#               "results/halsynteny_output/brassicaceae/card_ar_synteny_balanced.psl" \
#               "cardamine" \
#               "arabidopsis"
#
##-----------------------------------------
## Drosophila
##-----------------------------------------
#
## dmel-dsec 
#
#run_halsynteny "dmel" \
#               "dsec" \
#               4645 \
#               9290 \
#               "results/output_alignment_tmux.hal" \
#               "results/halsynteny_output/drosophila/dmel_dsec_synteny_new.psl" \
#               "larger_seq_dmel" \
#               "larger_seq_dsec"
#
## dsec-dmel
#
#run_halsynteny "dsec" \
#               "dmel" \
#               5817 \
#               11635 \
#               "results/output_alignment_tmux.hal" \
#               "results/halsynteny_output/drosophila/dsec_dmel_synteny_new.psl" \
#               "larger_seq_dsec" \
#               "larger_seq_dmel"
#
## dmel-dsim
#
#run_halsynteny "dmel" \
#               "dsim" \
#               4645 \
#               9290 \
#               "results/output_alignment_tmux.hal" \
#               "results/halsynteny_output/drosophila/dmel_dsim_synteny_new.psl" \
#               "larger_seq_dmel" \
#               "larger_seq_dsim"
#
## dsim-dmel
#
#run_halsynteny "dsim" \
#               "dmel" \
#               5067 \
#               10134 \
#               "results/output_alignment_tmux.hal" \
#               "results/halsynteny_output/drosophila/dsim_dmel_synteny_new.psl" \
#               "larger_seq_dsim" \
#               "larger_seq_dmel"               
#
## dsec-dsim
#
#run_halsynteny "dsec" \
#               "dsim" \
#               5817 \
#               11635 \
#               "results/output_alignment_tmux.hal" \
#               "results/halsynteny_output/drosophila/dsec_dsim_synteny_new.psl" \
#               "larger_seq_dsec" \
#               "larger_seq_dsim"
#
## dsim-dsec
#
#run_halsynteny "dsim" \
#               "dsec" \
#               5067 \
#               10134 \
#               "results/output_alignment_tmux.hal" \
#               "results/halsynteny_output/drosophila/dsim_dsec_synteny_new.psl" \
#               "larger_seq_dsim" \
#               "larger_seq_dsec"


#--------------------------------------------------
# Mammalia
#--------------------------------------------------

# can-mac

run_halsynteny "can" \
               "mac" \
               50000 \
               50000 \
               "results/output_alignment_sexspecies.hal" \
               "results/halsynteny_output/mammalia/can_mac_synteny_50kb.psl" \
               "Canis_lupus" \
               "Macaca_fascicularis"

# mac-can

run_halsynteny "mac" \
               "can" \
               50000 \
               50000 \
               "results/output_alignment_sexspecies.hal" \
               "results/halsynteny_output/mammalia/mac_can_synteny_50kb.psl" \
               "Macaca_fascicularis" \
               "Canis_lupus"

# can-rat

run_halsynteny "can" \
               "rat" \
               50000 \
               50000 \
               "results/output_alignment_sexspecies.hal" \
               "results/halsynteny_output/mammalia/can_rat_synteny_50kb.psl" \
               "Canis_lupus" \
               "Rattus_norvegicus"

# rat-can

run_halsynteny "rat" \
               "can" \
               50000 \
               50000 \
               "results/output_alignment_sexspecies.hal" \
               "results/halsynteny_output/mammalia/rat_can_synteny_50kb.psl" \
               "Rattus_norvegicus" \
               "Canis_lupus"

# can-mus 

run_halsynteny "can" \
               "mus" \
               50000 \
               50000 \
               "results/output_alignment_sexspecies.hal" \
               "results/halsynteny_output/mammalia/can_mus_synteny_50kb.psl" \
               "Canis_lupus" \
               "Mus_musculus"

# mus-can

run_halsynteny "mus" \
               "can" \
               50000 \
               50000 \
               "results/output_alignment_sexspecies.hal" \
               "results/halsynteny_output/mammalia/mus_can_synteny_50kb.psl" \
               "Mus_musculus" \
               "Canis_lupus"

# mac-rat 

run_halsynteny "mac" \
               "rat" \
               50000 \
               50000 \
               "results/output_alignment_sexspecies.hal" \
               "results/halsynteny_output/mammalia/mac_rat_synteny_50kb.psl" \
               "Macaca_fascicularis" \
               "Rattus_norvegicus"

# rat-mac

run_halsynteny "rat" \
               "mac" \
               50000 \
               50000 \
               "results/output_alignment_sexspecies.hal" \
               "results/halsynteny_output/mammalia/rat_mac_synteny_50kb.psl" \
               "Rattus_norvegicus" \
               "Macaca_fascicularis"

# mac-mus

run_halsynteny "mac" \
               "mus" \
               50000 \
               50000 \
               "results/output_alignment_sexspecies.hal" \
               "results/halsynteny_output/mammalia/mac_mus_synteny_50kb.psl" \
               "Macaca_fascicularis" \
               "Mus_musculus"

# mus-mac

run_halsynteny "mus" \
               "mac" \
               50000 \
               50000 \
               "results/output_alignment_sexspecies.hal" \
               "results/halsynteny_output/mammalia/mus_mac_synteny_50kb.psl" \
               "Mus_musculus" \
               "Macaca_fascicularis"

# rat-mus

run_halsynteny "rat" \
               "mus" \
               50000 \
               50000 \
               "results/output_alignment_sexspecies.hal" \
               "results/halsynteny_output/mammalia/rat_mus_synteny_50kb.psl" \
               "Rattus_norvegicus" \
               "Mus_musculus"

# mus-rat

run_halsynteny "mus" \
               "rat" \
               50000 \
               50000 \
               "results/output_alignment_sexspecies.hal" \
               "results/halsynteny_output/mammalia/mus_rat_synteny_50kb.psl" \
               "Mus_musculus" \
               "Rattus_norvegicus"



echo "Finished script."

