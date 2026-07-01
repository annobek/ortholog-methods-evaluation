#!/bin/bash
# Name of the job
#SBATCH --job-name=block_length
#SBATCH --output=slurm_outputs/block_length%j.out

# Check the length of the synteny blocks and compare forvard with reverse

echo "Starting the script"

mkdir -p results/data_check/block_length/

block_length() {
  local input=$1
  local output_length=$2
  local species1=$3
  local species2=$4   

  echo "Checking $species1-$species2"
  echo "---------------------"

  # Check for header and skip if present
  if head -n1 $input | grep -q "Query_Scaffold"; then
    skip_cmd="tail -n +2"
  else
    skip_cmd="cat"
  fi

  # Process the file and save lengths directly
  echo -e "Region_ID\tQuery_Scaffold\tQuery_Start\tQuery_End\tTarget_Scaffold\tTarget_Start\tTarget_End\tquery_length\ttarget_length" > $output_length

  $skip_cmd "$input" | \
    awk 'BEGIN{OFS="\t"} {
      query_len = $4 - $3
      target_len = $7 - $6
      print $1, $2, $3, $4, $5, $6, $7, query_len, target_len
    }' >> "$output_length"


  # Calculate combined summary stats (use last two columns)
 awk 'NR > 1 {
    # Extract query and target block lengths (last two columns)
    query_length  = $(NF-1)
    target_length = $NF

    # Store query stats
    all_lengths_query[++num_query_blocks] = query_length
    sum_query += query_length
    if (min_query == "" || query_length < min_query) min_query = query_length
    if (query_length > max_query) max_query = query_length

    # Store target stats
    all_target_lengths[++num_target_blocks] = target_length
    sum_target += target_length
    if (min_target == "" || target_length < min_target) min_target = target_length
    if (target_length > max_target) max_target = target_length
  }

  END {
    # Compute medians (requires GNU Awk for asort)
    asort(all_lengths_query)
    median_query = (num_query_blocks % 2 ? all_lengths_query[(num_query_blocks + 1) / 2] : (all_lengths_query[num_query_blocks / 2] + all_lengths_query[num_query_blocks / 2 + 1]) / 2)

    asort(all_target_lengths)
    median_target = (num_target_blocks % 2 ? all_target_lengths[(num_target_blocks + 1) / 2] : (all_target_lengths[num_target_blocks / 2] + all_target_lengths[num_target_blocks / 2 + 1]) / 2)

    # Print formatted summary
    printf "Query blocks:  num_blocks=%d\ttotal_coverage=%d\tmin_length=%d\tmax_length=%d\tmean_length=%.1f\tmedian_length=%.1f\n",
          num_query_blocks, sum_query, min_query, max_query,
          (num_query_blocks ? sum_query / num_query_blocks : 0), median_query

    printf "Target blocks: num_blocks=%d\ttotal_coverage=%d\tmin_length=%d\tmax_length=%d\tmean_length=%.1f\tmedian_length=%.1f\n",
          num_target_blocks, sum_target, min_target, max_target,
          (num_target_blocks ? sum_target / num_target_blocks : 0), median_target
  }' "$output_length"


  # Histogram: use last two columns (qlen = $(NF-1), tlen = $NF)
  #{
  #  echo -e "Bin\tQuery_Count\tTarget_Count"
  #  awk -v B="$BIN" '
  #    { q=int(($(NF-1))/B)*B; t=int(($NF)/B)*B; Q[q]++; T[t]++; seen[q]=1; seen[t]=1 }
  #    END{
  #      n=asorti(seen, keys);
  #      for(i=1;i<=n;i++){
  #        k=keys[i];
  #        printf "%d\t%d\t%d\n", k, (Q[k]+0), (T[k]+0)
  #      }
  #    }' "$output_length" | sort -k1,1n
  #}

  echo ""
  echo "======================================================="
  echo ""


}

compare_pairs() {
  local forward_file=$1
  local reverse_file=$2
  local species1=$3       # for example can
  local species2=$4       # for example mac

  echo ""
  echo "Comparing $species1-$species2 vs $species2-$species1" 
  echo "------------------------------"

  # Count number of blocks
  forward_count=$(( $(wc -l < "$forward_file") - 1 ))
  reverse_count=$(( $(wc -l < "$reverse_file") - 1 ))

  
  echo "Block counts:"
  echo "  ${species1}-${species2}: $forward_count blocks"
  echo "  ${species2}-${species1}: $reverse_count blocks"
  echo "  Difference: $((forward_count - reverse_count)) blocks"
  
  # Compare total coverage
  echo ""
  echo "Total coverage:"
  forward_query_total=$(awk '{sum+=$8} END {printf "%.0f", sum}' "$forward_file")
  forward_target_total=$(awk '{sum+=$9} END {printf "%.0f", sum}' "$forward_file")
  reverse_query_total=$(awk '{sum+=$8} END {printf "%.0f", sum}' "$reverse_file")
  reverse_target_total=$(awk '{sum+=$9} END {printf "%.0f", sum}' "$reverse_file")
  
  echo "  ${species1}-${species2}: ${species1}=$(( forward_query_total / 1000000 ))Mb, ${species2}=$(( forward_target_total / 1000000 ))Mb"
  echo "  ${species2}-${species1}: ${species2}=$(( reverse_query_total / 1000000 ))Mb, ${species1}=$(( reverse_target_total / 1000000 ))Mb"
  
  # Check if coverage is reciprocal (should be similar)
  echo ""
  echo "Coverage consistency check:"
  echo "  ${species1} in forward vs reverse: $(( forward_query_total / 1000000 ))Mb vs $(( reverse_target_total / 1000000 ))Mb"
  echo "  ${species2} in forward vs reverse: $(( forward_target_total / 1000000 ))Mb vs $(( reverse_query_total / 1000000 ))Mb"

  echo ""
  echo "----------------------------------------"

  # Create a summary table
  echo ""
  echo "SUMMARY TABLE"
  echo "----------------------------------"
  echo -e "Pair\tForward_Blocks\tReverse_Blocks\tDifference"
  diff=$((forward_count - reverse_count))
  echo -e "${species1}-${species2}\t${forward_count}\t${reverse_count}\t${diff}"
}



# Mammalia below


# Statistics:

#------------------------------------------------------
# Brassicaceae
#------------------------------------------------------

# ar-card
#
#block_length "results/halsynteny_output/brassicaceae/bed/ar_card_synteny_all_balanced.bed" \
#  "results/data_check/block_length/ar_card_block_length_balanced.tsv" \
#  "ar" \
#  "card"
#
## card-ar
#
#block_length "results/halsynteny_output/brassicaceae/bed/card_ar_synteny_all_balanced.bed" \
#  "results/data_check/block_length/card_ar_block_length_balanced.tsv" \
#  "card" \
#  "ar"
#
#
##--------------------------------------------------------------
## Drosophila
##--------------------------------------------------------------
#
## dmel-dsec
#
#block_length "results/halsynteny_output/drosophila/bed/dmel_dsec_synteny_all_new.bed" \
#  "results/data_check/block_length/dmel_dsec_block_length_new.tsv" \
#  "dmel" \
#  "dsec"
#
## dsec-dmel
#
#block_length "results/halsynteny_output/drosophila/bed/dsec_dmel_synteny_all_new.bed" \
#  "results/data_check/block_length/dsec_dmel_block_length_new.tsv" \
#  "dsec" \
#  "dmel"
#
## dmel-dsim
#
#block_length "results/halsynteny_output/drosophila/bed/dmel_dsim_synteny_all_new.bed" \
#  "results/data_check/block_length/dmel_dsim_block_length_new.tsv" \
#  "dmel" \
#  "dsim"
#
## dsim-dmel
#
#block_length "results/halsynteny_output/drosophila/bed/dsim_dmel_synteny_all_new.bed" \
#  "results/data_check/block_length/dsim_dmel_block_length_new.tsv" \
#  "dsim" \
#  "dmel"
#
## dsim-dsec
#
#block_length "results/halsynteny_output/drosophila/bed/dsim_dsec_synteny_all_new.bed" \
#  "results/data_check/block_length/dsim_dsec_block_length_new.tsv" \
#  "dsim" \
#  "dsec"
#
## dsec-dsim
#
#block_length "results/halsynteny_output/drosophila/bed/dsec_dsim_synteny_all_new.bed" \
#  "results/data_check/block_length/dsec_dsim_block_length_new.tsv" \
#  "dsec" \
#  "dsim"
#
#---------------------------------------------------------------
# Mammalia
#---------------------------------------------------------------

# can-mac

block_length "results/halsynteny_output/mammalia/bed/can_mac_synteny_blocks_50kb_all_fix.bed" \
  "results/data_check/block_length/can_mac_block_length_50kb_fix.tsv" \
  "can" \
  "mac"

# mac-can

block_length "results/halsynteny_output/mammalia/bed/mac_can_synteny_blocks_50kb_all_fix.bed" \
  "results/data_check/block_length/mac_can_block_length_50kb_fix.tsv" \
  "mac" \
  "can"

# can-rat

block_length "results/halsynteny_output/mammalia/bed/can_rat_synteny_blocks_50kb_all_fix.bed" \
  "results/data_check/block_length/can_rat_block_length_50kb_fix.tsv" \
  "can" \
  "rat"

# rat-can

block_length "results/halsynteny_output/mammalia/bed/rat_can_synteny_blocks_50kb_all_fix.bed" \
  "results/data_check/block_length/rat_can_block_length_50kb_fix.tsv" \
  "rat" \
  "can"

# can-mus

block_length "results/halsynteny_output/mammalia/bed/can_mus_synteny_blocks_50kb_all_fix.bed" \
  "results/data_check/block_length/can_mus_block_length_50kb_fix.tsv" \
  "can" \
  "mus"

# mus-can

block_length "results/halsynteny_output/mammalia/bed/mus_can_synteny_blocks_50kb_all_fix.bed" \
  "results/data_check/block_length/mus_can_block_length_50kb_fix.tsv" \
  "mus" \
  "can"

# mac-rat

block_length "results/halsynteny_output/mammalia/bed/mac_rat_synteny_blocks_50kb_all_fix.bed" \
  "results/data_check/block_length/mac_rat_block_length_50kb_fix.tsv" \
  "mac" \
  "rat"

# rat-mac

block_length "results/halsynteny_output/mammalia/bed/rat_mac_synteny_blocks_50kb_all_fix.bed" \
  "results/data_check/block_length/rat_mac_block_length_50kb_fix.tsv" \
  "rat" \
  "mac"

# mac-mus

block_length "results/halsynteny_output/mammalia/bed/mac_mus_synteny_blocks_50kb_all_fix.bed" \
  "results/data_check/block_length/mac_mus_block_length_50kb_fix.tsv" \
  "mac" \
  "mus"

# mus-mac

block_length "results/halsynteny_output/mammalia/bed/mus_mac_synteny_blocks_50kb_all_fix.bed" \
  "results/data_check/block_length/mus_mac_block_length_50kb_fix.tsv" \
  "mus" \
  "mac"

# rat-mus

block_length "results/halsynteny_output/mammalia/bed/rat_mus_synteny_blocks_50kb_all_fix.bed" \
  "results/data_check/block_length/rat_mus_block_length_50kb_fix.tsv" \
  "rat" \
  "mus"

# mus-rat

block_length "results/halsynteny_output/mammalia/bed/mus_rat_synteny_blocks_50kb_all_fix.bed" \
  "results/data_check/block_length/mus_rat_block_length_50kb_fix.tsv" \
  "mus" \
  "rat"


# Mammalia below


# Comparison forward vs reverse:

##--------------------
## Plants
##--------------------
#
## ar-card vs card-ar
#
#compare_pairs "results/data_check/block_length/ar_card_block_length_balanced.tsv" \
#              "results/data_check/block_length/card_ar_block_length_balanced.tsv" \
#              "ar" \
#              "card"
#
#
##-----------------
## Drosophila
##-----------------
#
## dmel-dsec vs dsec-dmel
#
#compare_pairs "results/data_check/block_length/dmel_dsec_block_length_new.tsv" \
#              "results/data_check/block_length/dsec_dmel_block_length_new.tsv" \
#              "dmel" \
#              "dsec"
#
## dmel-dsim vs dsim-dmel
#
#compare_pairs "results/data_check/block_length/dmel_dsim_block_length_new.tsv" \
#              "results/data_check/block_length/dsim_dmel_block_length_new.tsv" \
#              "dmel" \
#              "dsim"
#
## dsec-dsim vs dsim-dsec
#
#compare_pairs "results/data_check/block_length/dsec_dsim_block_length_new.tsv" \
#              "results/data_check/block_length/dsim_dsec_block_length_new.tsv" \
#              "dsec" \
#              "dsim"
#

#------------------
# Mammalia
#------------------

# can-mac vs mac-can

compare_pairs "results/data_check/block_length/can_mac_block_length_50kb_fix.tsv" \
              "results/data_check/block_length/mac_can_block_length_50kb_fix.tsv" \
              "can" \
              "mac"

# can-rat vs rat-can

compare_pairs "results/data_check/block_length/can_rat_block_length_50kb_fix.tsv" \
              "results/data_check/block_length/rat_can_block_length_50kb_fix.tsv" \
              "can" \
              "rat"

# can-mus vs mus-can

compare_pairs "results/data_check/block_length/can_mus_block_length_50kb_fix.tsv" \
              "results/data_check/block_length/mus_can_block_length_50kb_fix.tsv" \
              "can" \
              "mus"

# mac-rat vs rat-mac

compare_pairs "results/data_check/block_length/mac_rat_block_length_50kb_fix.tsv" \
              "results/data_check/block_length/rat_mac_block_length_50kb_fix.tsv" \
              "mac" \
              "rat"

# mac-mus vs mus-mac

compare_pairs "results/data_check/block_length/mac_mus_block_length_50kb_fix.tsv" \
              "results/data_check/block_length/mus_mac_block_length_50kb_fix.tsv" \
              "mac" \
              "mus"              

# rat-mus vs mus-rat

compare_pairs "results/data_check/block_length/rat_mus_block_length_50kb_fix.tsv" \
              "results/data_check/block_length/mus_rat_block_length_50kb_fix.tsv" \
              "rat" \
              "mus"

echo "Finished script"
