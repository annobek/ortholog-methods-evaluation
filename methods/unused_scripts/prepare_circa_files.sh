#!/bin/bash
#SBATCH --job-name=prepare_circa_files
#SBATCH --output=slurm_outputs/prepare_circa_files_%j.out

# Prepare Circa input files 

echo "Starting the script"

set -euo pipefail

THRESHOLD=${1:-1000000}  # default block length filter (1 Mb)
OUTDIR="results/circos_plots/circa"
FAI_DIR="material/sex_experiment/nucleotids"
BED_DIR="results/halsynteny_output/bed"

mkdir -p "$OUTDIR" slurm_outputs


echo "Starting Circa preparation with threshold: ${THRESHOLD}"

# Generate sizes and big-chromosome (filtered) files from .fai
generate_sizes_files() {
  local fai_input=$1
  local prefix=$2
  local outdir=$3

  local sizes_out="${outdir}/${prefix}.sizes.tsv"
  local big_out="${outdir}/${prefix}_big_chr.sizes.tsv"

  echo "Processing ${fai_input} -> ${sizes_out}"

  # Create .sizes.tsv from FAI
  awk 'BEGIN{OFS="\t"} {print $1,$2}' "$fai_input" > "$sizes_out"

  # Keep only major (NC_) chromosomes
  awk '$1 ~ /^NC_/' "$sizes_out" > "$big_out"
}

# Combine two species .sizes files into one setup.tsv for Circa
combine_setup_file() {
  local query_big=$1
  local target_big=$2
  local query_prefix=$3
  local target_prefix=$4
  local setup_out=$5

  echo "Creating setup file: ${setup_out}"

  {
    echo -e "chromosome\tsize\tflip"
    awk -v pfx="$query_prefix" 'BEGIN{OFS="\t"} {print pfx"_"$1,$2,"false"}' "$query_big"
    awk -v pfx="$target_prefix" 'BEGIN{OFS="\t"} {print pfx"_"$1,$2,"true"}'  "$target_big"
  } > "$setup_out"
}

# Prepare synteny block files -> filtered ribbons.tsv
process_blocks() {
  local bed_input=$1
  local prefix=$2
  local outdir=$3
  local threshold=$4
  local query_name=$5
  local target_name=$6

  local clean_tsv="${outdir}/${prefix}_blocks.tsv"
  local filtered_tsv="${outdir}/${prefix}_blocks_filtered.tsv"
  local ribbons_out="${outdir}/${prefix}_ribbons.tsv"

  echo "Processing blocks: ${bed_input}"

  # BED -> clean TSV (remove block ID, keep query and target coords)
  awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$5,$6,$7}' "$bed_input" > "$clean_tsv"

  # Filter by size threshold on both sides
  awk -v t=$threshold '($3-$2)>=t && ($6-$5)>=t' "$clean_tsv" > "$filtered_tsv"

  # Create ribbons file for Circa
  awk -v q="$query_name" -v t="$target_name" '
    BEGIN {OFS="\t"; print "chr1","start1","end1","chr2","start2","end2"}
    {print q"_"$1,$2,$3,t"_"$4,$5,$6}
  ' "$filtered_tsv" > "$ribbons_out"

}


# ---- Sizes files ----
generate_sizes_files "${FAI_DIR}/canis_lupus.fna.fai"         "canis"  "$OUTDIR"
generate_sizes_files "${FAI_DIR}/macaca_fascicularis.fna.fai" "macaca" "$OUTDIR"
generate_sizes_files "${FAI_DIR}/rattus_norvegicus.fna.fai"   "rat"    "$OUTDIR"
generate_sizes_files "${FAI_DIR}/mus_musculus.fna.fai"        "mouse"  "$OUTDIR"

# ---- Setup files for pairwise comparisons ----
combine_setup_file "$OUTDIR/canis_big_chr.sizes.tsv" "$OUTDIR/macaca_big_chr.sizes.tsv" "can" "mac" "$OUTDIR/can_mac_setup.tsv"
combine_setup_file "$OUTDIR/canis_big_chr.sizes.tsv" "$OUTDIR/rat_big_chr.sizes.tsv"    "can" "rat" "$OUTDIR/can_rat_setup.tsv"
combine_setup_file "$OUTDIR/canis_big_chr.sizes.tsv" "$OUTDIR/mouse_big_chr.sizes.tsv"  "can" "mus" "$OUTDIR/can_mus_setup.tsv"
combine_setup_file "$OUTDIR/macaca_big_chr.sizes.tsv" "$OUTDIR/rat_big_chr.sizes.tsv"   "mac" "rat" "$OUTDIR/mac_rat_setup.tsv"
combine_setup_file "$OUTDIR/macaca_big_chr.sizes.tsv" "$OUTDIR/mouse_big_chr.sizes.tsv" "mac" "mus" "$OUTDIR/mac_mus_setup.tsv"
combine_setup_file "$OUTDIR/rat_big_chr.sizes.tsv" "$OUTDIR/mouse_big_chr.sizes.tsv"    "rat" "mus" "$OUTDIR/rat_mus_setup.tsv"

# ---- Process block files (create filtered + ribbons) ----
process_blocks "$BED_DIR/can_mac_synteny_blocks_50kb_all_fix.bed" "can_mac" "$OUTDIR" "$THRESHOLD" "can" "mac"
#process_blocks "$BED_DIR/mac_can_synteny_blocks_50kb_all_fix.bed" "mac_can" "$OUTDIR" "$THRESHOLD"
process_blocks "$BED_DIR/can_rat_synteny_blocks_50kb_all_fix.bed" "can_rat" "$OUTDIR" "$THRESHOLD" "can" "rat"
#process_blocks "$BED_DIR/rat_can_synteny_blocks_50kb_all_fix.bed" "rat_can" "$OUTDIR" "$THRESHOLD"
process_blocks "$BED_DIR/can_mus_synteny_blocks_50kb_all_fix.bed" "can_mus" "$OUTDIR" "$THRESHOLD" "can" "mus"
#process_blocks "$BED_DIR/mus_can_synteny_blocks_50kb_all_fix.bed" "mus_can" "$OUTDIR" "$THRESHOLD"
process_blocks "$BED_DIR/mac_rat_synteny_blocks_50kb_all_fix.bed" "mac_rat" "$OUTDIR" "$THRESHOLD" "mac" "rat"
#process_blocks "$BED_DIR/rat_mac_synteny_blocks_50kb_all_fix.bed" "rat_mac" "$OUTDIR" "$THRESHOLD"
process_blocks "$BED_DIR/mac_mus_synteny_blocks_50kb_all_fix.bed" "mac_mus" "$OUTDIR" "$THRESHOLD" "mac" "mus"
#process_blocks "$BED_DIR/mus_mac_synteny_blocks_50kb_all_fix.bed" "mus_mac" "$OUTDIR" "$THRESHOLD"
process_blocks "$BED_DIR/rat_mus_synteny_blocks_50kb_all_fix.bed" "rat_mus" "$OUTDIR" "$THRESHOLD" "rat" "mus"
#process_blocks "$BED_DIR/mus_rat_synteny_blocks_50kb_all_fix.bed" "mus_rat" "$OUTDIR" "$THRESHOLD"

echo "Finished script."
