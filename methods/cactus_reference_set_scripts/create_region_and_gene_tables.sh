#!/bin/bash
# Name of the job
#SBATCH --job-name=create_region_and_gene_tables
#SBATCH --output=slurm_outputs/create_region_and_gene_tables%j.out

# Create region and gene tables from bedtools intersect output

echo "Starting the script"
echo "Processing bedtools output to create final tables..."

create_tables() {
  local query_file=$1
  local target_file=$2
  local synteny_blocks_query=$3
  local synteny_blocks_target=$4
  local species1=$5
  local species2=$6
  local species_for_table1=$7 # species full name
  local species_for_table2=$8
  local species_dir=$9 # mammalia, plants or drosophila

  local outdir="results/genes_in_blocks/${species_dir}/region_gene_tables"
  local output_prefix="${species1}_${species2}_synteny_tables"

  mkdir -p "$outdir"

  echo "Creating region table..."


  # Build the region table from the raw blocks (all blocks)
  awk '{print $4"\t"$1"\t"$2"\t"$3}' "$synteny_blocks_query" \
  | awk '{n=$1; gsub(/block_/, "", n); print n"\t"$0}' \
  | sort -n -k1,1 \
  | cut -f2- \
  | sort -u -k1,1 > "${outdir}/tmp_query_blocks.txt"

  awk '{print $4"\t"$1"\t"$2"\t"$3}' "$synteny_blocks_target" \
  | awk '{n=$1; gsub(/block_/, "", n); print n"\t"$0}' \
  | sort -n -k1,1 \
  | cut -f2- \
  | sort -u -k1,1 > "${outdir}/tmp_target_blocks.txt"

  # Header and join (inner join is fine here because the raw block files
  # are symmetrical; they both contain all block IDs)
  echo -e "Region_ID\tQuery_Scaffold\tQuery_Start\tQuery_End\tTarget_Scaffold\tTarget_Start\tTarget_End" \
  > "${outdir}/${output_prefix}_regions.tsv"

  join -t $'\t' "${outdir}/tmp_query_blocks.txt" "${outdir}/tmp_target_blocks.txt" \
  | awk 'BEGIN{OFS="\t"} {gsub(/block_/,"",$1); print $1,$2,$3,$4,$5,$6,$7}' \
  | sort -n -k1,1 >> "${outdir}/${output_prefix}_regions.tsv"
  echo "Region table created: ${outdir}/${output_prefix}_regions.tsv"

  echo "Creating genes table..."
  echo "Region_ID	Species	Gene-ID	Start	End	percentage" > "${outdir}/${output_prefix}_genes.tsv"

  # Process query species genes
  awk -v species="$species_for_table1" 'BEGIN{OFS="\t"} {
      # Extract region ID
      region = $4
      gsub(/block_/, "", region)
      
      # Calculate percentage (overlap / gene_length * 100)
      gene_length = $7 - $6
      if (gene_length > 0) {
          percentage = ($11 / gene_length) * 100
      } else {
          percentage = 0
      }
      
      printf "%s\t%s\t%s\t%d\t%d\t%.2f\n", region, species, $8, $6, $7, percentage
   }' "$query_file" >> "${outdir}/${output_prefix}_genes.tsv"

  # Process target species genes  
  awk -v species="$species_for_table2" 'BEGIN{OFS="\t"} {
      # Extract region ID
      region = $4
      gsub(/block_/, "", region)
      
      # Calculate percentage
      gene_length = $7 - $6
      if (gene_length > 0) {
          percentage = ($11 / gene_length) * 100
      } else {
          percentage = 0
      }
      
      printf "%s\t%s\t%s\t%d\t%d\t%.2f\n", region, species, $8, $6, $7, percentage
  }' "$target_file" >> "${outdir}/${output_prefix}_genes.tsv"

  # Sort genes table by Region_ID and Species
  sort -t$'\t' -k1,1n -k2,2 -k3,3 "${outdir}/${output_prefix}_genes.tsv" > "${outdir}/tmp_sorted.tsv"
  head -1 "${outdir}/${output_prefix}_genes.tsv" > "${outdir}/tmp_header.tsv"
  tail -n +2 "${outdir}/tmp_sorted.tsv" >> "${outdir}/tmp_header.tsv"
  mv "${outdir}/tmp_header.tsv" "${outdir}/${output_prefix}_genes.tsv"

  echo "Genes table created: ${outdir}/${output_prefix}_genes.tsv"

  echo ""
  echo "========================================="
  echo "SUMMARY STATISTICS"
  echo "========================================="
  echo ""
  echo "For $species1-$species2:"
  echo ""
  # Count regions
  total_regions=$(tail -n +2 "${outdir}/${output_prefix}_regions.tsv" | wc -l)
  echo "Total syntenic regions: $total_regions"

  # Count genes
  query_genes=$(awk -F'\t' -v s="$species_for_table1" 'NR>1 && $2==s' "${outdir}/${output_prefix}_genes.tsv" | wc -l)
  target_genes=$(awk -F'\t' -v s="$species_for_table2" 'NR>1 && $2==s' "${outdir}/${output_prefix}_genes.tsv" | wc -l)

  echo "Genes in syntenic regions:"
  echo "  $species1: $query_genes"
  echo "  $species2: $target_genes"

  # Regions with genes in query and target 
  regions_with_genes_query=$(awk -v s="$species_for_table1" 'NR>1 && index($0, s) > 0 {print $1}' "${outdir}/${output_prefix}_genes.tsv" | sort -u | wc -l)
  regions_with_genes_target=$(awk -v s="$species_for_table2" 'NR>1 && index($0, s) > 0 {print $1}' "${outdir}/${output_prefix}_genes.tsv" | sort -u | wc -l)

  # Regions that have at least one gene from both species
  regions_with_both=$(awk 'NR>1 {print $1, $2}' "${outdir}/${output_prefix}_genes.tsv" | \
    sort -u | awk '{count[$1]++} END {for (r in count) if (count[r]>1) c++; print c}')

  # Regions that have only genes in query or only genes in target
  regions_query_only=$((regions_with_genes_query - regions_with_both))
  regions_target_only=$((regions_with_genes_target - regions_with_both))

  # Total number of asymmetrical regions (one side has genes, other not and reverse)
  asymmetrical_regions=$((regions_query_only + regions_target_only))

  # Regions with no genes for both
  regions_no_genes=$((total_regions - regions_with_genes_query - regions_with_genes_target + regions_with_both))
  if (( regions_no_genes < 0 )); then regions_no_genes=0; fi

  echo "Regions with genes (query): $regions_with_genes_query"
  echo "Regions with genes (target): $regions_with_genes_target"
  echo "Regions with genes from both: $regions_with_both"
  echo "Total number of asymmetrical regions (one-sided): $asymmetrical_regions"
  echo "Regions with no genes at all: $regions_no_genes"

  # Clean up temp files
  rm -f "${outdir}"/tmp_*.txt "${outdir}"/tmp_*.tsv

  echo ""
  echo "Processing complete!"
  echo "Output files:"
  echo "${outdir}/${output_prefix}_regions.tsv"
  echo "${outdir}/${output_prefix}_genes.tsv"

}


# Mammalia below



## Plants
##-------
#
## ar-card
#create_tables "results/genes_in_blocks/brassicaceae/ar_card_genes_in_blocks_query.bed" \
#              "results/genes_in_blocks/brassicaceae/ar_card_genes_in_blocks_target.bed" \
#              "results/halsynteny_output/brassicaceae/bed/ar_card_synteny_query_balanced.bed" \
#              "results/halsynteny_output/brassicaceae/bed/ar_card_synteny_target_balanced.bed" \
#              "ar" \
#              "card" \
#              "Arabidopsis thaliana" \
#              "Cardamine hirsuta" \
#              "brassicaceae"
#
## card-ar
#create_tables "results/genes_in_blocks/brassicaceae/card_ar_genes_in_blocks_query.bed" \
#              "results/genes_in_blocks/brassicaceae/card_ar_genes_in_blocks_target.bed" \
#              "results/halsynteny_output/brassicaceae/bed/card_ar_synteny_query_balanced.bed" \
#              "results/halsynteny_output/brassicaceae/bed/card_ar_synteny_target_balanced.bed" \
#              "card" \
#              "ar" \
#              "Cardamine hirsuta" \
#              "Arabidopsis thaliana" \
#              "brassicaceae"
#
## Drosophila
##-----------
#
## dmel-dsec
#create_tables "results/genes_in_blocks/drosophila/dmel_dsec_genes_in_blocks_query.bed" \
#              "results/genes_in_blocks/drosophila/dmel_dsec_genes_in_blocks_target.bed" \
#              "results/halsynteny_output/drosophila/bed/dmel_dsec_synteny_query_new.bed" \
#              "results/halsynteny_output/drosophila/bed/dmel_dsec_synteny_target_new.bed" \
#              "dmel" \
#              "dsec" \
#              "Drosophila melanogaster" \
#              "Drosophila sechellia" \
#              "drosophila"
#
## dsec-dmel
#create_tables "results/genes_in_blocks/drosophila/dsec_dmel_genes_in_blocks_query.bed" \
#              "results/genes_in_blocks/drosophila/dsec_dmel_genes_in_blocks_target.bed" \
#              "results/halsynteny_output/drosophila/bed/dsec_dmel_synteny_query_new.bed" \
#              "results/halsynteny_output/drosophila/bed/dsec_dmel_synteny_target_new.bed" \
#              "dsec" \
#              "dmel" \
#              "Drosophila sechellia" \
#              "Drosophila melanogaster" \
#              "drosophila"
#
#
## dmel-dsim
#create_tables "results/genes_in_blocks/drosophila/dmel_dsim_genes_in_blocks_query.bed" \
#              "results/genes_in_blocks/drosophila/dmel_dsim_genes_in_blocks_target.bed" \
#              "results/halsynteny_output/drosophila/bed/dmel_dsim_synteny_query_new.bed" \
#              "results/halsynteny_output/drosophila/bed/dmel_dsim_synteny_target_new.bed" \
#              "dmel" \
#              "dsim" \
#              "Drosophila melanogaster" \
#              "Drosophila simulans" \
#              "drosophila"
#
## dsim-dmel
#create_tables "results/genes_in_blocks/drosophila/dsim_dmel_genes_in_blocks_query.bed" \
#              "results/genes_in_blocks/drosophila/dsim_dmel_genes_in_blocks_target.bed" \
#              "results/halsynteny_output/drosophila/bed/dsim_dmel_synteny_query_new.bed" \
#              "results/halsynteny_output/drosophila/bed/dsim_dmel_synteny_target_new.bed" \
#              "dsim" \
#              "dmel" \
#              "Drosophila simulans" \
#              "Drosophila melanogaster" \
#              "drosophila"
#
## dsec-dsim
#create_tables "results/genes_in_blocks/drosophila/dsec_dsim_genes_in_blocks_query.bed" \
#              "results/genes_in_blocks/drosophila/dsec_dsim_genes_in_blocks_target.bed" \
#              "results/halsynteny_output/drosophila/bed/dsec_dsim_synteny_query_new.bed" \
#              "results/halsynteny_output/drosophila/bed/dsec_dsim_synteny_target_new.bed" \
#              "dsec" \
#              "dsim" \
#              "Drosophila sechellia" \
#              "Drosophila simulans" \
#              "drosophila"
#
## dsim-dsec
#create_tables "results/genes_in_blocks/drosophila/dsim_dsec_genes_in_blocks_query.bed" \
#              "results/genes_in_blocks/drosophila/dsim_dsec_genes_in_blocks_target.bed" \
#              "results/halsynteny_output/drosophila/bed/dsim_dsec_synteny_query_new.bed" \
#              "results/halsynteny_output/drosophila/bed/dsim_dsec_synteny_target_new.bed" \
#              "dsim" \
#              "dsec" \
#              "Drosophila simulans" \
#              "Drosophila sechellia" \
#              "drosophila"
#
#
# Mammals
#--------

# can-mac
create_tables "results/genes_in_blocks/mammalia/can_mac_genes_in_blocks_50kb_query_fix.bed" \
              "results/genes_in_blocks/mammalia/can_mac_genes_in_blocks_50kb_target_fix.bed" \
              "results/halsynteny_output/bed/mammalia/can_mac_synteny_blocks_50kb_query_fix.bed" \
              "results/halsynteny_output/bed/mammalia/can_mac_synteny_blocks_50kb_target_fix.bed" \
              "can" \
              "mac" \
              "Canis lupus familiaris" \
              "Macaca fascicularis" \
              "mammalia"

# mac-can
create_tables "results/genes_in_blocks/mammalia/mac_can_genes_in_blocks_50kb_query_fix.bed" \
              "results/genes_in_blocks/mammalia/mac_can_genes_in_blocks_50kb_target_fix.bed" \
              "results/halsynteny_output/mammalia/bed/mac_can_synteny_blocks_50kb_query_fix.bed" \
              "results/halsynteny_output/mammalia/bed/mac_can_synteny_blocks_50kb_target_fix.bed" \
              "mac" \
              "can" \
              "Macaca fascicularis" \
              "Canis lupus familiaris" \
              "mammalia"

# can-rat
create_tables "results/genes_in_blocks/mammalia/can_rat_genes_in_blocks_50kb_query_fix.bed" \
              "results/genes_in_blocks/mammalia/can_rat_genes_in_blocks_50kb_target_fix.bed" \
              "results/halsynteny_output/mammalia/bed/can_rat_synteny_blocks_50kb_query_fix.bed" \
              "results/halsynteny_output/mammalia/bed/can_rat_synteny_blocks_50kb_target_fix.bed" \
              "can" \
              "rat" \
              "Canis lupus familiaris" \
              "Rattus norvegicus" \
              "mammalia"

# rat-can
create_tables "results/genes_in_blocks/mammalia/rat_can_genes_in_blocks_50kb_query_fix.bed" \
              "results/genes_in_blocks/mammalia/rat_can_genes_in_blocks_50kb_target_fix.bed" \
              "results/halsynteny_output/mammalia/bed/rat_can_synteny_blocks_50kb_query_fix.bed" \
              "results/halsynteny_output/mammalia/bed/rat_can_synteny_blocks_50kb_target_fix.bed" \
              "rat" \
              "can" \
              "Rattus norvegicus" \
              "Canis lupus familiaris" \
              "mammalia"

# can-mus
create_tables "results/genes_in_blocks/mammalia/can_mus_genes_in_blocks_50kb_query_fix.bed" \
              "results/genes_in_blocks/mammalia/can_mus_genes_in_blocks_50kb_target_fix.bed" \
              "results/halsynteny_output/mammalia/bed/can_mus_synteny_blocks_50kb_query_fix.bed" \
              "results/halsynteny_output/mammalia/bed/can_mus_synteny_blocks_50kb_target_fix.bed" \
              "can" \
              "mus" \
              "Canis lupus familiaris" \
              "Mus musculus" \
              "mammalia"

# mus-can
create_tables "results/genes_in_blocks/mammalia/mus_can_genes_in_blocks_50kb_query_fix.bed" \
              "results/genes_in_blocks/mammalia/mus_can_genes_in_blocks_50kb_target_fix.bed" \
              "results/halsynteny_output/mammalia/bed/mus_can_synteny_blocks_50kb_query_fix.bed" \
              "results/halsynteny_output/mammalia/bed/mus_can_synteny_blocks_50kb_target_fix.bed" \
              "mus" \
              "can" \
              "Mus musculus" \
              "Canis lupus familiaris" \
              "mammalia"

# mac-rat
create_tables "results/genes_in_blocks/mammalia/mac_rat_genes_in_blocks_50kb_query_fix.bed" \
              "results/genes_in_blocks/mammalia/mac_rat_genes_in_blocks_50kb_target_fix.bed" \
              "results/halsynteny_output/mammalia/bed/mac_rat_synteny_blocks_50kb_query_fix.bed" \
              "results/halsynteny_output/mammalia/bed/mac_rat_synteny_blocks_50kb_target_fix.bed" \
              "mac" \
              "rat" \
              "Macaca fascicularis" \
              "Rattus norvegicus" \
              "mammalia"    

# rat-mac
create_tables "results/genes_in_blocks/mammalia/rat_mac_genes_in_blocks_50kb_query_fix.bed" \
              "results/genes_in_blocks/mammalia/rat_mac_genes_in_blocks_50kb_target_fix.bed" \
              "results/halsynteny_output/mammalia/bed/rat_mac_synteny_blocks_50kb_query_fix.bed" \
              "results/halsynteny_output/mammalia/bed/rat_mac_synteny_blocks_50kb_target_fix.bed" \
              "rat" \
              "mac" \
              "Rattus norvegicus" \
              "Macaca fascicularis" \
              "mammalia"

# mac-mus
create_tables "results/genes_in_blocks/mammalia/mac_mus_genes_in_blocks_50kb_query_fix.bed" \
              "results/genes_in_blocks/mammalia/mac_mus_genes_in_blocks_50kb_target_fix.bed" \
              "results/halsynteny_output/mammalia/bed/mac_mus_synteny_blocks_50kb_query_fix.bed" \
              "results/halsynteny_output/mammalia/bed/mac_mus_synteny_blocks_50kb_target_fix.bed" \
              "mac" \
              "mus" \
              "Macaca fascicularis" \
              "Mus musculus" \
              "mammalia"

# mus-mac
create_tables "results/genes_in_blocks/mammalia/mus_mac_genes_in_blocks_50kb_query_fix.bed" \
              "results/genes_in_blocks/mammalia/mus_mac_genes_in_blocks_50kb_target_fix.bed" \
              "results/halsynteny_output/mammalia/bed/mus_mac_synteny_blocks_50kb_query_fix.bed" \
              "results/halsynteny_output/mammalia/bed/mus_mac_synteny_blocks_50kb_target_fix.bed" \
              "mus" \
              "mac" \
              "Mus musculus" \
              "Macaca fascicularis" \
              "mammalia"

# rat-mus
create_tables "results/genes_in_blocks/mammalia/rat_mus_genes_in_blocks_50kb_query_fix.bed" \
              "results/genes_in_blocks/mammalia/rat_mus_genes_in_blocks_50kb_target_fix.bed" \
              "results/halsynteny_output/mammalia/bed/rat_mus_synteny_blocks_50kb_query_fix.bed" \
              "results/halsynteny_output/mammalia/bed/rat_mus_synteny_blocks_50kb_target_fix.bed" \
              "rat" \
              "mus" \
              "Rattus norvegicus" \
              "Mus musculus" \
              "mammalia"

# mus-rat
create_tables "results/genes_in_blocks/mammalia/mus_rat_genes_in_blocks_50kb_query_fix.bed" \
              "results/genes_in_blocks/mammalia/mus_rat_genes_in_blocks_50kb_target_fix.bed" \
              "results/halsynteny_output/mammalia/bed/mus_rat_synteny_blocks_50kb_query_fix.bed" \
              "results/halsynteny_output/mammalia/bed/mus_rat_synteny_blocks_50kb_target_fix.bed" \
              "mus" \
              "rat" \
              "Mus musculus" \
              "Rattus norvegicus" \
              "mammalia"

echo "Finished script"              
