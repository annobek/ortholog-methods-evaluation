#!/bin/bash
# Name of the job
#SBATCH --job-name=gtf_to_bed
#SBATCH --output=slurm_outputs/gtf_to_bed%j.out

# Convert GTF files to BED for plants and drosophila

echo "Starting the script"


#########################
# DROSOPHILA (GTF)
#########################
# Expect attributes like:
#   gene_id "FBgn0031081"; gene_symbol "Nep3";
# We want: FBgn0031081 (exactly as is).

convert_gtf_to_bed_droso() {
  local input_gtf=$1
  local output_bed6=$2

  echo "Converting (Drosophila) $input_gtf -> $output_bed6 ..."

  awk -F'\t' 'BEGIN{OFS="\t"}
    $1 !~ /^#/ {
      chrom  = $1
      start  = $4 - 1    # BED is 0-based
      end    = $5
      strand = $7
      attrs  = $9        # full attribute column (now correct because FS = tab)

      id = ""

      # gene_id "FBgn0031081"
      if (match(attrs, /gene_id[[:space:]]*"([^"]+)"/, m)) {
        id = m[1]        # KEEP EXACT FBgn ID
      }

      if (id == "")
        id = "unknown_" NR

      print chrom, start, end, id, 0, strand
    }' "$input_gtf" > "$output_bed6"

  echo "Malformed strand values (should be 0):"
  awk '$6!="+" && $6!="-" && $6!="."' "$output_bed6" | wc -l
}


#########################
# PLANTS (GFF3-like)
#########################
# Arabidopsis example:
#   ID=gene:AT1G01010;Name=...
#   -> AT1G01010
#
# Cardamine example:
#   ID=CARHR000120;Description=...
#   -> CARHR000120

convert_gtf_to_bed_plants() {
  local input_gtf=$1
  local output_bed6=$2

  echo "Converting (Plant) $input_gtf -> $output_bed6 ..."

  awk -F'\t' 'BEGIN{OFS="\t"}
    $1 !~ /^#/ {
      chrom  = $1
      start  = $4 - 1
      end    = $5
      strand = $7
      attrs  = $9

      id = ""

      # ID=gene:AT1G01010 or ID=CARHR000120
      if (match(attrs, /ID=([^;]+)/, m)) {
        id = m[1]
        sub(/^gene:/, "", id)   # only strip "gene:" prefix if present
      }

      if (id == "")
        id = "unknown_" NR

      print chrom, start, end, id, 0, strand
    }' "$input_gtf" > "$output_bed6"

  echo "Malformed strand values (should be 0):"
  awk '$6!="+" && $6!="-" && $6!="."' "$output_bed6" | wc -l
}




# Plants
#-------

# Arabidopsis
convert_gtf_to_bed_plants "results/gene_only_annotation/brassicaceae/Arabidopsis_genes_only.gtf" \
                   "results/gene_only_annotation/brassicaceae/Arabidopsis_genes_only.bed6"

# Cardamine
convert_gtf_to_bed_plants "results/gene_only_annotation/brassicaceae/Cardamine_genes_only.gtf" \
                   "results/gene_only_annotation/brassicaceae/Cardamine_genes_only.bed6"


# Drosophila
#-----------

# dmel
convert_gtf_to_bed_droso "results/gene_only_annotation/drosophila/dmel_genes_only.gtf" \
                   "results/gene_only_annotation/drosophila/dmel_genes_only.bed6"

# dsim
convert_gtf_to_bed_droso "results/gene_only_annotation/drosophila/dsim_genes_only.gtf" \
                   "results/gene_only_annotation/drosophila/dsim_genes_only.bed6"

# dsec
convert_gtf_to_bed_droso "results/gene_only_annotation/drosophila/dsec_genes_only.gtf" \
                   "results/gene_only_annotation/drosophila/dsec_genes_only.bed6"


echo "Finished script"
