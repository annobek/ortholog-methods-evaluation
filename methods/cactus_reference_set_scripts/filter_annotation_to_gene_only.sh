#!/bin/bash
# Name of the job
#SBATCH --job-name=filter_ann
# output file name
#SBATCH --output=slurm_outputs/filter_ann%j.out

# Filter the GTF files to keep only gene entries

echo "Starting the script"

mkdir -p results/gene_only_annotation/mammalia

filter_ann_gene () {
  local original_ann=$1
  local filtered_ann=$2

  # Extracting only gene entries and excluding pseudogenes
  awk -F'\t' 'BEGIN{IGNORECASE=1}
  $3=="gene" && $0!~/^#/ {
    bt="";
    if (match($0,/gene_biotype "([^"]+)"/,m)) bt=m[1];
    else if (match($0,/gene_type "([^"]+)"/,m)) bt=m[1];
    if (bt !~ /pseudogene/)
      print $0;
  }' $original_ann > $filtered_ann

}


# Mammalia below


## Brassicaceae
#filter_ann_gene "material/cardamine/Arabidopsis.gff3" "results/gene_only_annotation/brassicaceae/Arabidopsis_genes_only.gtf"
#filter_ann_gene "material/cardamine/Cardamine.gff" "results/gene_only_annotation/brassicaceae/Cardamine_genes_only.gtf"

## Drosophila
#filter_ann_gene "material/dmel.gtf" "results/gene_only_annotation/drosophila/dmel_genes_only.gtf"
#filter_ann_gene "material/dsec.gtf" "results/gene_only_annotation/drosophila/dsec_genes_only.gtf"
#filter_ann_gene "material/dsim.gtf" "results/gene_only_annotation/drosophila/dsim_genes_only.gtf"

# Mammals
filter_ann_gene "material/sex_experiment/gtf_Canis.gtf" "results/gene_only_annotation/mammalia/Canis_genes_only_fix.gtf"
filter_ann_gene "material/sex_experiment/gtf_Macaca.gtf" "results/gene_only_annotation/mammalia/Macaca_genes_only_fix.gtf"
filter_ann_gene "material/sex_experiment/gtf_Rat.gtf" "results/gene_only_annotation/mammalia/Rat_genes_only_fix.gtf"
filter_ann_gene "material/sex_experiment/gtf_Mus.gtf" "results/gene_only_annotation/mammalia/Mus_genes_only_fix.gtf"


# Filter the GTF files to keep only CDS entries

# old approach
#awk -F'\t' '$3 == "CDS"' material/sex_experiment/gtf_Canis.gtf > results/ann_cds/Canis_cds_only.gtf
#awk -F'\t' '$3 == "CDS"' material/sex_experiment/gtf_Macaca.gtf > results/ann_cds/Macaca_cds_only.gtf
#awk -F'\t' '$3 == "CDS"' material/sex_experiment/gtf_Mus.gtf > results/ann_cds/Mus_cds_only.gtf
#awk -F'\t' '$3 == "CDS"' material/sex_experiment/gtf_Rat.gtf > results/ann_cds/Rat_cds_only.gtf

#awk -F'\t' '$3 == "CDS"' material/dmel.gtf > results/ann_cds/dmel_cds_only.gtf
#awk -F'\t' '$3 == "CDS"' material/dsec.gtf > results/ann_cds/dsec_cds_only.gtf
#awk -F'\t' '$3 == "CDS"' material/dsim.gtf > results/ann_cds/dsim_cds_only.gtf

echo "Finished script"
