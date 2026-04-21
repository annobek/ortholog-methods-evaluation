#!/bin/bash
#SBATCH --job-name=extract_geneid_biotype_from_annotation_mammalia
#SBATCH --output=slurm_outputs/extract_geneid_biotype_from_annotation_mammalia%j.out


echo "Starting GeneID extraction only..."

mkdir -p results/gene_only_annotation

extract_geneid () {
    local gtf_file=$1
    local out_file=$2

    echo "Processing $gtf_file - $out_file"
    grep -w "gene" "$gtf_file" | \
    awk -F'\t|;' '{
        g=""; b="";
        for (i=1; i<=NF; i++) {
            if ($i ~ /db_xref/ && $i ~ /GeneID:/)  {g=$i}
            if ($i ~ /gene_biotype/) {b=$i}
        }
        sub(/.*GeneID:/, "", g);
        sub(/".*/, "", g);
        sub(/.*gene_biotype "/, "", b);
        sub(/".*/, "", b);
        if (g != "" && b != "") print g "\t" b;
    }' > "$out_file"

    echo "Extracted $(wc -l < "$out_file") entries."
}

extract_geneid results/gene_only_annotation/Canis_genes_only_fix.gtf  results/gene_only_annotation/can_id_biotype_fix.tsv
extract_geneid results/gene_only_annotation/Macaca_genes_only_fix.gtf results/gene_only_annotation/mac_id_biotype_fix.tsv
extract_geneid results/gene_only_annotation/Rat_genes_only_fix.gtf    results/gene_only_annotation/rat_id_biotype_fix.tsv
extract_geneid results/gene_only_annotation/Mus_genes_only_fix.gtf    results/gene_only_annotation/mus_id_biotype_fix.tsv

echo "Finished extracting all GeneID-based annotations."
