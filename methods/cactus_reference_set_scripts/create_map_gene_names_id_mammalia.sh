#!/bin/bash
# Name of the job
#SBATCH --job-name=create_map_gene_names_id_mammalia
#SBATCH --output=slurm_outputs/create_map_gene_names_id_mammalia%j.out

# Creates a map where for each gene name is assigned its gene ID from db_xref 

echo "Starting the script."

mkdir -p results/gene_only_annotation/mammalia
#mkdir -p results/ann_cds

# Function for gene-only annotations
gene_only_extract() {
    local gtf_file=$1
    local map_output=$2
    
    echo "Processing gene-only: $gtf_file"

    grep -v "^#" "$gtf_file" | \
    awk -F'\t' '$3=="gene" {
        gene_symbol=""; gene_id="";
        # Parse gene_id or gene name
        if (match($9, /gene_id "([^"]+)"/, arr)) {
            gene_symbol = arr[1];
        } else if (match($9, /gene "([^"]+)"/, arr)) {
            gene_symbol = arr[1];
        }
        # Parse GeneID from db_xref
        if (match($9, /db_xref "GeneID:([0-9]+)"/, arr2)) {
            gene_id = arr2[1];
        }
        # Print mapping if both found
        if (gene_symbol != "" && gene_id != "") {
            print gene_symbol "\t" gene_id;
        }
    }' | sort -u > "$map_output"
    
    echo "Created: $map_output ($(wc -l < "$map_output") mappings)"
}


# Function for CDS-only annotations
cds_only_extract() {
    local gtf_file=$1
    local map_output=$2
    
    echo "Processing CDS-only: $gtf_file"

    awk -F'\t' -v OFS='\t' '
        # Keep only non-comment CDS rows
        $0 !~ /^#/ && $3=="CDS" {
            sym=""; gid="";
            # Parse attributes
            n=split($9, a, ";")
            for (i=1; i<=n; i++) {
                gsub(/^ +| +$/, "", a[i])
                if (sym=="" && match(a[i], /^gene_id "([^"]+)"/, m))
                    sym=m[1]
                if (gid=="" && match(a[i], /db_xref "GeneID:([0-9]+)"/, m2))
                    gid=m2[1]
            }
            if (sym!="" && gid!="") counts[sym OFS gid]++
        }
        END {
            # Collect symbols
            for (k in counts) {
                split(k, p, OFS); sym=p[1]; bysym[sym]=1
            }
            # For each symbol, pick GeneID with max count; tie-break by smallest numeric ID
            for (sym in bysym) {
                maxc=0; pick=""
                for (k in counts) {
                    split(k, p, OFS)
                    if (p[1]==sym) {
                        c=counts[k]; g=p[2]
                        if (c>maxc || (c==maxc && (pick=="" || g+0<pick+0))) {
                            maxc=c; pick=g
                        }
                    }
                }
                print sym, pick
            }
        }' "$gtf_file" | LC_ALL=C sort -u > "$map_output"
    
    echo "Created: $map_output ($(wc -l < "$map_output") mappings)"
}

# Mammals
#--------

# Canis
gene_only_extract "results/gene_only_annotation/mammalia/Canis_genes_only_fix.gtf" \
                  "results/gene_only_annotation/mammalia/Canis_gene_to_id_fix.tsv"

# Macaca
gene_only_extract "results/gene_only_annotation/mammalia/Macaca_genes_only_fix.gtf" \
                  "results/gene_only_annotation/mammalia/Macaca_gene_to_id_fix.tsv"

# Rat
gene_only_extract "results/gene_only_annotation/mammalia/Rat_genes_only_fix.gtf" \
                  "results/gene_only_annotation/mammalia/Rat_gene_to_id_fix.tsv"

# Mus
gene_only_extract "results/gene_only_annotation/mammalia/Mus_genes_only_fix.gtf" \
                  "results/gene_only_annotation/mammalia/Mus_gene_to_id_fix.tsv"

#---------------------------------------------------
# CDS-only

# Canis
#cds_only_extract "results/ann_cds/Canis_cds_only.gtf" \
#                 "results/ann_cds/Canis_cds_to_id.tsv"

# Macaca
#cds_only_extract "results/ann_cds/Macaca_cds_only.gtf" \
#                 "results/ann_cds/Macaca_cds_to_id.tsv"

# Rat
#cds_only_extract "results/ann_cds/Rat_cds_only.gtf" \
#                 "results/ann_cds/Rat_cds_to_id.tsv"

# Mus
#cds_only_extract "results/ann_cds/Mus_cds_only.gtf" \
#                 "results/ann_cds/Mus_cds_to_id.tsv"

#echo ""
#echo "=== Summary ==="
#echo "Gene mappings:"
#wc -l results/ann_genes/*_gene_to_id_fix.tsv
#echo ""
#echo "CDS mappings:"
#wc -l results/ann_cds/*_cds_to_id.tsv

#echo ""
echo "Finished script"
