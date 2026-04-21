#!/bin/bash
#SBATCH --job-name=get_alignments
#SBATCH --output=slurm_outputs/get_alignments%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

echo "Starting the script"
echo "Timestamp: $(date)"
echo ""

# Create output directories
mkdir -p results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/alignments_can_mac/
mkdir -p results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/alignments_can_rat/
mkdir -p results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/alignments_can_mus/
mkdir -p results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/alignments_mac_rat/
mkdir -p results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/alignments_mac_mus/
mkdir -p results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/alignments_rat_mus/

HAL=results/output_alignment_sexspecies.hal
PROTEIN_GENE_MAP="results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/gene_to_protein_map_FIXED.tsv"

# Verify files exist
if [ ! -f "$HAL" ]; then
    echo "ERROR: HAL file not found: $HAL"
    exit 1
fi

if [ ! -f "$PROTEIN_GENE_MAP" ]; then
    echo "ERROR: Protein-gene mapping file not found: $PROTEIN_GENE_MAP"
    exit 1
fi

echo "HAL file: $HAL"
echo "Protein-gene map: $PROTEIN_GENE_MAP"
echo ""

# Function to convert protein ID to gene ID
protein_to_gene() {
    local protein_id=$1
    local mapping_file=$2
    gene_id=$(awk -v prot="$protein_id" '$2 == prot {print $1; exit}' "$mapping_file")
    echo "$gene_id"
}

# Function to extract alignments for a species pair
extract_alignments() {
    local src_species=$1
    local tgt_species=$2
    local src_bed=$3
    local pairs_file=$4
    local outdir=$5
    local pair_name=$6
    
    echo "============================================"
    echo "Processing $pair_name: $src_species → $tgt_species"
    echo "============================================"
    
    # Verify input files exist
    if [ ! -f "$src_bed" ]; then
        echo "ERROR: BED file not found: $src_bed"
        return 1
    fi
    
    if [ ! -f "$pairs_file" ]; then
        echo "ERROR: Pairs file not found: $pairs_file"
        return 1
    fi
    
    # Counters
    local total=0
    local success=0
    local missing_gene=0
    local empty=0
    local unmapped_protein=0
    
    # Create missing lists
    missing_genes_list="${outdir}/../${pair_name}_missing_genes.txt"
    unmapped_proteins_list="${outdir}/../${pair_name}_unmapped_proteins.txt"
    > "$missing_genes_list"
    > "$unmapped_proteins_list"
    
    # Read each pair and extract alignments
    while read -r src_protein tgt_protein; do
        # Skip header
        if [ "$src_protein" == "Query_Protein" ]; then
            continue
        fi
    
        ((total++))
        
        # Convert protein IDs to gene IDs
        src_gene=$(protein_to_gene "$src_protein" "$PROTEIN_GENE_MAP")
        tgt_gene=$(protein_to_gene "$tgt_protein" "$PROTEIN_GENE_MAP")
        
        # Check if protein-to-gene mapping was successful
        if [ -z "$src_gene" ]; then
            echo "$src_protein" >> "$unmapped_proteins_list"
            ((unmapped_protein++))
            continue
        fi
        
        # Extract chr, start, end from source BED using GENE ID
        line=$(awk -v gene="$src_gene" '$4 == gene {print; exit}' "$src_bed")
        
        if [ -z "$line" ]; then
            echo "$src_gene" >> "$missing_genes_list"
            ((missing_gene++))
            continue
        fi
        
        chr=$(echo $line | awk '{print $1}')
        start=$(echo $line | awk '{print $2}')
        end=$(echo $line | awk '{print $3}')
        length=$((end - start))
        
        # Create filename - use PROTEIN IDs for filename
        filename="${src_protein}_vs_${tgt_protein}.maf"
        
        # Run hal2maf for this interval
        if hal2maf $HAL "${outdir}/${filename}" \
            --refGenome $src_species \
            --targetGenomes $tgt_species \
            --refSequence $chr \
            --start $start \
            --length $length \
            --noAncestors 2>/dev/null; then
            
            # Check if MAF has actual alignment content
            if [ -s "${outdir}/${filename}" ] && grep -q "^s" "${outdir}/${filename}"; then
                ((success++))
            else
                ((empty++))
            fi
        fi
        
        # Print progress every 100 genes
        if [ $((total % 100)) -eq 0 ]; then
            echo "  Progress: $total pairs processed ($success successful, $unmapped_protein unmapped, $missing_gene missing)..."
        fi
            
    done < "$pairs_file"
    
    # Print summary
    echo ""
    echo "Summary for $pair_name:"
    echo "  Total pairs:           $total"
    echo "  Successful:            $success"
    echo "  Unmapped proteins:     $unmapped_protein"
    echo "  Missing genes in BED:  $missing_gene"
    echo "  Empty alignments:      $empty"
    echo ""
    
    # Save unmapped proteins list
    if [ -s "$unmapped_proteins_list" ]; then
        sort -u "$unmapped_proteins_list" > "${unmapped_proteins_list}.tmp"
        mv "${unmapped_proteins_list}.tmp" "$unmapped_proteins_list"
        echo "  $(wc -l < "$unmapped_proteins_list") unique unmapped proteins saved to: $unmapped_proteins_list"
    else
        rm "$unmapped_proteins_list"
    fi
    
    # Save missing genes list
    if [ -s "$missing_genes_list" ]; then
        sort -u "$missing_genes_list" > "${missing_genes_list}.tmp"
        mv "${missing_genes_list}.tmp" "$missing_genes_list"
        echo "  $(wc -l < "$missing_genes_list") unique missing genes saved to: $missing_genes_list"
    else
        rm "$missing_genes_list"
    fi
    
    echo ""
}

# Process all species pairs
extract_alignments \
    "Canis_lupus" \
    "Macaca_fascicularis" \
    "results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6" \
    "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/fp_by_species_pair/missed_by_cactus_gm_can_mac.tsv" \
    "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/alignments_can_mac/" \
    "can_mac"

extract_alignments \
    "Canis_lupus" \
    "Rattus_norvegicus" \
    "results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6" \
    "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/fp_by_species_pair/missed_by_cactus_gm_can_rat.tsv" \
    "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/alignments_can_rat/" \
    "can_rat"

extract_alignments \
    "Canis_lupus" \
    "Mus_musculus" \
    "results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6" \
    "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/fp_by_species_pair/missed_by_cactus_gm_can_mus.tsv" \
    "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/alignments_can_mus/" \
    "can_mus"

extract_alignments \
    "Macaca_fascicularis" \
    "Rattus_norvegicus" \
    "results/gene_only_annotation/mammalia/Macaca_genes_only_stableid_fix.bed6" \
    "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/fp_by_species_pair/missed_by_cactus_gm_mac_rat.tsv" \
    "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/alignments_mac_rat/" \
    "mac_rat"

extract_alignments \
    "Macaca_fascicularis" \
    "Mus_musculus" \
    "results/gene_only_annotation/mammalia/Macaca_genes_only_stableid_fix.bed6" \
    "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/fp_by_species_pair/missed_by_cactus_gm_mac_mus.tsv" \
    "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/alignments_mac_mus/" \
    "mac_mus"

extract_alignments \
    "Rattus_norvegicus" \
    "Mus_musculus" \
    "results/gene_only_annotation/mammalia/Rat_genes_only_stableid_fix.bed6" \
    "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/fp_by_species_pair/missed_by_cactus_gm_rat_mus.tsv" \
    "results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/alignments_rat_mus/" \
    "rat_mus"

echo "============================================"
echo "Finished script"
echo "Timestamp: $(date)"
echo "============================================"