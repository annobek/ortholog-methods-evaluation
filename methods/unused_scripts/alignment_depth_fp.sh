#!/bin/bash
#SBATCH --job-name=alignment_depth
#SBATCH --output=slurm_outputs/alignment_depth%j.out
#SBATCH --error=slurm_outputs/alignment_depth%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

echo "Starting the script"
echo "Timestamp: $(date)"
echo ""


#==============================================================================
# CONFIGURATION
#==============================================================================

HAL="results/output_alignment_sexspecies.hal"

# Global variables
CLASSIFICATION_FILE=""
QUERY_BED=""
QUERY_GENOME=""
SPECIES_PAIR=""
OUTPUT_DIR=""
DEPTH_DIR=""
LOG_FILE=""

#==============================================================================
# UTILITY FUNCTIONS
#==============================================================================

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg"
    if [ -n "$LOG_FILE" ]; then
        echo "$msg" >> "$LOG_FILE"
    fi
}

error() {
    local msg="[ERROR] $*"
    echo "$msg" >&2
    if [ -n "$LOG_FILE" ]; then
        echo "$msg" >> "$LOG_FILE"
    fi
}

create_dir() {
    local dir=$1
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        log "Created directory: $dir"
    fi
}

#==============================================================================
# PARAMETER SETTING
#==============================================================================

set_params() {
    CLASSIFICATION_FILE=$1
    QUERY_BED=$2
    QUERY_GENOME=$3
    SPECIES_PAIR=$4
    OUTPUT_DIR=$5
    
    DEPTH_DIR="${OUTPUT_DIR}/depths_${SPECIES_PAIR}"
    LOG_FILE="${OUTPUT_DIR}/depth_extraction_${SPECIES_PAIR}.log"
}

#==============================================================================
# COORDINATE EXTRACTION
#==============================================================================

get_gene_coordinates() {
    local bed_file=$1
    local gene_id=$2
    
    awk -v gene="$gene_id" '
        BEGIN {OFS="\t"}
        /^#/ {next}
        /^track/ {next}
        /^browser/ {next}
        $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ {next}
        $4 == gene {
            print $1, $2, $3
            exit
        }
    ' "$bed_file"
}

#==============================================================================
# DEPTH EXTRACTION
#==============================================================================

run_hal_alignment_depth() {
    local query_genome=$1
    local query_chr=$2
    local query_start=$3
    local query_length=$4
    local output_wig=$5
    
    if [ "$query_length" -le 0 ]; then
        error "Invalid gene length: $query_length"
        return 1
    fi
    
    halAlignmentDepth "$HAL" "$query_genome" \
        --refSequence "$query_chr" \
        --start "$query_start" \
        --length "$query_length" \
        --noAncestors \
        --outWiggle "$output_wig" 2>/dev/null
    
    return $?
}

#==============================================================================
# MAIN PROCESSING
#==============================================================================

process_pair() {
    local query_gene=$1
    
    # Check if already processed
    local wig_file="${DEPTH_DIR}/${query_gene}_depth.wig"
    if [ -f "$wig_file" ]; then
        return 0
    fi
    
    # Get coordinates
    local coords
    coords=$(get_gene_coordinates "$QUERY_BED" "$query_gene")
    
    if [ -z "$coords" ]; then
        error "Query gene $query_gene not found in BED"
        return 1
    fi
    
    local chr start end length
    read -r chr start end <<< "$coords"
    length=$((end - start))
    
    if [ "$length" -le 0 ]; then
        error "Invalid coordinates for gene $query_gene"
        return 1
    fi
    
    # Run halAlignmentDepth
    if ! run_hal_alignment_depth "$QUERY_GENOME" "$chr" "$start" "$length" "$wig_file"; then
        error "halAlignmentDepth failed for gene $query_gene"
        return 1
    fi
    
    return 0
}

#==============================================================================
# MAIN
#==============================================================================

main() {
    log "=========================================="
    log "Starting depth extraction"
    log "=========================================="
    
    create_dir "$OUTPUT_DIR"
    create_dir "$DEPTH_DIR"
    touch "$LOG_FILE"
    
    log "Configuration:"
    log "  Classification file: $CLASSIFICATION_FILE"
    log "  Query BED: $QUERY_BED"
    log "  Query genome: $QUERY_GENOME"
    log "  Species pair: $SPECIES_PAIR"
    log "  Output: $DEPTH_DIR"
    log ""
    
    local total=0
    local success=0
    local skipped=0
    local errors=0
    
    # Get unique query genes from C_contradictory and C_anomaly pairs
    log "Extracting unique query genes from C_contradictory and C_anomaly pairs..."
    
    local unique_genes
    unique_genes=$(awk -F'\t' 'NR>1 && ($5=="C_contradictory" || $5=="C_anomaly") {print $3}' "$CLASSIFICATION_FILE" | sort -u)
    
    local num_genes=$(echo "$unique_genes" | wc -l)
    log "Found $num_genes unique query genes to process"
    log ""
    
    # Process each unique gene
    while read -r query_gene; do
        if [ -z "$query_gene" ]; then
            continue
        fi
        
        total=$((total + 1))
        
        local wig_file="${DEPTH_DIR}/${query_gene}_depth.wig"
        
        if [ -f "$wig_file" ]; then
            skipped=$((skipped + 1))
        elif process_pair "$query_gene"; then
            success=$((success + 1))
        else
            errors=$((errors + 1))
        fi
        
        if [ $((total % 50)) -eq 0 ]; then
            log "Progress: $total genes ($success extracted, $skipped skipped, $errors errors)"
        fi
        
    done <<< "$unique_genes"
    
    # Final progress
    if [ $((total % 50)) -ne 0 ] && [ $total -gt 0 ]; then
        log "Final: $total genes ($success extracted, $skipped skipped, $errors errors)"
    fi
    
    log ""
    log "=========================================="
    log "Depth extraction complete"
    log "=========================================="
    log "Total genes: $total"
    log "Extracted: $success"
    log "Skipped (already exist): $skipped"
    log "Errors: $errors"
    log "WIG files saved to: $DEPTH_DIR"
}

#==============================================================================
# EXECUTE FOR MULTIPLE SPECIES PAIRS
#==============================================================================

mkdir -p slurm_outputs

echo ""
echo "============================================================"
echo "DEPTH EXTRACTION - MULTIPLE SPECIES PAIRS"
echo "============================================================"
echo "Start time: $(date)"
echo ""

# can-mac
#echo "=========================================="
#echo "STARTING SPECIES PAIR: can-mac"
#echo "=========================================="
#set_params \
#    results/evaluation/mammalia/vs_prot_syn/nscore_approach/disagreement_investigation/false_positives/alignment_classification_can_mac.tsv \
#    results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6 \
#    Canis_lupus \
#    can_mac \
#    results/evaluation/mammalia/vs_prot_syn/nscore_approach/disagreement_investigation/false_positives/alignment_depth
#main


# can-rat
echo ""
echo "=========================================="
echo "STARTING SPECIES PAIR: can-rat"
echo "=========================================="
set_params \
    results/evaluation/mammalia/vs_prot_syn/nscore_approach/disagreement_investigation/false_positives/alignment_classification_can_rat.tsv \
    results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6 \
    Canis_lupus \
    can_rat \
    results/evaluation/mammalia/vs_prot_syn/nscore_approach/disagreement_investigation/false_positives/alignment_depth
main

# can-mus
echo ""
echo "=========================================="
echo "STARTING SPECIES PAIR: can-mus"
echo "=========================================="
set_params \
    results/evaluation/mammalia/vs_prot_syn/nscore_approach/disagreement_investigation/false_positives/alignment_classification_can_mus.tsv \
    results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6 \
    Canis_lupus \
    can_mus \
    results/evaluation/mammalia/vs_prot_syn/nscore_approach/disagreement_investigation/false_positives/alignment_depth
main

# mac-rat
echo ""
echo "=========================================="
echo "STARTING SPECIES PAIR: mac-rat"
echo "=========================================="
set_params \
    results/evaluation/mammalia/vs_prot_syn/nscore_approach/disagreement_investigation/false_positives/alignment_classification_mac_rat.tsv \
    results/gene_only_annotation/mammalia/Macaca_genes_only_stableid_fix.bed6 \
    Macaca_fascicularis \
    mac_rat \
    results/evaluation/mammalia/vs_prot_syn/nscore_approach/disagreement_investigation/false_positives/alignment_depth
main

# mac-mus
echo ""
echo "=========================================="
echo "STARTING SPECIES PAIR: mac-mus"
echo "=========================================="
set_params \
    results/evaluation/mammalia/vs_prot_syn/nscore_approach/disagreement_investigation/false_positives/alignment_classification_mac_mus.tsv \
    results/gene_only_annotation/mammalia/Macaca_genes_only_stableid_fix.bed6 \
    Macaca_fascicularis \
    mac_mus \
    results/evaluation/mammalia/vs_prot_syn/nscore_approach/disagreement_investigation/false_positives/alignment_depth
main

# rat-mus
echo ""
echo "=========================================="
echo "STARTING SPECIES PAIR: rat-mus"
echo "=========================================="
set_params \
    results/evaluation/mammalia/vs_prot_syn/nscore_approach/disagreement_investigation/false_positives/alignment_classification_rat_mus.tsv \
    results/gene_only_annotation/mammalia/Rat_genes_only_stableid_fix.bed6 \
    Rattus_norvegicus \
    rat_mus \
    results/evaluation/mammalia/vs_prot_syn/nscore_approach/disagreement_investigation/false_positives/alignment_depth
main



echo ""
echo "============================================================"
echo "ALL SPECIES PAIRS COMPLETED"
echo "============================================================"
echo "End time: $(date)"
echo ""
