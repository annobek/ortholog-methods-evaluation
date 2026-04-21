#!/bin/bash
#SBATCH --job-name=get_alignments
#SBATCH --output=slurm_outputs/get_alignments%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G


#hal2maf results/output_alignment_sexspecies.hal results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/test_output.maf \
#  --refGenome Canis_lupus \
#  --targetGenomes Macaca_fascicularis \
#  --refSequence NC_049260.1 \
#  --start 46813739 \
#  --length 9149 \
#  --noAncestors

# Check if it worked:
#ls -lh results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/test_output.maf
#cat results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/test_output.maf



#filename="NP_001002931.1_vs_XP_005589777.3.maf"
#outdir="results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/alignments_can_mac"
#
#echo "Testing: ${outdir}/${filename}"
#echo ""
#
## Show what the s lines actually look like
#echo "First few 's' lines in file:"
#grep "^s" "${outdir}/${filename}" | head -3 | cat -A
#echo ""
#
## Test with space
#if grep -q "^s " "${outdir}/${filename}"; then
#    echo "✓ Found with 's ' (with space)"
#else
#    echo "✗ Not found with 's ' (with space)"
#fi
#
## Test without space
#if grep -q "^s" "${outdir}/${filename}"; then
#    echo "✓ Found with 's' (no space)"
#else
#    echo "✗ Not found with 's' (no space)"
#fi
#
## Count them
#echo ""
#echo "Lines starting with 's ' (with space): $(grep -c "^s " "${outdir}/${filename}")"
#echo "Lines starting with 's' (no space): $(grep -c "^s" "${outdir}/${filename}")"


# Show first few lines with column numbers
head -5 results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6 | awk '{print "Col1:"$1, "Col2:"$2, "Col3:"$3, "Col4:"$4, "Col5:"$5, "Col6:"$6}'

# Check what one of the "missing" genes actually looks like in the mapping
head -1 results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/can_mac_missing_genes.txt

# Now try to find that gene in the BED file
MISSING_GENE=$(head -1 results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/cactus_support_score/false_positives/can_mac_missing_genes.txt)
echo "Looking for gene: $MISSING_GENE"
grep "$MISSING_GENE" results/gene_only_annotation/mammalia/Canis_genes_only_stableid_fix.bed6 | head -1

# Also check what's in the mapping for this gene
grep "$MISSING_GENE" results/evaluation/mammalia/vs_neighborhood/disagreement_investigation/gene_to_protein_map_FIXED.tsv
