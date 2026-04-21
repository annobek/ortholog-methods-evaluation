#!/bin/bash
# Name of the job
#SBATCH --job-name=debug_alignment_depth
#SBATCH --output=slurm_outputs/debug_alignment_depth%j.out
#SBATCH --error=slurm_outputs/debug_alignment_depth%j.err

echo "Starting the script"

# Check if this region exists in HAL
#halStats --sequenceStats Canis_lupus results/output_alignment_sexspecies.hal | grep "NC_049229.1"

# 3. Try to get ANY depth data for this region
halAlignmentDepth results/output_alignment_sexspecies.hal Canis_lupus \
  --refSequence NC_049229.1 \
  --start 63782314 \
  --length 6018 \
  --outWiggle results/evaluation/mammalia/vs_prot_syn/nscore_approach/disagreement_investigation/false_positives/alignment_depth/test_612383.wig

echo "Finished script"
