#!/bin/bash
# Name of the job
#SBATCH --job-name=get_parameters
#SBATCH --output=get_parameters%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# Computes 90th percentile of intergenic distances and parameters for halSynteny
# from full GFF files (not prefiltered).


set -euo pipefail

mkdir -p results/halsynteny_parameters

# List of GFF3 files (adjust filenames/paths if needed)

# Plants
#GFFS=("material/cardamine/Arabidopsis.gff3" "material/cardamine/Cardamine.gff")

# Drosophila
GFFS=("material/dmel.gtf" "material/dsec.gtf" "material/dsim.gtf")

# Function to process one GFF/GTF
process_gtf () {
    local gff=$1
    local base=$(basename "$gff")

    echo "Processing $gff ..."

    # 1. Extract gene coordinates
    awk '!/^#/ && $3 == "gene" {print $1"\t"($4-1)"\t"$5}' "$gff" > "results/halsynteny_parameters/${base}.genes.bed"

    # 2. Sort by chromosome and start position
    sort -k1,1 -k2,2n "results/halsynteny_parameters/${base}.genes.bed" > "results/halsynteny_parameters/${base}.genes.sorted.bed"

    # 3. Compute intergenic distances
    awk 'BEGIN{OFS="\t"} {
        chr=$1; start=$2; end=$3;
        if(chr==prev_chr){
            dist = start - prev_end;
            if(dist < 0) dist = 0;
            print dist;
        }
        prev_chr=chr; prev_end=end;
    }' "results/halsynteny_parameters/${base}.genes.sorted.bed" > "results/halsynteny_parameters/${base}.intergenic.distances.txt"

    # 4. Compute 90th percentile and derived parameters
    python3 - "$base" <<'PY'
import sys, numpy as np, statistics

base = sys.argv[1]
vals = [int(x.strip()) for x in open(f"results/halsynteny_parameters/{base}.intergenic.distances.txt") if x.strip()]
if not vals:
    print(f"{base}\tNA\tNA\tNA")
else:
    p90 = int(np.percentile(vals, 90))
    minBlock = int(p90 * 0.5)
    print(f"{base}\tP90={p90}\tmaxAnchorDistance={p90}\tminBlockSize={minBlock}")
PY
}

echo -e "Species\tP90_intergenic_bp\tmaxAnchorDistance_bp\tminBlockSize_bp"
for gff in "${GFFS[@]}"; do
    process_gtf "$gff"
done



