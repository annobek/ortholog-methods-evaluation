#!/bin/bash
# Name of the job
#SBATCH --job-name=get_parameters
#SBATCH --output=get_parameters%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# Computes intergenic median distance and suggested parameters for halSynteny
# from full GTF files (not prefiltered).

set -euo pipefail

# List of GTFs (adjust if filenames differ)
GTFS=("material/gtf_Canis.gtf" "material/gtf_Macaca.gtf" "material/gtf_Mus.gtf" "material/gtf_Rat.gtf")


# Function to process a GTF
process_gtf () {
    local gtf=$1
    local base=$(basename "$gtf")

    echo "Processing $gtf ..."

    # 1. Extract gene coordinates to BED
    awk '$3 == "gene" {print $1"\t"($4-1)"\t"$5}' "$gtf" > "${base}.genes.bed"

    # 2. Sort
    sort -k1,1 -k2,2n "${base}.genes.bed" > "${base}.genes.sorted.bed"

    # 3. Compute intergenic distances
    awk 'BEGIN{OFS="\t"} {
        chr=$1; start=$2; end=$3;
        if(chr==prev_chr){
            dist = start - prev_end;
            if(dist < 0) dist = 0;
            print dist;
        }
        prev_chr=chr; prev_end=end;
    }' "${base}.genes.sorted.bed" > "${base}.intergenic.distances.txt"

    # 4. Compute median and W=20*median (using Python for robustness)
    python3 - "$base" <<PY
import sys, statistics
base = sys.argv[1]
vals = [int(x.strip()) for x in open(f"{base}.intergenic.distances.txt") if x.strip()]
if not vals:
    print(f"{base}\tNA\tNA\tNA")
else:
    med = int(statistics.median(vals))
    W = med * 20
    minBlock = med*2
    print(f"{base}\tmedian={med}\tW={W}\tminBlockSize={minBlock}")
PY
}

echo -e "Species\tMedian_interg_bp\tWindow20genes_bp\tminBlockSize_suggested"
for gtf in "${GTFS[@]}"; do
    process_gtf "$gtf"
done

rm *.bed
rm *.txt
