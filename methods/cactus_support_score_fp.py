import sys
from pathlib import Path
from collections import defaultdict


"""
Calculate Cactus Support Score for FP pairs.

Processes C_contradictory, C_anomaly, and C_no_result pairs.
Calculates PID, Overlap, Depth normalization, and final Cactus Score.

Output format:
Query_Gene  Target_Gene  Classification  Alt_Gene  PID  Overlap  Depth  Cactus_Score  Alt_Cactus_Score
"""


#==============================================================================
# MAF PARSING - PID & OVERLAP
#==============================================================================

def parse_maf_sequences(maf_file, target_genome):
    """
    Parse MAF file and extract aligned sequences.
    
    Returns:
        Tuple of (query_seq, target_seq, query_aligned_length, target_aligned_length)
    """
    query_seq_parts = []
    target_seq_parts = []
    query_aligned = 0
    target_aligned = 0
    
    if not Path(maf_file).exists():
        return None, None, 0, 0
    
    current_query = None
    current_target = None
    
    with open(maf_file) as f:
        for line in f:
            if line.startswith('a'):
                # New alignment block - save previous
                if current_query and current_target:
                    query_seq_parts.append(current_query)
                    target_seq_parts.append(current_target)
                current_query = None
                current_target = None
            
            elif line.startswith('s'):
                parts = line.split()
                genome_chr = parts[1]
                length = int(parts[3])
                sequence = parts[6]
                
                # First sequence is query
                if current_query is None:
                    current_query = sequence
                    query_aligned += length
                # Check if this is our target genome
                elif genome_chr.startswith(target_genome + '.'):
                    current_target = sequence
                    target_aligned += length
        
        # Save last block
        if current_query and current_target:
            query_seq_parts.append(current_query)
            target_seq_parts.append(current_target)
    
    if not query_seq_parts:
        return None, None, 0, 0
    
    # Concatenate
    query_full = ''.join(query_seq_parts)
    target_full = ''.join(target_seq_parts)
    
    return query_full, target_full, query_aligned, target_aligned


def calculate_pid(seq1, seq2):
    """Calculate percent identity (ignoring gaps)."""
    if not seq1 or not seq2:
        return None
    
    if len(seq1) != len(seq2):
        return None
    
    matches = 0
    total = 0
    
    for base1, base2 in zip(seq1, seq2):
        # Skip gaps
        if base1 == '-' or base2 == '-':
            continue
        
        total += 1
        if base1.upper() == base2.upper():
            matches += 1
    
    if total == 0:
        return None
    
    return matches / total


def calculate_overlap(query_aligned, target_aligned, query_length, target_length):
    """Calculate overlap fraction."""
    if query_length == 0 or target_length == 0:
        return None
    
    overlap = (query_aligned + target_aligned) / (query_length + target_length)
    return overlap


#==============================================================================
# DEPTH PARSING & NORMALIZATION
#==============================================================================

def parse_wig_depth(wig_file):
    """Parse WIG file and return list of depth values."""
    depths = []
    
    if not Path(wig_file).exists():
        return None
    
    with open(wig_file) as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('fixedStep') or line.startswith('variableStep') or \
               line.startswith('track') or line.startswith('#') or not line:
                continue
            
            try:
                depths.append(int(line))
            except ValueError:
                continue
    
    return depths if depths else None


def normalize_depth(depths, n_species):
    """
    Normalize alignment depth to [0, 1].
    
    Formula: DepthNorm = (MeanDepth - 1) / (N_species - 1)
    """
    if not depths or n_species <= 1:
        return None
    
    mean_depth = sum(depths) / len(depths)
    depth_norm = (mean_depth - 1) / (n_species - 1)
    
    # Clamp to [0, 1]
    return max(0.0, min(1.0, depth_norm))


#==============================================================================
# CACTUS SCORE CALCULATION
#==============================================================================

def calculate_cactus_score(pid, overlap, depth_norm):
    """
    Calculate Cactus Support Score.
    
    Score = (PID * Overlap * DepthNorm)^(1/3)
    """
    if pid is None or overlap is None or depth_norm is None:
        return None
    
    score = (pid * overlap * depth_norm) ** (1/3)
    return score


#==============================================================================
# GENE LENGTH EXTRACTION
#==============================================================================

def load_gene_lengths(bed_file):
    """Load gene lengths from BED file."""
    lengths = {}
    
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            
            try:
                gene_id = parts[3]
                start = int(parts[1])
                end = int(parts[2])
                lengths[gene_id] = end - start
            except ValueError:
                continue
    
    return lengths


#==============================================================================
# MAIN PROCESSING
#==============================================================================

def process_pair(
    pair,
    alignment_dir,
    depth_dir,
    query_lengths,
    target_lengths,
    target_genome,
    n_species
):
    """
    Process one pair and calculate all metrics.
    
    Returns updated pair dict with metrics.
    """
    query_gene = pair['query_gene']
    target_gene = pair['target_gene']
    classification = pair['classification']
    alt_gene = pair.get('alt_gene', '-')
    
    # Initialize metrics
    metrics = {
        'pid': None,
        'overlap': None,
        'depth_norm': None,
        'cactus_score': None,
        'alt_cactus_score': None
    }
    
    # C_no_result: No alignment
    if classification == 'C_no_result':
        metrics['cactus_score'] = 0.0
        pair.update(metrics)
        return pair
    
    # Get MAF file
    maf_file = alignment_dir / f"{query_gene}_vs_{target_gene}.maf"
    
    # DEBUG
    # print(f"  Processing {query_gene} -> {target_gene}", file=sys.stderr)
    # print(f"    MAF: {maf_file}", file=sys.stderr)
    
    # Parse MAF for PID and overlap
    query_seq, target_seq, query_aligned, target_aligned = \
        parse_maf_sequences(maf_file, target_genome)
    
    if query_seq:
        metrics['pid'] = calculate_pid(query_seq, target_seq)
        # print(f"    PID: {metrics['pid']}", file=sys.stderr)
    
    # Get gene lengths
    query_length = query_lengths.get(query_gene)
    
    # For contradictory: use alternative gene length
    # For anomaly: use predicted target gene length
    if classification == 'C_contradictory' and alt_gene != '-':
        target_length = target_lengths.get(alt_gene)
    else:
        target_length = target_lengths.get(target_gene)
    
    if query_length and target_length:
        metrics['overlap'] = calculate_overlap(
            query_aligned, target_aligned, query_length, target_length
        )
        # print(f"    Overlap: {metrics['overlap']}", file=sys.stderr)
    
    # Get depth
    wig_file = depth_dir / f"{query_gene}_depth.wig"
    
    # DEBUG - UNCOMMENT THESE
    #print(f"    Looking for WIG: {wig_file}", file=sys.stderr)
    #print(f"    WIG exists: {wig_file.exists()}", file=sys.stderr)
    
    depths = parse_wig_depth(wig_file)
    
    if depths:
        metrics['depth_norm'] = normalize_depth(depths, n_species)
        #print(f"    DepthNorm: {metrics['depth_norm']}", file=sys.stderr)
    else:
        print(f"    No depth data!", file=sys.stderr)
    
    # Calculate Cactus score
    metrics['cactus_score'] = calculate_cactus_score(
        metrics['pid'], metrics['overlap'], metrics['depth_norm']
    )
    
    # For contradictory pairs: calculate score for alternative gene alignment
    if classification == 'C_contradictory':
        metrics['alt_cactus_score'] = metrics['cactus_score']
        # Set main score to 0 (alignment to predicted target doesn't exist)
        metrics['cactus_score'] = 0.0
    
    pair.update(metrics)
    return pair


def load_classification_pairs(classification_file):
    """Load pairs from classification file."""
    pairs = []
    
    with open(classification_file) as f:
        next(f)  # Skip header
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            classification = parts[4]
            
            # Only process these classifications
            if classification not in ['C_contradictory', 'C_anomaly', 'C_no_result']:
                continue
            
            pairs.append({
                'query_protein': parts[0],
                'target_protein': parts[1],
                'query_gene': parts[2],
                'target_gene': parts[3],
                'classification': classification,
                'alt_gene': parts[5] if len(parts) > 5 else '-',
                'aligned_chr': parts[6] if len(parts) > 6 else '-',
                'aligned_start': parts[7] if len(parts) > 7 else '-',
                'aligned_end': parts[8] if len(parts) > 8 else '-',
                'overlap_info': parts[9] if len(parts) > 9 else '-'
            })
    
    return pairs


def write_results(output_file, pairs):
    """Write results to TSV file."""
    with open(output_file, 'w') as out:
        out.write('\t'.join([
            'Query_Protein', 'Target_Protein', 'Query_Gene', 'Target_Gene',
            'Classification', 'Alternative_Gene',
            'PID', 'Overlap', 'DepthNorm', 'CactusScore', 'AltCactusScore'
        ]) + '\n')
        
        for pair in pairs:
            # Format numbers
            pid = f"{pair['pid']:.4f}" if pair['pid'] is not None else '-'
            overlap = f"{pair['overlap']:.4f}" if pair['overlap'] is not None else '-'
            depth = f"{pair['depth_norm']:.4f}" if pair['depth_norm'] is not None else '-'
            score = f"{pair['cactus_score']:.4f}" if pair['cactus_score'] is not None else '-'
            alt_score = f"{pair['alt_cactus_score']:.4f}" if pair['alt_cactus_score'] is not None else '-'
            
            out.write('\t'.join([
                pair['query_protein'],
                pair['target_protein'],
                pair['query_gene'],
                pair['target_gene'],
                pair['classification'],
                pair['alt_gene'],
                pid, overlap, depth, score, alt_score
            ]) + '\n')


def process_species_pair(
    classification_file,
    alignment_dir,
    depth_dir,
    query_bed,
    target_bed,
    target_genome,
    output_file,
    species_pair,
    n_species=3
):
    """Process one species pair."""
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"Processing species pair: {species_pair}", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    
    # Convert to Path objects
    classification_file = Path(classification_file)
    alignment_dir = Path(alignment_dir)
    depth_dir = Path(depth_dir)
    query_bed = Path(query_bed)
    target_bed = Path(target_bed)
    output_file = Path(output_file)
    
    # Load gene lengths
    print(f"  Loading gene lengths...", file=sys.stderr)
    query_lengths = load_gene_lengths(query_bed)
    target_lengths = load_gene_lengths(target_bed)
    print(f"    Query genes: {len(query_lengths)}", file=sys.stderr)
    print(f"    Target genes: {len(target_lengths)}", file=sys.stderr)
    
    # Load pairs
    print(f"  Loading pairs...", file=sys.stderr)
    pairs = load_classification_pairs(classification_file)
    print(f"    Found {len(pairs)} pairs to process", file=sys.stderr)
    
    if not pairs:
        print(f"  No pairs to process. Skipping.", file=sys.stderr)
        return
    
    # Process each pair
    print(f"  Calculating Cactus scores...", file=sys.stderr)
    for i, pair in enumerate(pairs, 1):
        pair = process_pair(
            pair, alignment_dir, depth_dir,
            query_lengths, target_lengths,
            target_genome, n_species
        )
        
        if i % 50 == 0:
            print(f"    Processed {i}/{len(pairs)} pairs...", file=sys.stderr)
    
    # Write results
    print(f"  Writing results...", file=sys.stderr)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    write_results(output_file, pairs)
    
    # Summary statistics
    scores = [p['cactus_score'] for p in pairs if p['cactus_score'] is not None]
    alt_scores = [p['alt_cactus_score'] for p in pairs if p['alt_cactus_score'] is not None]
    
    print(f"\n  Summary:", file=sys.stderr)
    print(f"  {'-'*40}", file=sys.stderr)
    print(f"    Total pairs: {len(pairs)}", file=sys.stderr)
    print(f"    Scores calculated: {len(scores)}", file=sys.stderr)
    if scores:
        print(f"    Mean Cactus Score: {sum(scores)/len(scores):.4f}", file=sys.stderr)
    if alt_scores:
        print(f"    Mean Alt Score: {sum(alt_scores)/len(alt_scores):.4f}", file=sys.stderr)
    print(f"  Results saved to: {output_file}", file=sys.stderr)


def main():
    print("Starting Cactus Support Score calculation", file=sys.stderr)
    print(f"Timestamp: {__import__('datetime').datetime.now()}", file=sys.stderr)
    
    # Paths
    base_dir = "results/evaluation/mammalia/vs_prot_syn/nscore_approach"
    fp_dir = f"{base_dir}/disagreement_investigation/false_positives"
    depth_base = f"{base_dir}/disagreement_investigation/false_positives/alignment_depth"
    bed_dir = "results/gene_only_annotation/mammalia"
    output_dir = f"{base_dir}/disagreement_investigation/false_positives"
    
    # Number of species in alignment (excluding reference)
    n_species = 3  # Adjust based on your HAL file
    
    # can-mac
    process_species_pair(
        classification_file=f"{fp_dir}/alignment_classification_can_mac.tsv",
        alignment_dir=f"{fp_dir}/alignments_can_mac",
        depth_dir=f"{depth_base}/depths_can_mac",
        query_bed=f"{bed_dir}/Canis_genes_only_stableid_fix.bed6",
        target_bed=f"{bed_dir}/Macaca_genes_only_stableid_fix.bed6",
        target_genome="Macaca_fascicularis",
        output_file=f"{output_dir}/cactus_scores_can_mac.tsv",
        species_pair="can-mac",
        n_species=n_species
    )
    
    
    # can-rat
    process_species_pair(
        classification_file=f"{fp_dir}/alignment_classification_can_rat.tsv",
        alignment_dir=f"{fp_dir}/alignments_can_rat",
        depth_dir=f"{depth_base}/depths_can_rat",
        query_bed=f"{bed_dir}/Canis_genes_only_stableid_fix.bed6",
        target_bed=f"{bed_dir}/Rat_genes_only_stableid_fix.bed6",
        target_genome="Rattus_norvegicus",
        output_file=f"{output_dir}/cactus_scores_can_rat.tsv",
        species_pair="can-rat",
        n_species=n_species
    )
    
    # can-mus
    process_species_pair(
        classification_file=f"{fp_dir}/alignment_classification_can_mus.tsv",
        alignment_dir=f"{fp_dir}/alignments_can_mus",
        depth_dir=f"{depth_base}/depths_can_mus",
        query_bed=f"{bed_dir}/Canis_genes_only_stableid_fix.bed6",
        target_bed=f"{bed_dir}/Mus_genes_only_stableid_fix.bed6",
        target_genome="Mus_musculus",
        output_file=f"{output_dir}/cactus_scores_can_mus.tsv",
        species_pair="can-mus",
        n_species=n_species
    )
    
    # mac-rat
    process_species_pair(
        classification_file=f"{fp_dir}/alignment_classification_mac_rat.tsv",
        alignment_dir=f"{fp_dir}/alignments_mac_rat",
        depth_dir=f"{depth_base}/depths_mac_rat",
        query_bed=f"{bed_dir}/Macaca_genes_only_stableid_fix.bed6",
        target_bed=f"{bed_dir}/Rat_genes_only_stableid_fix.bed6",
        target_genome="Rattus_norvegicus",
        output_file=f"{output_dir}/cactus_scores_mac_rat.tsv",
        species_pair="mac-rat",
        n_species=n_species
    )
    
    # mac-mus
    process_species_pair(
        classification_file=f"{fp_dir}/alignment_classification_mac_mus.tsv",
        alignment_dir=f"{fp_dir}/alignments_mac_mus",
        depth_dir=f"{depth_base}/depths_mac_mus",
        query_bed=f"{bed_dir}/Macaca_genes_only_stableid_fix.bed6",
        target_bed=f"{bed_dir}/Mus_genes_only_stableid_fix.bed6",
        target_genome="Mus_musculus",
        output_file=f"{output_dir}/cactus_scores_mac_mus.tsv",
        species_pair="mac-mus",
        n_species=n_species
    )
    
    # rat-mus
    process_species_pair(
        classification_file=f"{fp_dir}/alignment_classification_rat_mus.tsv",
        alignment_dir=f"{fp_dir}/alignments_rat_mus",
        depth_dir=f"{depth_base}/depths_rat_mus",
        query_bed=f"{bed_dir}/Rat_genes_only_stableid_fix.bed6",
        target_bed=f"{bed_dir}/Mus_genes_only_stableid_fix.bed6",
        target_genome="Mus_musculus",
        output_file=f"{output_dir}/cactus_scores_rat_mus.tsv",
        species_pair="rat-mus",
        n_species=n_species
    )
    
    print("\n" + "="*60, file=sys.stderr)
    print("ALL SPECIES PAIRS COMPLETED", file=sys.stderr)
    print("="*60, file=sys.stderr)
    print(f"Finished: {__import__('datetime').datetime.now()}", file=sys.stderr)


if __name__ == '__main__':
    main()