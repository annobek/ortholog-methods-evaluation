EVALUATION OF METHODS FOR ORTHOLOG IDENTIFICATION
================================================

Author: Anna Beketova

This README contains information about evaluation of multiple 1:1 ortholog
identification methods using reference-based metrics across mammalian species.


# REFERENCE ORTHOLOG SETS

Following methods were used as reference:

  - Cactus-based approach for 1:1 ortholog identification
  - Union set of BUSCO and Compleasm 1:1 orthologs (BUSComp-set)
  - Union set of Cactus and BUSComp-set orthologs (CactBUSComp-set)


# EVALUATED METHODS

1:1 orthologs of following methods were evaluated against the references:

  - Prot-syn - N-score and Sum mods, no intermediate transitivity
  - Tree-algorithm - majority, max_score (standard) and whitelist mods
  - PhyloTreePruner (PTP)
  - OrthoFinder2

Methods were also evaluated against each other (Jaccard Similarity only).

Evaluation metrics:
  Precision, Recall, F1-score, Jaccard Similarity

For the evaluation of prot-syn vs cactus, False Negatives (FN) and False Positives (FP)
were analysed in details.


# EVALUATION OUTPUTS

Two heatmaps were generated:

  1) Methods vs reference sets
     results/evaluation/mammalia/metrics_methods_vs_reference_sets_v2.png

  2) Jaccard similarity: methods vs each other and reference sets
     results/evaluation/mammalia/all_vs_all/combined_jaccard_heatmap_v3.png

Plots were copied to:
  gene_relationship_classifier/plots/mammals/evaluation/


# PROJECT NAVIGATION

Only directories used for the evaluation are shown.
```
cactus/
|-- cactus-bin-v2.9.2/ - virtual environment directory for halTools
|-- venv-analysis/ - virtual environment directory for Python scripts
|-- material/
|   |-- sex_experiment/ - some files needed for evaluation
|-- gene_relationship_classifier/
|   |-- methods/ - scripts needed for the evaluation (except for BUSCO vs compleasm and FP-scripts)
|   |-- plots/
|       |-- mammals/
|           |-- evaluation/ - plots displayed in the gene_relationship_classifier repository
|-- methods/ - scripts needed for FP analysis
|-- results/
|   |-- evaluation/
|       |-- mammalia/
|           |-- all_vs_all/ - jaccard similarity of methods vs each other
|           |-- all_vs_cactbuscomp/ - methods vs CactBUSComp-set
|           |   |-- vs_orthofinder/
|           |   |-- vs_prot_syn/
|           |   |   |-- nscore_approach/
|           |   |   |   |-- test_no_transitivity/
|           |   |   |-- sum_approach/
|           |   |       |-- test_no_transitivity/
|           |   |-- vs_ptp/
|           |   |-- vs_tree/
|           |       |-- majority/
|           |       |-- max_score/
|           |       |-- whitelist/
|           |-- vs_buscomp_set/ - cactus vs BUSComp-set and methods vs BUSComp-set
|           |   |-- vs_orthofinder/
|           |   |-- vs_prot_syn/
|           |   |   |-- nscore_approach/
|           |   |   |   |-- test_no_transitivity/
|           |   |   |-- sum_approach/
|           |   |       |-- test_no_transitivity/
|           |   |-- vs_ptp/
|           |   |-- vs_tree/
|           |       |-- majority/
|           |       |-- max_score/
|           |       |-- whitelist/
|           |-- vs_orthofinder/ - cactus vs OrthoFinder2
|           |-- vs_ptp/ - cactus vs PTP
|           |-- vs_prot_syn/ - cactus vs prot-syn
|           |   |-- nscore_approach/
|           |   |   |-- test_no_transitivity/ - metrics, FN and FP
|           |   |       |-- disagreement_investigation/ - detailed analysis of FN and FP
|           |   |           |-- false_negatives/
|           |   |           |-- false_positives/
|           |   |-- sum_approach/
|           |       |-- test_no_transitivity/ - metrics, FN and FP
|           |           |-- disagreement_investigation/ - detailed analysis of FN and FP
|           |               |-- false_negatives/
|           |               |-- false_positives/
|           |-- vs_tree/ - cactus vs tree
|               |-- majority/
|               |-- max_score/
|               |-- whitelist/
|-- slurm_outputs/ - SLURM output and error scripts
```

# DATASET

Following species were analyzed:

  - Canis lupus familiaris (can)
  - Macaca fascicularis (mac)
  - Rattus norvegicus (rat)
  - Mus musculus (mus)

Analyzed pairs:

  - can-mac
  - can-rat
  - can-mus
  - mac-rat
  - mac-mus
  - rat-mus


# SCRIPTS AND EXECUTION ENVIRONMENT

Scripts (except for BUSCO vs compleasm and for FP investigation) are located in:
  - `gene_relationship_classifier/methods/`

For BUSCO vs compleasm:
  - `/storage/EasyVectorOmics/busco/methods/evaluate_busco_vs_compleasm_v2.py`
  - `/storage/EasyVectorOmics/busco/methods/evaluate_busco_compleasm.sh`

For FP investigation:
  - `methods/`

The analysis was done on the Bioserver (Linux Ubuntu server). For each step,
a corresponding bash script was executed with SLURM.

For all Python scripts virtual environment for Python was used:
  - `source venv-analysis/bin/activate`

For FP investigation virtual environment for HALTools was used:
  - `source cactus-bin-v2.9.2/venv-cactus-v2.9.2/bin/activate`


# FILES USED

Cactus 1:1 orthologs (dog-macaca pair example):
  - `results/reciprocal_pairs/can_mac_one2one_fix.tsv`

BUSCO (dog-macaca pair example):
  - `/storage/EasyVectorOmics/busco/results/mammalia/busco_proteome/mammalia_lineage/
  mapped_1to1_busco_can_mac.tsv`

Compleasm (dog-macaca pair example):
  - `/storage/EasyVectorOmics/busco/results/mammalia/compleasm_proteome/longest_isoform/mapped_1to1_compleasm_can_mac.tsv`

BUSComp-set:
  - `/storage/EasyVectorOmics/busco/results/mammalia/evaluation/BUSComp_set.tsv`

CactBUSComp-set:
  - `results/evaluation/mammalia/vs_buscomp_set/CactBUSComp_set.tsv`

Map geneID <-> ProteinID:
  - `material/sex_experiment/gene_to_protein_map_FIXED.tsv`

Prot-syn (nscore approach):
  - `/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_nscore.tsv`

Prot-syn (sum approach):
  - `/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_sum.tsv`

Tree - majority:
  - `/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/majority/
  complete_tree_conserved_orthologs_2.tsv`

Tree - max_score:
  - `/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/standard/
  complete_tree_conserved_orthologs_2.tsv`

Tree - whitelist:
  - `/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/whitelist/
  complete_tree_conserved_orthologs_2.tsv`

PTP:
  - `/storage/EasyVectorOmics/phylotreepruner/results/mammalia/ptp_one_to_one_pairs_new.tsv`

OrthoFinder:
  - `results/evaluation/mammalia/vs_orthofinder/orthofinder_1to1.tsv`

Dictionaries from Prot-syn:
  - `/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/`

HAL graph (Cactus raw output):
  - `results/output_alignment_sexspecies.hal`

# EVALUATION STEPS

## 1) BUSCO vs Compleasm + BUSComp-set creation

Task:
  To evaluate BUSCO (reference) against Compleasm (test) and create BUSComp-set

Script:
  
  `/storage/EasyVectorOmics/busco/methods/evaluate_busco_vs_compleasm.sh`
  (executes evaluate_busco_vs_compleasm_v2.py)

Input:
  BUSCO and Compleasm 1:1 ortholog files and map geneID<->proteinID.

Output location:
  `/storage/EasyVectorOmics/busco/results/mammalia/evaluation/`

Creates:
  - busco_compleasm_metrics.tsv (and .txt) - evaluation metrics
  - all_false_negatives.tsv - pairs that are found in reference, but not found in test
    (all species pairs)
  - all_false_positives.tsv - pairs that are found in test, but not found in reference
    (all species pairs)
  - per-species-pair false negatives and false positives
  - BUSComp_set.tsv (and .pkl) - union set of BUSCO and compleasm orthologs


## 2) BUSComp-set vs Cactus + CactBUSComp-set creation


Task:
  To evaluate BUSComp-set (reference) against Cactus (test)

Script:
  
  `evaluate_cactus_vs_buscomp.sh`
  (executes evaluate_cactus_vs_buscomp_v2.py)

Input:
  Cactus 1:1 ortholog files, map geneID<->proteinID, BUSComp-set.

Output location:
  `results/evaluation/mammalia/vs_buscomp_set/`

Creates:
  - cactus_BUSComp_metrics.tsv (and .txt) - evaluation metrics
  - all_false_negatives_cactus.tsv - pairs that are found in reference, but not found
    in test (all species pairs)
  - all_false_positives_cactus.tsv - pairs that are found in test, but not found
    in reference (all species pairs)
  - per-species-pair false negatives and false positives
  - CactBUSComp_set.tsv (and .pkl) - union set of Cactus and BUSComp-set orthologs


## 3) Methods vs Cactus 

Task:
  To evaluate methods (test) against cactus (reference)

Script:
  
  `evaluate_methods_vs_cactus.sh` 
  (executes evaluate_methods_vs_cactus_v2.py)

Input:
  cactus 1:1 ortholog files, map geneID<->proteinID, and 1:1 ortholog files of the tested
  methods.

Output location:
  `results/evaluation/mammalia/` (in corresponding methods directories)

Creates:
  - cactus_<method>_metrics.tsv (and .txt) - Statistical data and evaluation metrics:
    Precision, Recall, F1-score and Jaccard Similarity
  - all_false_negatives_<method>.tsv - pairs that are found in reference, but not found
    in test (all species pairs)
  - all_false_positives_<method>.tsv - pairs that are found in test, but not found
    in reference (all species pairs)

Notes:
  See the script for exact paths.


## 4) Methods vs BUSComp-set

Task:
  To evaluate methods (test) against BUSComp-set (reference)

Script:
  
  `evaluate_methods_vs_buscomp.sh`
  (executes evaluate_methods_vs_buscomp_v2.py)

Input:
  BUSComp-set file, map geneID<->proteinID, and 1:1 ortholog files of the tested methods.

Output location:
  `results/evaluation/mammalia/vs_buscomp_set/` (in corresponding methods directories)

Creates:
  - BUSComp_<method>_metrics.tsv (and .txt) - Statistical data and evaluation metrics:
    Precision, Recall, F1-score and Jaccard Similarity
  - all_false_negatives_<method>.tsv - pairs that are found in reference, but not found
    in test (all species pairs)
  - all_false_positives_<method>.tsv - pairs that are found in test, but not found
    in reference (all species pairs)
  - per-species-pair false negatives and false positives

Notes:
  See the script for exact paths.


## 5) Methods vs CactBUSComp-set

Task:
  To evaluate methods (test) against CactBUSComp-set (reference)

Script:
  
  `evaluate_methods_vs_cactbuscomp.sh`
  (executes evaluate_methods_vs_cactbuscomp_v2.py)

Input:
  CactBUSComp-set file, map geneID<->proteinID, and 1:1 ortholog files of the tested
  methods.

Output location:
  `results/evaluation/mammalia/all_vs_cactbuscomp/` (in corresponding methods directories)

Creates:
  - CactBUSComp_<method>_metrics.tsv (and .txt) - Statistical data and evaluation metrics:
    Precision, Recall, F1-score and Jaccard Similarity
  - all_false_negatives_<method>.tsv - pairs that are found in reference, but not found
    in test (all species pairs)
  - all_false_positives_<method>.tsv - pairs that are found in test, but not found
    in reference (all species pairs)
  - per-species-pair false negatives and false positives

Notes:
  See the script for exact paths.


## 6) Methods vs each other + Jaccard similarity heatmap

Task:
  To evaluate methods against each other and plot Jaccard similarity

Script:
  
  `evaluate_methods_vs_each_other_and_plot.sh`
  (executes evaluate_methods_vs_each_other_and_plot_v3.py)

Input:
  Ortholog files of methods Prot-syn, Tree, PTP and OrthoFinder, and already
  evaluated TSV files with metrics.

Output location:
  `results/evaluation/mammalia/all_vs_all/`

Creates:
  - combined_jaccard_heatmap_v3.png - heatmap of Jaccard similarities
  - combined_jaccard_matrix_v3.txt - Jaccard as text file


## 7) Heatmap: Methods vs reference sets

Task:
  To visualise evaluation of methods against reference sets (Cactus, CactBUSComp and
  BUSComp) as heatmap

Script:
  
  `visualize_evaluation.sh`
  (executes visualize_evaluation_metrics_v2.py)

Input:
  TSV with metrics for each method-reference pair

Output location:
  `results/evaluation/mammalia/`

Creates:
  - metrics_methods_vs_reference_sets_v2.png


# DISAGREEMENT INVESTIGATION

Disagreements between Cactus and prot-syn ortholog sets (both nscore and sum
approaches) were investigated in details.
Here examples for nscore approach are shown, for sum approach see the directory
sum_approach/ in results/evaluation/mammalia/vs_prot_syn/


## A) INVESTIGATION OF FN

### A1) neighborhood_score_cactus.sh (executes neighborhood_score_cactus_fn_v4.py)

Analyzes FN ortholog pairs (Pairs found by Cactus, and missed by prot-syn) by
computing neighborhood scores (N-scores, number of shared neighbors between genes) for the Cactus pair and (when available)
the prot-syn assigned alternative, classifies the FN into prot-syn_no_result and
contradictory (Cactus and prot-syn find different candidates to a gene).

The analysis is at both gene-level (two rows per FN pair, each row for the gene from
the FN pair and it's partner gene assigned by Cactus and prot-syn, if exists) and pair-level
(one row per FN with mean/max N-score for multiple ProtSyn candidates for a FN pair and score
fractions).

Requires:
  - neighborhood dictionary (neighborhoods.pkl)
  - relationship dictionary (geometric_mean.pkl)
  - geneID<->proteinID mapping TSV
  - prot-syn 1:1 ortholog TSV
  - FN pairs TSV

Input dictionaries:
  `/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/`

Prot-syn orthologs:
  `/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/
  one_to_one_nscore.tsv and one_to_one_sum.tsv`

FN table:
  `results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/
  all_false_negatives_nscore.tsv`

Outputs:
  `results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/
  disagreement_investigation/false_negatives/`

Outputs include:
  - false_negatives_with_ids_pair_level.tsv - FN pairs with assigned FN_IDs for each pair
  - false_negatives_classified_gene_level.tsv - classified gene-level table (one gene per row with ID of its FN pair, 
    N-scores, fractions, relative-score columns)
  - false_negatives_edge_affected_gene_rows.tsv - genes affected by max n-score < 20
    ("edge-affected" genes)
  - false_negatives_classified_pair_level.tsv - classified pair-level table (one FN pair per row, N-scores, fractions)
  - fn_distributions_2x2_gene_level.png and fn_distributions_2x2_pair_level.png -
    two 2x2 distribution plots (gene-level and pair-level)
  - contradictory_cactus_higher_nscore_than_protsyn_gene_level.tsv and
    anomalous_protsyn_nscore_zero_gene_level.tsv - two gene-level subsets
    (contradictory cases, where cactus chose a candidate with better N-score and anomalous ProtSyn score=0 cases)


### A2) check_in_excluded_fn.sh (executes check_in_excluded_fn.py)

Checks whether FN gene IDs appear in species-specific gene lists
(dog, rat, mouse) excluded by prot-syn.

Each row in the FN table represents one gene from an FN pair.

Requires:
  - FN table
  - excluded-gene tables located in `material/sex_experiment/`
    excluded_genes_<species>.tsv (Provided files for Canis, Rat and Mouse;
    no excluded genes for Macaca)

Outputs:
  `results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/
  disagreement_investigation/false_negatives/`

Outputs include:

  - `fn_genes_found_in_excluded.tsv`
    - Subset of FN rows where Focus_Gene is found in the excluded list
    - Includes all original FN columns + helper columns

  - `fn_genes_found_in_excluded_summary.tsv`
    - Per-species summary: species_key | total_fn_rows | excluded_hits
    - NOTE: contains Unknown species (Macaca), because no Macaca excluded-genes table
      was provided


## B) INVESTIGATION OF FP

### B1) find_c_no_results_fp.sh (executes find_c_no_results_fp.py)

Among FP pairs finds those that have no result in Cactus, and saves them for further
investigation.

Logic:
  - Check Gene1: if found anywhere in cactus -> ignore the whole FP pair
  - If Gene1 not found -> check Gene2: if found anywhere -> ignore
  - Only if BOTH genes are not found anywhere -> C_no_result

Input:
  - all_false_positives_nscore.tsv - False positives table:
    `results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/`

  - can_mac_one2one_fix.tsv - Cactus ortholog file (Dog-Macaca example):
    `results/reciprocal_pairs/`

Output:
  `all_C_no_results.tsv` - FP pairs where both genes have no Cactus result:
  
  - `results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/
  disagreement_investigation/false_positives/`


### B2) build_pair_coordinate_table.sh (located in methods/)

From a TSV file of gene pairs, extract:
  - Gene1, Species1
  - Gene2_prot_syn, Species2_prot_syn

Then, for each gene and species:
  - look up genomic coordinates in the corresponding species GTF
  - convert them to HAL-friendly coordinates (0-based start, half-open end)

Input TSV requirements (column headers):
  - Gene1
  - Species1
  - Gene2_prot_syn
  - Species2_prot_syn

GTF files:

  Location:
    `/storage/EasyVectorOmics/cactus/material/sex_experiment`

  Names:
  - gtf_Canis.gtf
  - gtf_Macaca.gtf
  - gtf_Mus.gtf
  - gtf_Rat.gtf

Output TSV columns:

  `gene1 gene2 species1 species2 chr1 chr2 start1 start2 end1 end2 length1 length2`

Coordinate conversion:
  - GTF: 1-based inclusive [start1..end1]
  - HAL: start0 = start1-1, length = end1-start1+1, end0 = start0+length


### B3) get_cactus_score.sh (located in methods/)

For every pair performs reciprocal liftover to estimate aligned overlap, measures
sequence identity and alignment depth, and then combines them into:
  ```
  Score = (PID × overlap × DepthNorm)^(1/3)
  (all coordinates are expected to be 0-based, half-open)
  ```

Tools used:
  - HAL tools: halLiftover, halSnps, halAlignmentDepth
  - BED tools: bedtools intersect, bedtools sort, bedtools merge
  - Python 3 + awk for calculations and summarization

What is calculated (per pair):
  - aligned_A, aligned_B: aligned bases after reciprocal liftover overlap with the
    gene interval
  - o (overlap): `(aligned_A + aligned_B) / (len1 + len2)`
  - PID: from halSnps summary (identity based on mismatches vs totalClean)
  - MeanDepth: mean alignment depth over the gene interval from halAlignmentDepth
  - DepthNorm: `(MeanDepth - 1) / (N_SPECIES - 1)` (no clamping in the current script)
  - Score: `cube root of PID * o * DepthNorm` (signed cube root)

Input:
  - Pairs with coordinates (TSV):

    `/storage/EasyVectorOmics/cactus/results/evaluation/mammalia/vs_prot_syn/
    nscore_approach/test_no_transitivity/disagreement_investigation/false_positives/
    pairs_with_coords.tsv`

  - HAL alignment:

    `/storage/EasyVectorOmics/cactus/results/output_alignment_sexspecies.hal`

Output:
  - Pairs with computed Cactus score (TSV):
  
    `/storage/EasyVectorOmics/cactus/results/evaluation/mammalia/vs_prot_syn/
    nscore_approach/test_no_transitivity/disagreement_investigation/false_positives/
    pairs_with_cactus_score.tsv`


# REFERENCES

- Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov,
  BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper
  Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes,
  Molecular Biology and Evolution, Volume 38, Issue 10, October 2021, Pages 4647–4654,
  https://doi.org/10.1093/molbev/msab199

- Neng Huang, Heng Li, compleasm: a faster and more accurate reimplementation of BUSCO.
  Bioinformatics, 39, btad595, 2023. doi:10.1093/bioinformatics/btad595

- Glenn Hickey, Benedict Paten, Dent Earl, Daniel Zerbino, and David Haussler.
  HAL: A Hierarchical Format for Storing and Analyzing Multiple Genome Alignments.
  Bioinformatics. 2013. https://doi.org/10.1093/bioinformatics/btt128

- Armstrong, J., Hickey, G., Diekhans, M. et al. Progressive Cactus is a multiple-genome
  aligner for the thousand-genome era. Nature 587, 246–251 (2020).
  https://doi.org/10.1038/s41586-020-2871-y

- Kocot KM, Citarella MR, Moroz LL, Halanych KM. PhyloTreePruner: A Phylogenetic Tree-Based
  Approach for Selection of Orthologous Sequences for Phylogenomics. Evol Bioinform Online.
  2013 Oct 29;9:429-35. doi: 10.4137/EBO.S12813. PMID: 24250218; PMCID: PMC3825643.

- David M Emms, Yi Liu, Laurence Belcher, Jonathan Holmes, Steven Kelly. OrthoFinder:
  scalable phylogenetic orthology inference for comparative genomics. bioRxiv 2025.07.15.664860;
  doi: https://doi.org/10.1101/2025.07.15.664860

- Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic
  features. Bioinformatics. 26, 6, pp. 841–842.
