Comparative Evaluation of One-to-One Ortholog Identification Methods
=====================================================================

Bachelor thesis of Anna Beketova.

This repository documents the workflows and scripts used to identify and evaluate 1:1 orthologs across multiple mammalian species using different approaches. The primary goal is to compare different ortholog identification methods against multiple reference sets (Cactus, BUSComp-set, and CactBUSComp-set) using standard evaluation metrics.

## Project navigation

* Cactus workflow
* BUSCO/Compleasm workflow
* PTP workflow
* OrthoFinder2 workflow
* Evaluation workflow

All paths in this README are in relation to the Bioserver. The repository contains only the scripts and plots that were used for the thesis. 

In case some scripts are missing on the repository they are to be found on the server:

- Cactus-related scripts, evaluation and 1:1 ortholog extraction for OrthoFinder: `/storage/EasyVectorOmics/cactus/gene_relationship_classifier/methods/`

- BUSCO and Compleasm: `/storage/EasyVectorOmics/busco/methods/`

- PhyloTreePruner: `/storage/EasyVectorOmics/phylotreepruner/methods/`

In case there are some mammals-related paths left commented, uncomment them before reproducing the results.
Change the input/output paths if needed. 

## Analysed species and species pairs

Provided mammalian data were obtained from this study:  
Naqvi, S. et al. (2019). *Conservation, acquisition, and functional impact of sex-biased gene expression in mammals*. **Science**, 365(6450).  
[https://doi.org/10.1126/science.aaw7317](https://doi.org/10.1126/science.aaw7317)


- Canis lupus familiaris: can
- Macaca fascicularis: mac
- Rattus norvegicus: rat
- Mus musculus: mus

Analyzed pairs: can-mac, can-rat, can-mus, mac-rat, mac-mus, rat-mus


## Software used

### Python environment

- **Python** v3.12.3  
  with the following packages:
  - biopython==1.85
  - contourpy==1.3.2
  - cycler==0.12.1
  - fonttools==4.58.1
  - kiwisolver==1.4.8
  - matplotlib==3.10.3
  - matplotlib-venn==1.1.2
  - natsort==8.4.0
  - ncls==0.0.68
  - networkx==3.5
  - numpy==2.2.6
  - packaging==25.0
  - pandas==2.2.3
  - pillow==11.2.1
  - ptitprince==0.3.1
  - pyparsing==3.2.3
  - pyranges==0.1.4
  - python-dateutil==2.9.0.post0
  - pytz==2025.2
  - scipy==1.16.3
  - seaborn==0.13.2
  - setuptools==80.9.0
  - six==1.17.0
  - sorted_nearest==0.0.39
  - tabulate==0.9.0
  - tzdata==2025.2
  - wheel==0.45.1

### Tools used for analysis

- Cactus v7.0.0
- halSynteny v2.2
- bedtools intersect v2.30.0
- convert2bed (bedops v2.4.39)
- DIAMOND v2.1.9
- Compleasm v0.2.6
- BUSCO v5.8.3
- gffread v0.12.7
- AGAT v1.4.3
- OrthoFinder (MSA mode) v2.5.5
- PhyloTreePruner v1.0

### Additional tools for virtual environment management

- conda v25.5.1
- mamba v2.3.0



All analysis was done on a Linux Ubuntu server, for each step corresponding bash script was executed with SLURM.

The scripts were written with the help of AI for robustness, debugging and readability. 

## Preparing Cactus-based reference set

This reference-free and synteny-based pipeline performs pairwise analysis of protein-coding genes in given species pairs using alignments done with Cactus and synteny identification with halSynteny, and extracts 1:1 orthologs (for each source gene exactly one target candidate) using custom Python scripts.

There are optional scripts located in `methods/data_check/` that can be used for debugging and additional data checks. They are not required for further evaluation steps. 

On the server the scripts are located in `/storage/EasyVectorOmics/cactus/gene_relationship_classifier/methods/`

The results are located in `/storage/EasyVectorOmics/cactus/results/`

Note: in the scripts GRASS is called GRIS or geometric mean

### Material

- `/storage/EasyVectorOmics/cactus/results/output_alignment_sexspecies.hal` - Cactus alignment
- `/storage/EasyVectorOmics/cactus/material/sex_experiment/gtf_* ` - genome annotations

### Workflow (see detailed instructions below)

1. `process_hal.sh` – Runs halSynteny on cactus results  
2. `psl_to_bed.sh` – Converts halSynteny results to BED  
3. Optional: `data_check/check_regions_after_halsynteny.sh` – Prints statistics about number and length of synteny regions  
4. `filter_annotation_to_gene_only.sh` – Filters genome annotations keeping only gene entries  
5. `create_map_gene_names_id_mammalia.sh` (for mammals) – Creates a map of gene symbols <-> unique gene IDs  
6. `gtf_to_bed_mammalia.sh` – Converts annotations from GTF/GFF to BED  
7. `map_gene_names_to_ids_mammalia.sh` – Maps gene symbols from the BED files to unique IDs  
8. `identify_genes_in_regions.sh` – Runs `bedtools intersect` to identify genes in synteny regions  
9. `create_region_and_gene_tables.sh` – Creates tables of synteny regions and genes located in those regions  
10. Optional: `data_check/plot_genes_per_region.sh` (executes `data_check/plot_genes_per_region.py`) – Plots distribution of number of genes per synteny region  
11. Optional: `data_check/plot_total_number_genes_mammals.sh` (executes `data_check/plot_total_number_genes_mammals.py`) – Plots pie charts showing how much of the total number of genes lies within synteny regions  
12. Optional: `data_check/plot_total_coverage.sh` (executes `data_check/plot_total_coverage.py`) – Plots percentage of the whole genome within synteny regions  
13. Optional: `data_check/plot_region_conservativeness.sh` (executes `data_check/plot_region_conservativeness.py`) – Plots percentage of how much of the synteny between two genomes is reciprocally conserved  
14. `run_linear_synteny.sh` (executes `run_linear_synteny.py` and `synteny_with_legend.py`) – Creates linear synteny plots  
15. `calculate_region_gris.sh` (executes `calculate_region_gris.py`) – Calculates region GRIS score  
16. `reciprocity_and_resolve_one_to_many.sh` (executes `reciprocity_and_resolve_one_to_many.py`) – Keeps only those pairs that are present in forward and reverse, extracts 1:1 pairs and resolves 1:many cases, placing duplicates in a separate file  
17. Optional: `data_check/check_one_to_one.sh` (executes `data_check/check_one_to_one.py`) – Ensures the extracted 1:1 pairs are truly 1:1



### Instructions

Before running the pipeline, make sure to install all required software.
Prepare a `material/` directory with genome annotation files (GTF or GFF) and
nucleotide and amino acid FASTA files.

This analysis was performed on the server using the following virtual environments:

**venv for Cactus and halTools:**

```bash
source /storage/EasyVectorOmics/cactus/cactus-bin-v2.9.2/venv-cactus-v2.9.2/bin/activate
```

**bedtools venv (for `bedtools intersect`):**

```bash
eval "$(/storage/EasyVectorOmics/busco/miniconda3/bin/conda shell.bash hook)"
eval "$(mamba shell hook --shell bash)"
mamba activate bedtools-env
```

**bedops venv (for `convert2bed`):**

```bash
eval "$(/storage/EasyVectorOmics/busco/miniconda3/bin/conda shell.bash hook)"
eval "$(mamba shell hook --shell bash)"
mamba activate bedops-env
```

**venv for Python scripts:**

```bash
source /storage/EasyVectorOmics/cactus/venv-analysis/bin/activate
```

The scripts were executed using SLURM.




### Detailed workflow

1. **halSynteny:**
   halTools venv required.
   To execute halSynteny run the script `process_hal.sh`; this will create the directory
   `results/halsynteny_output/mammalia/` and `.psl` files inside.

2. **Convert PSL files to BED:**
   Run the script `psl_to_bed.sh`; this will create the directory
   `results/halsynteny_output/mammalia/bed/` and `.bed` files inside.

3. **Optional: Check number of synteny regions and their length:**
   Run the script `data_check/check_regions_after_halsynteny.sh`; this prints statistics to stdout and creates the directory
   `results/data_check/block_length/` with distributions of synteny region lengths.
   The latest results are in SLURM output:
   `/storage/EasyVectorOmics/cactus/slurm_outputs/block_length825.out`

4. **Filter genome annotation to gene-only:**
   Run the script `filter_annotation_to_gene_only.sh`; this will create the directory
   `results/gene_only_annotation/mammalia/` and `*_genes_only_fix.gtf` files inside.

5. **Create a map gene name ↔ unique gene ID:**
   Run the script `create_map_gene_names_id_mammalia.sh`; this will create
   `results/gene_only_annotation/mammalia/*_gene_to_id_fix.tsv` files.

6. **Convert GTF files to BED with bedops:**
   bedops venv required.
   Run the script `gtf_to_bed_mammalia.sh`; this will create
   `results/gene_only_annotation/mammalia/*.bed6` files.

7. **Map gene symbols from BED files to unique IDs:**
   Run the script `map_gene_names_to_ids_mammalia.sh`; this will create
   `results/gene_only_annotation/mammalia/*_stableid_fix.bed6` files.

8. **Identify genes in synteny regions with bedtools intersect:**
   bedtools venv required.
   Run the script `identify_genes_in_regions.sh`; this will create directories under
   `results/genes_in_blocks/mammalia/` and `*query_fix.bed` / `*target_fix.bed` files for query and target species synteny regions.

   The latest results are in SLURM output:
   `/storage/EasyVectorOmics/cactus/slurm_outputs/find_genes_in_blocks823.out`

9. **Create tables of synteny regions and genes located in those regions:**
   Run the script `create_region_and_gene_tables.sh`; this will create
   `results/genes_in_blocks/mammalia/region_gene_tables/` and the files
   `*_region_fix.tsv` and `*_genes_fix.tsv`.

10. **Optional: Plot distribution of number of genes per synteny region:**
    Python venv required.
    Run the script `data_check/plot_genes_per_region.sh` (executes `data_check/plot_genes_per_region.py`).
    This creates the following files in `results/data_check/mammalia/gene_distribution_in_blocks/`:

    * `*_number_genes_in_blocks.tsv`
    * `*_length_vs_gene_number_query.png`
    * `*_length_vs_gene_number_target.png`
    * `*_gene_distribution_overview.png`
    * `*_gene_distribution_0_100.png`

11. **Optional: Plot pie charts showing how much of the total number of genes lies within synteny regions:**
    Python venv required.
    Run the script `data_check/plot_total_number_genes_mammalia.sh` (executes `data_check/plot_total_number_genes_mammalia.py`).
    This creates files in `results/data_check/mammalia/total_number_genes/`:

    Per-pair pie charts:

    * `*piechart_{source}_fix.png`
    * `*piechart_{target}_fix.png`

    Summary files:

    * `summary_points_fix.tsv`
    * `summary_boxplot_fix.png`

12. **Optional: Plot percentage of whole genome within synteny regions:**
    Python venv required.
    Run the script `data_check/plot_total_coverage.sh` (executes `data_check/plot_total_coverage.py`).
    This creates files in `results/data_check/mammalia/coverage/whole_genome/`:

    Per-pair bar plots:

    * `{query}_{target}_whole_genome_coverage_bar_fix.png`

    Summary files:

    * `summary_points_fix.tsv`
    * `summary_boxplot_fix.png`

13. **Optional: Plot percentage of reciprocal synteny conservation:**
    Python venv required.
    Run the script `data_check/plot_region_conservativeness.sh` (executes `data_check/plot_region_conservativeness.py`).
    The threshold of overlap can be adjusted if needed.
    Results are written to `results/data_check/mammalia/conservativeness/`, with directories `overlap_{threshold}/`.

    Per-pair bar plots:

    * `{query}_{target}_conservativeness_bar_overlap{threshold}_fix.png`

    Summary files:

    * `summary_points_fix.tsv`
    * `summary_boxplot_fix.png`

14. **Visualize synteny regions across chromosomes and scaffolds:**
    Python venv required.
    Run the script `run_linear_synteny.sh` (executes `run_linear_synteny.py`).
    The script `run_linear_synteny.py` defines input paths, species pairs, and filtering parameters, and calls `synteny_with_legend.py` for plotting.

    Parameters used:

    ```
    TOP_SCAFFOLDS = 10
    MIN_REGION_LENGTH = 2000000
    ```

    Input files (TSV):

    ```
    results/genes_in_blocks/mammalia/region_gene_tables/*_synteny_tables_regions_fix.tsv
    ```

    Output:

    ```
    results/linear_synteny_plots/linear_synteny_plot_forward.png
    ```

15. **Calculate region GRIS:**
    Python venv required.
    Run the script `calculate_region_gris.sh` (executes `calculate_region_gris.py`); this creates
    `results/region_gris/` and `*_region_gris_mean_fix.tsv`.

16. **Extract reciprocal pairs and resolve 1:many cases:**
    Python venv required.
    Run the script `reciprocity_and_resolve_one_to_many.sh` (executes `reciprocity_and_resolve_one_to_many.py`); this creates
    `results/reciprocal_pairs/` with:

    * `*_one2one_fix.tsv`
    * `*_inparalogs_fix.tsv`

17. **Check 1:1 relationships of extracted pairs:**
    Python venv required.
    Run the script `data_check/check_one_to_one.sh` (executes `data_check/check_one_to_one.py`); results are printed to stdout.


### Results navigation

After the analysis the results will have following structure:

```
results/
  |-- data_check/   <- plots and additional statistics are here
  |-- gene_only_annotation/ <- gene-only annotation
  |-- genes_in_blocks/ <- intermediate results: genes identified in synteny blocks           
  |-- halsynteny_output/ <- intermediate results: raw halsynteny output
  |-- reciprocal_pairs/   <- 1:1 ortholog pairs are here
  |-- region_gris/ <- interediate results: region GRASS
```

## BUSCO and Compleasm

This step prepares BUSCO and Compleasm 1:1 orthologs for another reference set.

**Note**:
All BUSCO and Compleasm preparation and evaluation steps are organized under the
`busco/` directory on the server. This includes extraction of 1:1 orthologs,
BUSCO–Compleasm comparison, and creation of the BUSComp-set.
These steps are separate from the Cactus-based workflow and from the main
evaluation directories under `results/evaluation/`.


### Material

- `/storage/EasyVectorOmics/cactus/material/sex_experiment/gtf_*` - genome annotations
- `/storage/EasyVectorOmics/cactus/material/sex_experiment/nucleotids/*.fna ` - nucleotide FASTA files




### INSTRUCTIONS

**Virtual environments**

**For AGAT:**

```bash
eval "$(/storage/EasyVectorOmics/busco/miniconda3/bin/conda shell.bash hook)"
eval "$(mamba shell hook --shell bash)"
mamba activate agat_env
```

**For gffread and BUSCO**
*(not for Compleasm, although the environment name suggests otherwise; for Compleasm use `compleasm026_env`):*

```bash
eval "$(/storage/EasyVectorOmics/busco/miniconda3/bin/conda shell.bash hook)"
eval "$(mamba shell hook --shell bash)"
mamba activate compleasm_env
```

**For Compleasm:**

```bash
eval "$(/storage/EasyVectorOmics/busco/miniconda3/bin/conda shell.bash hook)"
eval "$(mamba shell hook --shell bash)"
mamba activate compleasm026_env
```

**For Python scripts:**

```bash
source /storage/EasyVectorOmics/cactus/venv-analysis/bin/activate
```

On the server, the scripts are located in
`/storage/EasyVectorOmics/busco/methods/`,
and during the analysis they were executed relative to that location.
In the repository, they are stored in the `buscomp_scripts/` directory for display purposes only.

The results are located in:
`/storage/EasyVectorOmics/busco/results/`


#### To filter GTF files with AGAT, keeping the longest isoform

1. Start the virtual environment for AGAT.

2. Go to `methods/` directory

3. Run the script `agat.sh` with SLURM to filter GTF files and keep longest isoform representatives only. 
This will create `*longest_iso.gtf` files in `results/` directory.

#### To extract proteomes with gffread based on GTF file and FASTA

1. Start the virtual environment for gffread.

2. Go to the `methods/` directory

3. Run the script `gffread.sh` with SLURM to create proteomes based on the GTFs and genome FASTA files. 
This will create files `*proteome_longest.fa` in `results/` directory.


#### To run compleasm

1. Start the virtual environment for compleasm.

2. Create `canis/`, `macaca/`, `rat/` and `mus/` directories in `results/compleasm_proteome/`.

3. Go to the `methods/` directory

4. Run the script `download_lineages.sh` with SLURM to download `mammalia_odb10` lineage for compleasm.
All files will be saved in `material/compleasm_lineage/lineages/new_odb10/mammalia_odb10/` directory.

5. Run the script `compleasm.sh` with SLURM to perform compleasm.
This will generate files in corresponding species directories in `results/compleasm_proteome/` directory.


#### To run BUSCO:

1. Start the virtual environment for BUSCO.

2. Go to the `methods/` directory

3. Run the script `download_lineage_busco.sh` with SLURM to download `mammalia_odb10` lineage for BUSCO.
All files will be saved in `/storage/EasyVectorOmics/busco/methods/busco_downloads/lineages/mammalia_odb10/` directory.

4. Go to the `results/busco_proteome/mammalia_lineage/` directory.

5. Run the script `busco.sh` with SLURM to perform BUSCO.

This will generate files in corresponding species directories.
**NOTE**: Don't create the directories for the species because the script creates it automatically, otherwise there will be double directories, for example `canis/canis/`

Run the jobs separately for each species, change the corresponding paths in the script.


#### To extract 1:1 orthologs from BUSCO results

1. Rename the `full_table.tsv` files for each species in their directories in `mammalia_lineage/` so that they have following pattern: `species_busco_output.tsv`
Move the files to the `mammalia_lineage/` directory.

2. Deactivate the BUSCO and conda base environment.
`conda deactivate` (to deactivate compleasm venv)
`conda deactivate` (once again, to deactivate conda venv)


3. Activate Python virtual environment for python scripts.

4. Go to the `methods/` directory.

5. Uncomment lines needed for BUSCO in the script `extract_sco.py`. 
Run the script with SLURM to extract single-copy orthologs from the BUSCO output tables. 
This will generate in `results/busco_proteome/mammalia_lineage/` directory files with the pattern `sco_*_busco.tsv`.

6. Uncomment lines needed for BUSCO in the script `one_to_one.py`.
Run the script with SLURM to extract 1:1 orthologs of the analysed pairs from single-copy-othologs tables.
This will generate in `results/busco_proteome/mammalia_lineage/` directory files with the pattern `one_to_one_busco_*.tsv`.

7. Uncomment lines needed for BUSCO in the script `map_transc_prot_ids.py`.
Run the script with SLURM to map `transcript_ids` to protein and gene IDs.
This will generate in `results/busco_proteome/mammalia_lineage/` directory files with the pattern `mapped_1to1_busco_*.tsv`


#### To extract 1:1 orthologs from Compleasm results

1. Rename the `full_table.tsv` files for each species so that they follow the pattern
   `species_compleasm_output.tsv`.
   Move the files to the `longest_isoform/` directory.

2. Deactivate the Compleasm and conda base environments:

   ```bash
   conda deactivate   # deactivate compleasm venv
   conda deactivate   # deactivate conda base venv
   ```

3. Activate the Python virtual environment for Python scripts.

4. Go to the `methods/` directory.

5. Uncomment the lines required for Compleasm in the script `extract_sco.py`.
   Run the script with SLURM to extract single-copy orthologs from the Compleasm output tables.
   This generates files in `results/compleasm_proteome/longest_isoform/` with the pattern
   `sco_*_compleasm.tsv`.

6. Uncomment the lines required for Compleasm in the script `one_to_one.py`.
   Run the script with SLURM to extract 1:1 orthologs of the analyzed species pairs from the
   single-copy ortholog tables.
   This generates files in `results/compleasm_proteome/longest_isoform/` with the pattern
   `one_to_one_compleasm_*.tsv`.

7. Uncomment the lines required for Compleasm in the script `map_transc_prot_ids.py`.
   Run the script with SLURM to map transcript IDs to protein and gene IDs.
   This generates files in `results/compleasm_proteome/longest_isoform/` with the pattern
   `mapped_1to1_compleasm_*.tsv`.

### Results navigation

After the BUSCO and Compleasm preparation steps, the results have the following structure
under `/storage/EasyVectorOmics/busco/results/`:

```
results/
  |-- busco_proteome/
  |   |-- mammalia_lineage/
  |       |-- canis/
  |       |-- macaca/
  |       |-- rat/
  |       |-- mus/
  |       |-- sco_*_busco.tsv              <- single-copy orthologs per species
  |       |-- one_to_one_busco_*.tsv       <- extracted 1:1 ortholog pairs
  |       |-- mapped_1to1_busco_*.tsv      <- 1:1 pairs with mapped gene and protein IDs
  |
  |-- compleasm_proteome/
  |   |-- longest_isoform/
  |       |-- canis/
  |       |-- macaca/
  |       |-- rat/
  |       |-- mus/
  |       |-- sco_*_compleasm.tsv          <- single-copy orthologs per species
  |       |-- one_to_one_compleasm_*.tsv   <- extracted 1:1 ortholog pairs
  |       |-- mapped_1to1_compleasm_*.tsv  <- 1:1 pairs with mapped gene and protein IDs
  |
  |-- evaluation/
      |-- BUSComp_set.tsv                  <- union of BUSCO and Compleasm 1:1 orthologs
      |-- busco_compleasm_metrics.tsv      <- evaluation metrics for BUSCO vs Compleasm
      |-- all_false_negatives.tsv          <- BUSCO pairs missing in Compleasm
      |-- all_false_positives.tsv          <- Compleasm pairs missing in BUSCO
```



## PTP

### Material

Orthofinder2 results:

- `/storage/EasyVectorOmics/FastQ_GSE125483_JK/proteom/OrthoFinder/Results_Aug15/MultipleSequenceAlignments` - MSA results
- `/storage/EasyVectorOmics/FastQ_GSE125483_JK/proteom/OrthoFinder/Results_Aug15/Resolved_Gene_Trees` - Trees results

### Workflow (detailed instruction below)

1. `orthofinder_msa.sh` - Runs OrthoFinder2 with `-msa` flag
2. `fix_labels.sh` - changes lables after OrthoFinder to the PTP input format
3. `phylotreepruner.sh` - executes PTP
4. `extract_one_to_one.sh` (runs `extract_one_to_one.py`) - extracts 1:1 orthologs
5. `check_one_to_one_ptp.sh` (runs `check_one_to_one_ptp.py`) - check if 1:1 relationship is kept



### INSTRUCTIONS

Before running the pipeline, make sure to install all required software.
Prepare the `material/input_orthofinder/` directory with proteome FASTA files.

The scripts are located on the server in
`/storage/EasyVectorOmics/phylotreepruner/methods/`,
and the analysis was performed relative to that location.

The results are located in
`/storage/EasyVectorOmics/phylotreepruner/results/`.

---

### PTP installation

PTP was installed on the server as follows:

```bash
mamba create -n ptp-java -c conda-forge openjdk=11
mamba activate ptp-java
wget https://downloads.sourceforge.net/project/phylotreepruner/src_and_wrapper_scripts.zip
unzip src_and_wrapper_scripts.zip
```

---

### Virtual environments

**venv for PTP:**

```bash
eval "$(/storage/EasyVectorOmics/busco/miniconda3/bin/conda shell.bash hook)"
eval "$(mamba shell hook --shell bash)"
mamba activate ptp-java
```

**venv for Python scripts:**

```bash
source /storage/EasyVectorOmics/cactus/venv-analysis/bin/activate
```

---

### Workflow

1. **Run OrthoFinder2 with MSA option:**
   If needed, run the script `orthofinder_msa.sh`.
   This creates the directory `material/input_orthofinder/OrthoFinder/` with results inside.
   For PTP, the directories `MultipleSequenceAlignments/` and `Resolved_Gene_Trees/` are required.

2. **Change headers of FASTA files and trees:**
   OrthoFinder results contain sequence names in the form `species_ID`, but PTP requires the form `species|ID`.
   Run the script `fix_labels.sh`.
   This creates the directories:

   * `results/species_type/MSA_fixed/`
   * `results/species_type/Trees_fixed/`

3. **PhyloTreePruner:**
   Run the script `phylotreepruner.sh`.
   This creates FASTA files with the pattern `*.fa_pruned.fa` in
   `results/species_type/MSA_fixed/`,
   and log files in
   `results/species_type/ptp_logs/`.

4. **Extract 1:1 ortholog pairs:**
   Uncomment lines based on the input.
   Run the script `extract_one_to_one.sh`.
   This creates
   `results/species_type/ptp_one_to_one_pairs.tsv`.

5. **Check 1:1 relationships:**
   Uncomment lines based on the input.
   Run the script `check_one_to_one_ptp.sh`.
   Results are printed to stdout.

### Results navigation

After the PTP-based ortholog inference, the results have the following structure
under `/storage/EasyVectorOmics/phylotreepruner/results/`:

```
results/
  |-- mammalia/
      |-- MSA_fixed/
      |   |-- *.fa                         <- fixed MSAs after label conversion
      |   |-- *.fa_pruned.fa               <- pruned MSAs after PhyloTreePruner
      |
      |-- Trees_fixed/
      |   |-- *.tree                       <- gene trees with fixed labels
      |
      |-- ptp_logs/
      |   |-- *.log                        <- PhyloTreePruner log files
      |
      |-- ptp_one_to_one_pairs.tsv         <- extracted 1:1 ortholog pairs
```

**Notes:**

* Only ortholog pairs that pass PTP pruning and subsequent 1:1 filtering are included in
  `ptp_one_to_one_pairs.tsv`.
* These results are later used as input for reference-based evaluation against Cactus,
  BUSComp-set, and CactBUSComp-set.


## Preparing Orthofinder2 1:1 orthologs

### Material

- `/storage/EasyVectorOmics/FastQ_GSE125483_JK/proteom/OrthoFinder/Results_Aug15/Orthologues` OrthoFinder2 results

### Instructions

Run the script `extract_one_to_one_orthofinder.sh` (executes `extract_one_to_one_orthofinder.py`) to extract ready-to-use 1:1 ortholog pairs from OrthoFinder2. This will create `results/evaluation/mammalia/vs_orthofinder/orthofinder_1to1.tsv` file. 

## Evaluation of methods against reference sets 

This chapter contains information about evaluation of multiple 1:1 ortholog
identification methods using reference-based metrics across mammalian species.


### REFERENCE ORTHOLOG SETS

Following methods were used as reference:

  - Cactus-based approach for 1:1 ortholog identification
  - Union set of BUSCO and Compleasm 1:1 orthologs (BUSComp-set)
  - Union set of Cactus and BUSComp-set orthologs (CactBUSComp-set)


### EVALUATED METHODS

1:1 orthologs of following methods were evaluated against the references:

  - prot-syn - N-score and Sum mods, no intermediate transitivity
  - Tree-algorithm - majority, max_score (standard) and whitelist mods
  - PhyloTreePruner (PTP)
  - OrthoFinder2

Methods were also evaluated against each other (Jaccard Similarity only).

Evaluation metrics:
  Precision, Recall, F1-score, Jaccard Similarity

For the evaluation of prot-syn vs cactus, False Negatives (FN) and False Positives (FP)
were analysed in details.


### EVALUATION OUTPUTS

Two heatmaps were generated:

  1) Methods vs reference sets
     `results/evaluation/mammalia/metrics_methods_vs_reference_sets_v2.png`

  2) Jaccard similarity: methods vs each other and reference sets
     `results/evaluation/mammalia/all_vs_all/combined_jaccard_heatmap_v3.png`

Plots were copied to:
  `gene_relationship_classifier/plots/mammals/evaluation/`


### PROJECT NAVIGATION

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

### SCRIPTS AND EXECUTION ENVIRONMENT

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


### FILES USED

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

prot-syn (nscore approach):
  - `/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/one_to_one_nscore.tsv`

prot-syn (sum approach):
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

Dictionaries from prot-syn:
  - `/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/`

HAL graph (Cactus raw output):
  - `results/output_alignment_sexspecies.hal`

### EVALUATION STEPS

#### 1) BUSCO vs Compleasm + BUSComp-set creation

Task:
  To evaluate BUSCO (reference) against Compleasm (test) and create BUSComp-set

Script:
  
  `/storage/EasyVectorOmics/busco/methods/evaluate_busco_vs_compleasm.sh`
  (executes `evaluate_busco_vs_compleasm_v2.py`)

Input:
  BUSCO and Compleasm 1:1 ortholog files and map geneID<->proteinID.

Output location:
  `/storage/EasyVectorOmics/busco/results/mammalia/evaluation/`

Creates:
  - `busco_compleasm_metrics.tsv` (and `.txt`) - evaluation metrics
  - `all_false_negatives.tsv` - pairs that are found in reference, but not found in test
    (all species pairs)
  - `all_false_positives.tsv` - pairs that are found in test, but not found in reference
    (all species pairs)
  - per-species-pair false negatives and false positives
  - `BUSComp_set.tsv` (and `.pkl`) - union set of BUSCO and compleasm orthologs


#### 2) BUSComp-set vs Cactus + CactBUSComp-set creation


Task:
  To evaluate BUSComp-set (reference) against Cactus (test)

Script:
  
  `evaluate_cactus_vs_buscomp.sh`
  (executes `evaluate_cactus_vs_buscomp_v2.py`)

Input:
  Cactus 1:1 ortholog files, map geneID<->proteinID, BUSComp-set.

Output location:
  `results/evaluation/mammalia/vs_buscomp_set/`

Creates:
  - `cactus_BUSComp_metrics.tsv` (and `.txt`) - evaluation metrics
  - `all_false_negatives_cactus.tsv` - pairs that are found in reference, but not found
    in test (all species pairs)
  - `all_false_positives_cactus.tsv` - pairs that are found in test, but not found
    in reference (all species pairs)
  - per-species-pair false negatives and false positives
  - `CactBUSComp_set.tsv` (and `.pkl`) - union set of Cactus and BUSComp-set orthologs


#### 3) Methods vs Cactus 

Task:
  To evaluate methods (test) against cactus (reference)

Script:
  
  `evaluate_methods_vs_cactus.sh` 
  (executes `evaluate_methods_vs_cactus_v2.py`)

Input:
  cactus 1:1 ortholog files, map geneID<->proteinID, and 1:1 ortholog files of the tested
  methods.

Output location:
  `results/evaluation/mammalia/` (in corresponding methods directories)

Creates:
  - `cactus_<method>_metrics.tsv` (and `.txt`) - Statistical data and evaluation metrics:
    Precision, Recall, F1-score and Jaccard Similarity
  - `all_false_negatives_<method>.tsv` - pairs that are found in reference, but not found
    in test (all species pairs)
  - `all_false_positives_<method>.tsv` - pairs that are found in test, but not found
    in reference (all species pairs)

Notes:
  See the script for exact paths.


#### 4) Methods vs BUSComp-set

Task:
  To evaluate methods (test) against BUSComp-set (reference)

Script:
  
  `evaluate_methods_vs_buscomp.sh`
  (executes `evaluate_methods_vs_buscomp_v2.py`)

Input:
  BUSComp-set file, map geneID<->proteinID, and 1:1 ortholog files of the tested methods.

Output location:
  `results/evaluation/mammalia/vs_buscomp_set/` (in corresponding methods directories)

Creates:
  - `BUSComp_<method>_metrics.tsv` (and `.txt`) - Statistical data and evaluation metrics:
    Precision, Recall, F1-score and Jaccard Similarity
  - `all_false_negatives_<method>.tsv` - pairs that are found in reference, but not found
    in test (all species pairs)
  - `all_false_positives_<method>.tsv` - pairs that are found in test, but not found
    in reference (all species pairs)
  - per-species-pair false negatives and false positives

Notes:
  See the script for exact paths.


#### 5) Methods vs CactBUSComp-set

Task:
  To evaluate methods (test) against CactBUSComp-set (reference)

Script:
  
  `evaluate_methods_vs_cactbuscomp.sh`
  (executes `evaluate_methods_vs_cactbuscomp_v2.py`)

Input:
  CactBUSComp-set file, map geneID<->proteinID, and 1:1 ortholog files of the tested
  methods.

Output location:
  `results/evaluation/mammalia/all_vs_cactbuscomp/` (in corresponding methods directories)

Creates:
  - `CactBUSComp_<method>_metrics.tsv` (and `.txt`) - Statistical data and evaluation metrics:
    Precision, Recall, F1-score and Jaccard Similarity
  - `all_false_negatives_<method>.tsv` - pairs that are found in reference, but not found
    in test (all species pairs)
  - `all_false_positives_<method>.tsv` - pairs that are found in test, but not found
    in reference (all species pairs)
  - per-species-pair false negatives and false positives

Notes:
  See the script for exact paths.


#### 6) Methods vs each other + Jaccard similarity heatmap

Task:
  To evaluate methods against each other and plot Jaccard similarity

Script:
  
  `evaluate_methods_vs_each_other_and_plot.sh`
  (executes `evaluate_methods_vs_each_other_and_plot_v3.py`)

Input:
  Ortholog files of methods prot-syn, Tree, PTP and OrthoFinder, and already
  evaluated TSV files with metrics.

Output location:
  `results/evaluation/mammalia/all_vs_all/`

Creates:
  - `combined_jaccard_heatmap_v3.png` - heatmap of Jaccard similarities
  - `combined_jaccard_matrix_v3.txt` - Jaccard as text file


#### 7) Heatmap: Methods vs reference sets

Task:
  To visualise evaluation of methods against reference sets (Cactus, CactBUSComp and
  BUSComp) as heatmap

Script:
  
  `visualize_evaluation.sh`
  (executes `visualize_evaluation_metrics_v2.py`)

Input:
  TSV with metrics for each method-reference pair

Output location:
  `results/evaluation/mammalia/`

Creates:
  - `metrics_methods_vs_reference_sets_v2.png`


## DISAGREEMENT INVESTIGATION

Disagreements between Cactus and prot-syn ortholog sets (both nscore and sum
approaches) were investigated in details.
Here examples for nscore approach are shown, for sum approach see the directory
`sum_approach/` in `results/evaluation/mammalia/vs_prot_syn/`


#### A) INVESTIGATION OF FN

##### A1) neighborhood_score_cactus.sh (executes neighborhood_score_cactus_fn_v4.py)

Analyzes FN ortholog pairs (Pairs found by Cactus, and missed by prot-syn) by
computing neighborhood scores (N-scores, number of shared neighbors between genes) for the Cactus pair and (when available)
the prot-syn assigned alternative, classifies the FN into prot-syn_no_result and
contradictory (Cactus and prot-syn find different candidates to a gene).

The analysis is at both gene-level (two rows per FN pair, each row for the gene from
the FN pair and it's partner gene assigned by Cactus and prot-syn, if exists) and pair-level
(one row per FN with mean/max N-score for multiple prot-syn candidates for a FN pair and score
fractions).

For thesis the file `/storage/EasyVectorOmics/cactus/slurm_outputs/neighborhood_score_cactus_fn_v4_1307.out`  with results for nscore approach and `neighborhood_score_cactus_fn_v4_1318.out` for SUM approach was used.



Requires:
  - neighborhood dictionary (`neighborhoods.pkl`)
  - relationship dictionary (`geometric_mean.pkl`)
  - geneID<->proteinID mapping TSV
  - prot-syn 1:1 ortholog TSV
  - FN pairs TSV

Input dictionaries:
  `/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/`

prot-syn orthologs:
  `/storage/EasyVectorOmics/synteny_algorithm/results/mammalian/
  one_to_one_nscore.tsv and one_to_one_sum.tsv`

FN table:
  `results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/
  all_false_negatives_nscore.tsv`

Outputs:
  `results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/
  disagreement_investigation/false_negatives/`

Outputs include:
  - `false_negatives_with_ids_pair_level.tsv` - FN pairs with assigned FN_IDs for each pair
  - `false_negatives_classified_gene_level.tsv` - classified gene-level table (one gene per row with ID of its FN pair, 
    N-scores, fractions, relative-score columns)
  - `false_negatives_edge_affected_gene_rows.tsv` - genes affected by max n-score < 20
    ("edge-affected" genes)
  - `false_negatives_classified_pair_level.tsv` - classified pair-level table (one FN pair per row, N-scores, fractions)
  - `fn_distributions_2x2_gene_level.png` and `fn_distributions_2x2_pair_level.png` -
    two 2x2 distribution plots (gene-level and pair-level)
  - `contradictory_cactus_higher_nscore_than_protsyn_gene_level.tsv` and
    `anomalous_protsyn_nscore_zero_gene_level.tsv` - two gene-level subsets
    (contradictory cases, where cactus chose a candidate with better N-score and anomalous prot-syn score=0 cases)

Note: - `false_negatives_edge_affected_gene_rows.tsv` - genes classified as edge-affected
  with respect to the relative N-score, where the maximum achievable theoretical
  N-score is low (max N-score < 20). These cases reflect limitations imposed by
  neighborhood size rather than poor synteny support.
    

##### A2) check_in_excluded_fn.sh (executes check_in_excluded_fn.py)

Checks whether FN gene IDs appear in species-specific gene lists
(dog, rat, mouse) excluded by prot-syn.

Each row in the FN table represents one gene from an FN pair.

Requires:
  - FN table
  - excluded-gene tables located in `material/sex_experiment/`
    `excluded_genes_<species>.tsv` (Provided files for Canis, Rat and Mouse;
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

##### A3) max_chr_distance.sh (executes max_chr_distance.py)      

Analyzes chromosomal distances for **contradictory false negative pairs** by comparing
Cactus synteny region length with the maximum genomic distance to scored neighbors
from the prot-syn neighborhood.

For thesis the file `/storage/EasyVectorOmics/cactus/slurm_outputs/max_chr_distance_1326.out` with results for nscore approach and for SUM approach was used

For each contradictory FN pair, the script computes:

* **Cactus region length** (geometric mean of query and target region lengths)
* **Maximum neighborhood distance**: distance between the focus gene and its furthest
  scored neighbor on the same chromosome

Only pairs where **both values are available** are retained.

Requires:

* gene-level FN table
* neighborhood dictionary (`neighborhoods.pkl`)
* prot-syn score dictionary (`geometric_mean.pkl`)
* species-specific GTF files
* Cactus gene→region tables
* Cactus region tables

Outputs:

* TSV file with per-pair distances
* Scatter plot comparing Cactus region length vs neighborhood max distance



##### A4) split_contradictory_fn.sh (split_contradictory_fn.py)

For thesis the file `/storage/EasyVectorOmics/cactus/slurm_outputs/max_chr_distance_1326.out` with results for nscore approach and for SUM approach was used

Splits **contradictory false negatives** into two complementary representations:

* **Cactus-focused pairs**
  Focus_Gene → Cactus_Ortholog_Gene
* **prot-syn-focused pairs**
  Focus_Gene → ProtSyn_Ortholog_Gene

The script produces:

* Pair tables suitable for coordinate extraction and Cactus scoring
* A focus–candidate mapping table used for downstream score comparison plots

Outputs are written with **Unix newlines** to avoid downstream awk/bedtools issues.

Outputs:

* `*_c_pairs.tsv`
* `*_n_pairs.tsv`
* `*_focus_candidate_map.tsv`


##### A5) Calculate Cactus Support Score for FN

See the detailed script information in FP section. 

The paths were changed according to the prot-syn approach, and whether it's a focus-cactus or focus-prot-syn pair.

Contradictory FN pairs are further analyzed by computing **Cactus Support Scores**
for both candidate relationships:

* Focus → Cactus ortholog
* Focus → prot-syn ortholog

Workflow:

1. `split_contradictory_fn.sh` produces pair tables
2. `build_pair_coordinate_table.sh` extracts genomic coordinates
3. `get_cactus_score.sh` computes alignment-based support scores

Paths are adjusted automatically depending on:

* prot-syn approach (nscore / sum)
* focus type (cactus / prot-syn)

Note:
In the output tables used for Cactus Support Score calculation, both cactus-focused
and prot-syn-focused split tables retain column names with the `_protsyn` suffix.
This naming reflects the original table schema and does not indicate the origin
of the candidate; cactus-focused pairs are identified by the focus type and input
path, not by column names.



##### A6) plot_cactus_score_contradictory_fn.sh (executes plot_cactus_score_contradictory_fn.py)

Creates a **scatter plot** comparing Cactus Support Scores for contradictory FN pairs:

* X-axis: CactusSupportScore for **Cactus pair**
* Y-axis: CactusSupportScore for **prot-syn pair**

Each point represents **one contradictory focus gene**.

Outputs:

* merged TSV with both scores
* scatter plot PNG

Missing scores are excluded by default.


### B) INVESTIGATION OF FP

#### B1) find_c_no_results_fp.sh (executes find_c_no_results_fp.py)

Among FP pairs finds those that have no result in Cactus, and saves them for further
investigation.

For thesis the file `/storage/EasyVectorOmics/cactus/slurm_outputs/find_c_no_results_fp1302.out with results for nscore approach and `find_c_no_results_fp1346.out` for SUM approach was used

Logic:
  - Check Gene1: if found anywhere in cactus -> ignore the whole FP pair
  - If Gene1 not found -> check Gene2: if found anywhere -> ignore
  - Only if BOTH genes are not found anywhere -> C_no_result

Input:
  - `all_false_positives_nscore.tsv` - False positives table:
    `results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/`

  - `can_mac_one2one_fix.tsv` - Cactus ortholog file (Dog-Macaca example):
    `results/reciprocal_pairs/`

Output:
  `all_C_no_results.tsv` - FP pairs where both genes have no Cactus result:
  
  - `results/evaluation/mammalia/vs_prot_syn/nscore_approach/test_no_transitivity/
  disagreement_investigation/false_positives/`


#### B2) build_pair_coordinate_table.sh (located in methods/)

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
  - `gtf_Canis.gtf`
  - `gtf_Macaca.gtf`
  - `gtf_Mus.gtf`
  - `gtf_Rat.gtf`

Output TSV columns:

  `gene1 gene2 species1 species2 chr1 chr2 start1 start2 end1 end2 length1 length2`

Coordinate conversion:
  - GTF: 1-based inclusive [start1..end1]
  - HAL: start0 = start1-1, length = end1-start1+1, end0 = start0+length


#### B3) get_cactus_score.sh (located in methods/)

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


#### B4) plot_cactus_score_no_result_fp.sh (executes plot_cactus_score_no_result_fp.py)   


Visualizes the distribution of **Cactus Support Scores** for false positive pairs
classified as **C_no_result** (genes not found anywhere in Cactus).

The script:

* loads the score TSV produced by `get_cactus_score.sh`
* plots a histogram of the Score column

Output:

* histogram PNG showing score distribution for C_no_result FP pairs

### NEIGHBORHOOD LENGTHS FOR FN AND FP

For thesis the file `/storage/EasyVectorOmics/cactus/slurm_outputs/neighborhood_lengths1329.out` with results for nscore approach and `find_c_no_results_fp1347.out` for SUM approach was used

Script:
`neighborhood_lengths.sh` (executes `neighborhood_lengths.py`)

Analyzes **genomic neighborhood lengths** for three categories:

1. All genes (baseline)
2. False Negatives
3. False Positives

For each gene:

* neighborhood span is calculated as
  `max(end) − min(start)` across all neighbors

For FN and FP pairs:

* neighborhood lengths are computed for both genes
* the **mean neighborhood length** is used as pair-level value

Requires:

* neighborhood dictionary (`neighborhoods.pkl`)
* gene-level FN table
* FP table
* species-specific GTF files

Outputs:

* Raincloud (or violin) plots per species pair
* Summary statistics TSV
* Detailed per-gene and per-pair TSV

## REFERENCES

- Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov.  
  *BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper
  Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes*.  
  **Molecular Biology and Evolution**, Volume 38, Issue 10, October 2021, Pages 4647–4654.  
  [https://doi.org/10.1093/molbev/msab199](https://doi.org/10.1093/molbev/msab199)

- Neng Huang, Heng Li.  
  *compleasm: a faster and more accurate reimplementation of BUSCO*.  
  **Bioinformatics**, 39, btad595, 2023.  
  [https://doi.org/10.1093/bioinformatics/btad595](https://doi.org/10.1093/bioinformatics/btad595)

- Glenn Hickey, Benedict Paten, Dent Earl, Daniel Zerbino, David Haussler.  
  *HAL: A Hierarchical Format for Storing and Analyzing Multiple Genome Alignments*.  
  **Bioinformatics**, 2013.  
  [https://doi.org/10.1093/bioinformatics/btt128](https://doi.org/10.1093/bioinformatics/btt128)

- Armstrong, J., Hickey, G., Diekhans, M. et al.  
  *Progressive Cactus is a multiple-genome aligner for the thousand-genome era*.  
  **Nature**, 587, 246–251 (2020).  
  [https://doi.org/10.1038/s41586-020-2871-y](https://doi.org/10.1038/s41586-020-2871-y)

- Kocot KM, Citarella MR, Moroz LL, Halanych KM.  
  *PhyloTreePruner: A Phylogenetic Tree-Based Approach for Selection of Orthologous Sequences for Phylogenomics*.  
  **Evolutionary Bioinformatics Online**, 2013;9:429–435.  
  [https://doi.org/10.4137/EBO.S12813](https://doi.org/10.4137/EBO.S12813)

- Krasheninnikova K, Diekhans M, Armstrong J, Dievskii A, Paten B, O’Brien S.  
  *halSynteny: a fast, easy-to-use conserved synteny block construction method for multiple whole-genome alignments*.  
  **GigaScience**, 2020;9(6).  
  [https://doi.org/10.1093/gigascience/giaa047](https://doi.org/10.1093/gigascience/giaa047)

- Buchfink B, Reuter K, Drost HG.  
  *Sensitive protein alignments at tree-of-life scale using DIAMOND*.  
  **Nature Methods**, 2021;18(4):366–368.  
  [https://doi.org/10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x)

- Quinlan AR, Hall IM.  
  *BEDTools: a flexible suite of utilities for comparing genomic features*.  
  **Bioinformatics**, 2010;26(6):841–842.  
  [https://doi.org/10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)

- Neph S, Kuehn MS, Reynolds AP, et al.  
  *BEDOPS: high-performance genomic feature operations*.  
  **Bioinformatics**, 2012;28(14):1919–1920.  
  [https://doi.org/10.1093/bioinformatics/bts277](https://doi.org/10.1093/bioinformatics/bts277)

- Pertea G, Pertea M.  
  *GFF Utilities: GffRead and GffCompare*.  
  **F1000Research**, 2020;9:304.  
  [https://doi.org/10.12688/f1000research.23297.2](https://doi.org/10.12688/f1000research.23297.2)

- Dainat J, Cannoodt R, Soares A, et al.  
  *NBISweden/AGAT: AGAT V1.6.1*.  
  **Zenodo**, CERN, published January 13, 2026.  
  [https://doi.org/10.5281/zenodo.3552717](https://doi.org/10.5281/zenodo.3552717)

- David M Emms, Yi Liu, Laurence Belcher, Jonathan Holmes, Steven Kelly.  
  *OrthoFinder: scalable phylogenetic orthology inference for comparative genomics*.  
  **bioRxiv**, 2025.07.15.664860.  
  [https://doi.org/10.1101/2025.07.15.664860](https://doi.org/10.1101/2025.07.15.664860)


