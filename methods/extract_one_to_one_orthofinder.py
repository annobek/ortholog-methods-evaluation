import pandas as pd
from pathlib import Path
from itertools import combinations

# Mammals
#OF_DIR = Path("/storage/EasyVectorOmics/FastQ_GSE125483_JK/proteom/OrthoFinder/Results_Aug15/Orthologues")
# Plants
#OF_DIR = Path("/storage/EasyVectorOmics/cardamine/material/OrthoFinder/Results_Oct31/Orthologues")
# Drosophila
OF_DIR = Path("/media/BioNAS/ag_hallab/EasyVectorOmics/analyze_sequences/material/dsim_dmel_dsec/OrthoFinder/Results_Dec02/Orthologues")

# for mammals
'''
species_names = [
    "Canis_lupus",
    "Macaca_fascicularis",
    "Rattus_norvegicus",
    "Mus_musculus",
    "arabidopsis",
    "cardamine"
]

species_filenames = [
    "Canis_lupus_protein",
    "Macaca_fascicularis_protein",
    "Rattus_norvegicus_protein",
    "Mus_musculus_protein",
    "arabidopsis",
    "cardamine"
]
'''

# for plants
'''
species_names = [
    "arabidopsis",
    "cardamine"
]

species_filenames = [
    "arabidopsis",
    "cardamine"
]
'''

# for drosophila
species_names = [
    "Drosophila melanogaster",
    "Drosophila sechellia",
    "Drosophila simulans"
]

species_filenames = [
    "larger_seq_dmel",
    "larger_seq_dsec",
    "larger_seq_dsim"
]




def extract_1to1(species1, species2, output=None):
    """
    Extracts 1:1 orthologs between two species from OrthoFinder's nested Orthologues directories.
    Args:
        species1 (str): query species (e.g. "Canis_lupus_protein")
        species2 (str): target species (e.g. "Macaca_fascicularis_protein")
        output (str, optional): if given, saves the result as TSV
    Returns:
        pd.DataFrame: with Species1, Protein1, Species2, Protein2 columns
    """

    subdir = OF_DIR / f"Orthologues_{species1}"
    file_path = subdir / f"{species1}__v__{species2}.tsv"

    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    df = pd.read_csv(file_path, sep="\t")

    # Get actual column names
    col1, col2 = df.columns[1], df.columns[2]

    # Split genes by comma+space
    df["genes1"] = df[col1].str.split(", ")
    df["genes2"] = df[col2].str.split(", ")

    # Keep only 1:1 orthologs
    df_1to1 = df[df["genes1"].str.len().eq(1) & df["genes2"].str.len().eq(1)]

    result = pd.DataFrame({
        "Species1": species1.replace("_protein", ""),
        "Protein1": df_1to1["genes1"].str[0],
        "Species2": species2.replace("_protein", ""),
        "Protein2": df_1to1["genes2"].str[0],
    })

    if output:
        result.to_csv(output, sep="\t", index=False)
        print(f"Saved {len(result)} 1:1 orthologs to {output}")

    return result




all_results = []
for s1, s2 in combinations(species_filenames, 2):
    df = extract_1to1(s1, s2)
    all_results.append(df)

all_1to1 = pd.concat(all_results, ignore_index=True)
# mammals
#all_1to1.to_csv("results/evaluation/vs_orthofinder/orthofinder_1to1.tsv", sep="\t", index=False)
# plants
#all_1to1.to_csv("results/evaluation/plants/vs_orthofinder/orthofinder_1to1.tsv", sep="\t", index=False)
# drosophila
all_1to1.to_csv("results/evaluation/drosophila/vs_orthofinder/orthofinder_1to1.tsv", sep="\t", index=False)

