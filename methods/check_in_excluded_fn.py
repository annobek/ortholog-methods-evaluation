import argparse
import os
import pandas as pd

"""
This script checks whether False Negative (FN) gene IDs appear in
species-specific excluded gene lists (dog, rat, mouse).

Each row in the FN table represents one gene from an FN pair.
The script:
  - determines the species of each FN gene
  - searches the corresponding excluded-gene table
  - flags FN rows whose gene IDs were excluded

Inputs:
  1. FN table (TSV)
     Required columns:
       - Focus_Gene (gene ID, string or int)
       - Focus_Species (species name, string)

  2. Excluded gene tables (TSV), one per species (dog, rat, mouse)
     Required columns:
       - gene_id (gene ID, string or int)

Outputs:
  Two TSV files written to the output directory:
    1. fn_rows_found_in_excluded.tsv
       - Subset of FN rows where Focus_Gene is found in the excluded list
       - Includes all original FN columns + helper columns

    2. summary.tsv
       - Per-species summary:
           species_key | total_fn_rows | excluded_hits

Expected file format:
    - All input files must be TAB-separated (.tsv)
    - gene IDs are compared as strings (whitespace stripped)

Usage example:
    python check_fn_vs_excluded.py \
      --fn FN_table.tsv \
      --excluded-dog excluded_dog.tsv \
      --excluded-rat excluded_rat.tsv \
      --excluded-mouse excluded_mouse.tsv \
      --out-dir results/false_negatives/

NOTE: For now hardcoded for 3 species: Dog, rat and mouse. If more species required, need to change some code parts.      
"""


def normalize_species_name(s: str) -> str:
    """
    Map species names from the FN table to internal species keys.

    Parameters
    ----------
    s : str
        Species name from the FN table (Focus_Species column).

    Returns
    -------
    str
        One of:
          - "dog"
          - "rat"
          - "mouse"
          - "unknown" (if species cannot be mapped)
    """
    s = (s or "").lower()

    if "canis" in s:
        return "dog"
    if "rattus" in s:
        return "rat"
    if "mus" in s:
        return "mouse"

    return "unknown"


def load_excluded_table(path: str) -> set:
    """
    Load an excluded gene table and return its gene IDs as a set.

    Parameters
    ----------
    path : str
        Path to a TSV file containing an excluded gene table.
        The file must have a column named 'gene_id'.

    Returns
    -------
    set[str]
        Set of gene IDs (as strings) for fast membership testing.
    """

    df = pd.read_csv(path, sep="\t", dtype={"gene_id": "string"})
    return set(df["gene_id"].astype("string").str.strip())


def main(fn_path, dog_path, rat_path, mouse_path, out_dir):

    # Create output directory
    os.makedirs(out_dir, exist_ok=True)

    # Load FN table
    fn = pd.read_csv(fn_path, sep="\t", dtype={"Focus_Gene": "string"})
    fn["Focus_Gene"] = fn["Focus_Gene"].str.strip()
    fn["species_key"] = fn["Focus_Species"].apply(normalize_species_name)

    # Load excluded gene sets
    excluded = {
        "dog": load_excluded_table(dog_path),
        "rat": load_excluded_table(rat_path),
        "mouse": load_excluded_table(mouse_path),
    }

    # Check exclusion
    fn["is_excluded"] = fn.apply(
        lambda r: r["Focus_Gene"] in excluded.get(r["species_key"], set()),
        axis=1,
    )

    # Save excluded FN rows
    fn_excluded = fn[fn["is_excluded"]]
    fn_excluded.to_csv(
        os.path.join(out_dir, "fn_genes_found_in_excluded.tsv"),
        sep="\t",
        index=False,
    )

    # Save summary
    summary = (
        fn.groupby("species_key")["is_excluded"]
        .agg(total_fn_rows="size", excluded_hits="sum")
        .reset_index()
    )
    summary.to_csv(
        os.path.join(out_dir, "fn_genes_found_in_excluded_summary.tsv"),
        sep="\t",
        index=False,
    )

    print("Done.")
    print(f"Results written to: {out_dir}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fn", required=True, help="FN table TSV")
    ap.add_argument("--excluded-dog", required=True, help="Dog excluded genes TSV")
    ap.add_argument("--excluded-rat", required=True, help="Rat excluded genes TSV")
    ap.add_argument("--excluded-mouse", required=True, help="Mouse excluded genes TSV")
    ap.add_argument(
        "--out-dir",
        required=True,
        help="Output directory for result files",
    )
    args = ap.parse_args()

    main(
        args.fn,
        args.excluded_dog,
        args.excluded_rat,
        args.excluded_mouse,
        args.out_dir,
    )
