#!/usr/bin/env python3
"""
Score distribution visualization script.

This script:
1. Loads a tab-separated table containing genomic alignment information
2. Extracts the 'Score' column
3. Plots a histogram showing the distribution of scores
4. Saves the plot to an output file

Arguments:
    --input   Path to input TSV file
    --output  Path to output image file (e.g. score_distribution.png)

Output:
    A histogram plot saved to disk

Returns:
    None (side effects: file I/O and plotting)
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt


def parse_arguments():
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Object containing input and output file paths.
    """
    parser = argparse.ArgumentParser(description="Plot distribution of Score values.")
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input TSV file containing a 'Score' column"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output image file (e.g. score_histogram.png)"
    )
    return parser.parse_args()


def load_table(path):
    """
    Load the input table as a pandas DataFrame.

    Parameters
    ----------
    path : str
        Path to TSV file.

    Returns
    -------
    pandas.DataFrame
        Loaded table.
    """
    return pd.read_csv(path, sep="\t")


def extract_scores(df):
    """
    Extract the Score column from the DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        Input table.

    Returns
    -------
    pandas.Series
        Series containing score values.
    """
    return df["Score"]


def plot_score_distribution(scores, output_path):
    """
    Plot and save a histogram of score values.

    Parameters
    ----------
    scores : pandas.Series
        Numeric score values.
    output_path : str
        File path where the plot will be saved.

    Returns
    -------
    None
    """
    plt.figure(figsize=(8, 6))
    plt.hist(scores, bins=30)
    plt.xlabel("Cactus Support Score")
    plt.ylabel("Frequency")
    plt.title("Distribution of Cactus Support Scores, C-no-results")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def main():
    """
    Main workflow function.

    Orchestrates argument parsing, data loading,
    score extraction, and plotting.
    """
    args = parse_arguments()
    df = load_table(args.input)
    scores = extract_scores(df)
    plot_score_distribution(scores, args.output)


if __name__ == "__main__":
    main()
