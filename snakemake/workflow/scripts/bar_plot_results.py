import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import json
import argparse

from peptonizer_rust import get_names_for_taxa_py

"""
Script that takes PepGM .csv output, translates taxIDS to scientific names, and barplots the *number of results* highest
scoring taxa.
"""

def plot_peptonizer_results(input_file: str, output_file: str, number_of_taxa: int = 25):
    """
    Read the results of a Peptonizer run from a CSV-file (denoted by the input_file argument) and write bar charts
    representing these results to a PNG-file.
    """
    assert input_file.lower().endswith(".csv"), "Input file should be a CSV."
    assert output_file.lower().endswith(".png"), "Output file should be a PNG."

    # read csv using pandas
    ids = pd.read_csv(input_file, names=["ID", "score", "type"])
    tax_ids = ids.loc[ids["type"] == "taxon"]
    tax_ids = tax_ids.dropna()

    tax_ids.loc[:, "score"] = pd.to_numeric(tax_ids["score"], downcast="float")
    tax_ids = tax_ids.sort_values("score")
    taxa_check = tax_ids.ID.tolist()

    taxa_name_str = get_names_for_taxa_py([int(x) for x in tax_ids["ID"]])
    taxa_name_dict = json.loads(taxa_name_str)
    taxa_names = [taxa_name_dict[tax] for tax in taxa_check]
    scores = tax_ids["score"]

    # make the barplot
    fig, ax = plt.subplots()
    fig.set_size_inches(30, 15)
    bars = ax.barh(
        range(len(taxa_names[-number_of_taxa:])),
        scores[-number_of_taxa:],
        color="#283593",
    )

    ax.set_yticks(range(len(taxa_names[-number_of_taxa:])))
    ax.set_yticklabels(taxa_names[-number_of_taxa:], fontsize=24, color="#283593", fontweight="bold")
    ax.tick_params(axis='y', which='major', pad=15)
    plt.xlim((0, 1))
    plt.xlabel("Probability score", fontsize=35, fontweight="bold")
    ax.xaxis.set_ticks(np.arange(0, 1.2, 0.2))
    ax.xaxis.set_ticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.bar_label(bars, fmt='{:,.3f}', fontsize=24, fontweight='bold', color='black', padding=20)

    fig.tight_layout()

    plt.savefig(output_file)
    plt.close()

matplotlib.use("Agg")

parser = argparse.ArgumentParser(description="Generate BarPlot of PepGM results.")

parser.add_argument(
    "--results-file", type=str, help="Path(s) to your PepGM results CSV."
)
parser.add_argument(
    "--number-of-results",
    type=int,
    default=12,
    help="How many taxa you want to show up on the results plot.",
)
parser.add_argument(
    "--out",
    type=str,
    help="Path(s) to where the generated BarPlot PNG's should be stored.",
)

args = parser.parse_args()
plot_peptonizer_results(args.results_file, args.results_file.replace(".csv", ".png"), args.number_of_results)
