import argparse

from peptonizer_rust import generate_pepgm_graph_py


parser = argparse.ArgumentParser(
    description="Run the PepGM algorithm from command line"
)

parser.add_argument(
    "--sequence-scores-dataframe-file",
    type=str,
    required=True,
    help="Dataframe file containing the taxa weights that have been computed before.",
)
parser.add_argument(
    "--out",
    type=str,
    required=True,
    help="Path to output file where GraphML will be saved.",
)

args = parser.parse_args()

with open(args.sequence_scores_dataframe_file, "r") as f:
    csv_str = f.read()

ct_factor_graph = generate_pepgm_graph_py(csv_str)

with open(args.out, 'w') as graphml_file:
    graphml_file.write(ct_factor_graph)
