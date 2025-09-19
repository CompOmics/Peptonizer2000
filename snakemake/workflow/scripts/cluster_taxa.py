import argparse

from peptonizer_rust import cluster_taxa_py

parser = argparse.ArgumentParser(description = 'cluster Taxa based on peptidome similarity and weight attributed')

parser.add_argument(
    '--full-graphml-path',
    type = str,
    help = 'Path(s) to the full Peptonizer graphml file for which you wish to cluster taxa (not containing communities).'
)
parser.add_argument(
    '--taxa-weights-dataframe-file',
    type = str,
    help = 'Path to file with weighted taxa computed in the taxa weighing step.'
)
parser.add_argument(
    '--similarity-threshold',
    type = float,
    help = 'Threshold for the peptidome similarity at which two taxa should belong to the same cluster.'
)
parser.add_argument(
    '--out',
    type = str,
    help= 'Path to clustered taxa output csv file.'
)

args = parser.parse_args()

with open(args.full_graphml_path, 'r') as graph_file, open(args.taxa_weights_dataframe_file, 'r') as taxa_weights_file:
    graph_xml = graph_file.read()
    taxa_weights_csv = taxa_weights_file.read()

clustered_taxa_csv = cluster_taxa_py(
    graph_xml,
    taxa_weights_csv,
    args.similarity_threshold
)

with open(args.out, 'w') as clustered_taxa_file:
    clustered_taxa_file.write(clustered_taxa_csv)
