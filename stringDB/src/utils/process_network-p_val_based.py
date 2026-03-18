# Processes network:
    # - reads GML file
    # - clusters using InfoMap
    # - reads p-values from all GWAS files in the given directory
    # - selects "relevant clusters" based on avg p-values

import argparse

def handle_args():
    parser = argparse.ArgumentParser(
        description='Process network by detecting relevant clusters based on p-values.'
    )

    parser.add_argument(
        'network_file',
        type=str,
        help='Path to the network GML file, eg: gene_networks/dementia001.gml'
    )

    parser.add_argument(
        'network_file',
        type=str,
        help='Path to the network GML file, eg: gene_networks/dementia001.gml'
    )