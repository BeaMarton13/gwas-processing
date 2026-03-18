#!/bin/bash

# Check arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 tsv_file_w_path (e.g. gene_networks/dementia001.tsv)"
    exit 1
fi

tsv_file_w_path="$1"
python3 src/utils/export_network_from_tsv.py "$tsv_file_w_path"

graph_file="${tsv_file_w_path%.tsv}.graphml"

if [ -f "$graph_file" ]; then
    echo "Graph file found: $graph_file"
else
    echo "Error: Graph file not found: $graph_file"
    exit 1
fi

