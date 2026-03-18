#!/bin/bash

# set -euo pipefail

# Check for bedtools
if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools not found in PATH."
    echo "Install with: conda install -c bioconda bedtools"
    exit 1
fi

# Check arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 mental_disorder_name (e.g. dementia) processed_folder_name (e.g. processed or processed_0_01) p_value_limit (e.g. 0.01 or 1e-5) whole_word_match (True/False)"
    exit 1
fi

disease="$1"
processed_folder="$2"
p_value_limit="$3"
whole_world_match="$4"

( python3 main.py "$disease" "$processed_folder" "$p_value_limit" "$whole_word_match" )
wait
( python3 src/utils/convert_json_to_bed.py "$disease" "$processed_folder" )
wait

# for file in `ls data/"$disease"/"$processed_folder"/*.bed`; do
#     echo $file
#     base=$(basename "$file")  
#     out="bed_comparison/$disease/$processed_folder/new_${base}"
#     mkdir -p "bed_comparison/$disease/$processed_folder"
#     bedtools intersect -a data/homo_sapiens/homo_sapiens.bed -b "$file" -wa -wb > "$out"
# done

# Put all files into a Bash array
files=(bed_comparison/"$disease"/"$processed_folder"/new*.bed)

# Make sure there is at least one file
if [ ${#files[@]} -eq 0 ]; then
    echo "No files found!"
    exit 1
fi


# Start with the first file
bed_files=(bed_comparison/"$disease"/"$processed_folder"/new*.bed)

str="${p_value_limit//./_}"   # replace . with _

# NOTE run separately 
    # python3 src/utils/intersect_files.py dementia_0_01 data/homo_sapiens/homo_sapiens.bed data/dementia/processed_new_0_01/*.bed
    # because this:
# ( python3 src/utils/intersect_files.py "${disease}_${str}" data/homo_sapiens/homo_sapiens.bed "${bed_files[@]}" )
# wait
    # gives ./process_disease.sh: line 51:  9375 Killed: 9               ( python3 src/utils/intersect_files.py "${disease}_${str}" data/homo_sapiens/homo_sapiens.bed "${bed_files[@]}" )