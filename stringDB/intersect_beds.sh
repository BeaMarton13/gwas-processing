#!/usr/bin/env bash
# ---------------------------------------------------------
# intersect_beds.sh
# Intersect multiple BED files using bedtools.
#
# Usage:
#   ./intersect_beds.sh output.bed file1.bed file2.bed [file3.bed ...]
#
# Requirements:
#   - bedtools must be installed and in your PATH
# ---------------------------------------------------------

# Exit on errors
set -euo pipefail

# Check for bedtools
if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools not found in PATH."
    echo "Install with: conda install -c bioconda bedtools"
    exit 1
fi

# Check arguments
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 output.bed file1.bed file2.bed [file3.bed ...]"
    exit 1
fi

# Parse arguments
output="$1"
shift
bed_files=("$@")

# Start with the first file
temp_file="${bed_files[0]}"

# Loop through the rest and intersect iteratively
for (( i=1; i<${#bed_files[@]}; i++ )); do
    echo "Intersecting: ${bed_files[i]} ..."
    bedtools intersect -a "$temp_file" -b "${bed_files[i]}"  -wa -wb  | grep -v "transcript"  | uniq > "${output}.tmp"
    mv "${output}.tmp" "$output"
    temp_file="$output"
done

echo "✅ Intersection complete. Result saved to: $output"
