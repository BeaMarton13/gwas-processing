#!/bin/bash

# Usage check
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 file1.bed file2.bed new_bed_name_w_path"
    exit 1
fi

BED1=$1
BED2=$2
BED_NAME=$3

# Check if bedtools is installed
if ! command -v bedtools &> /dev/null
then
    echo "bedtools could not be found. Please install it first."
    exit 1
fi

# Create output directory
mkdir -p bed_comparison

# 1️⃣ Intersect (common regions)
bedtools intersect -a "$BED1" -b "$BED2" -wa -wb \
  | grep -v "transcript" \
  | uniq > ${BED_NAME}
  
# # 2️⃣ Unique to BED1 (include BED2 filename)
# bedtools subtract -a "$BED1" -b "$BED2" > bed_comparison/unique_to_${BED1_NAME}_vs_${BED2_NAME}.bed
# echo "Regions unique to $BED1 written to bed_comparison/unique_to_${BED1_NAME}_vs_${BED2_NAME}.bed"

# # 3️⃣ Unique to BED2
# bedtools subtract -a "$BED2" -b "$BED1" > bed_comparison/unique_to_${BED2_NAME}_vs_${BED1_NAME}.bed
# echo "Regions unique to $BED2 written to bed_comparison/unique_to_${BED2_NAME}_vs_${BED1_NAME}.bed"