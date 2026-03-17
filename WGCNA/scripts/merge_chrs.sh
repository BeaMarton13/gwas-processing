#!/bin/bash

# Start from chr1 as the base
cp ./magma/ref/plink_chr/chr1.bed ./magma/ref/REF_GRCh38.bed
cp ./magma/ref/plink_chr/chr1.bim ./magma/ref/REF_GRCh38.bim
cp ./magma/ref/plink_chr/chr1.fam ./magma/ref/REF_GRCh38.fam

for chr in {2..22}; do
  echo "Merging chr${chr}..."
  plink2 \
    --bfile ./magma/ref/REF_GRCh38 \
    --bmerge ./magma/ref/plink_chr/chr${chr} \
    --make-bed \
    --out ./magma/ref/REF_GRCh38_tmp

  # Replace REF_GRCh38 with the merged output for the next iteration
  mv ./magma/ref/REF_GRCh38_tmp.bed ./magma/ref/REF_GRCh38.bed
  mv ./magma/ref/REF_GRCh38_tmp.bim ./magma/ref/REF_GRCh38.bim
  mv ./magma/ref/REF_GRCh38_tmp.fam ./magma/ref/REF_GRCh38.fam
done