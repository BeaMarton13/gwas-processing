#!/bin/bash

# Start from chr1 as the base
cp project/magma/ref/plink_chr/chr1.bed project/magma/ref/REF_GRCh38.bed
cp project/magma/ref/plink_chr/chr1.bim project/magma/ref/REF_GRCh38.bim
cp project/magma/ref/plink_chr/chr1.fam project/magma/ref/REF_GRCh38.fam

for chr in {2..22}; do
  echo "Merging chr${chr}..."
  plink2 \
    --bfile project/magma/ref/REF_GRCh38 \
    --bmerge project/magma/ref/plink_chr/chr${chr} \
    --make-bed \
    --out project/magma/ref/REF_GRCh38_tmp

  # Replace REF_GRCh38 with the merged output for the next iteration
  mv project/magma/ref/REF_GRCh38_tmp.bed project/magma/ref/REF_GRCh38.bed
  mv project/magma/ref/REF_GRCh38_tmp.bim project/magma/ref/REF_GRCh38.bim
  mv project/magma/ref/REF_GRCh38_tmp.fam project/magma/ref/REF_GRCh38.fam
done