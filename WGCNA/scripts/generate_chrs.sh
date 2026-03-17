#!/bin/bash

mkdir -p ./magma/ref/plink_chr

for chr in {1..22}; do
  vcf="./ld_vcf/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"

  if [[ ! -f "$vcf" ]]; then
    echo "WARNING: $vcf not found, skipping"
    continue
  fi

  echo "Processing $vcf"
  plink2 \
    --vcf "$vcf" \
    --max-alleles 2 \
    --snps-only just-acgt \
    --make-bed \
    --out ./magma/ref/plink_chr/chr${chr}
done