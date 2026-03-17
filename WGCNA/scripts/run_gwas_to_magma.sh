for f in ./gwas_raw/*.tsv.gz; do
  base=$(basename "$f" .tsv.gz)
  python ./scripts/gwas_to_magma_pval.py \
    "$f" "./gwas_clean/${base}.magma.pval.tsv.gz"
done