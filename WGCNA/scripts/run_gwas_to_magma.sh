for f in project/gwas_raw/*.tsv.gz; do
  base=$(basename "$f" .tsv.gz)
  python project/scripts/gwas_to_magma_pval.py \
    "$f" "project/gwas_clean/${base}.magma.pval.tsv.gz"
done