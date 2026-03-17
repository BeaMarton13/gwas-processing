#!/usr/bin/env python3
import gzip
import re
import sys

def parse_attrs(attr_str: str) -> dict:
    # GTF attributes look like: key "value"; key2 "value2";
    attrs = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_str):
        attrs[m.group(1)] = m.group(2)
    return attrs

def main(gtf_gz: str, out_loc: str):
    # Output columns: GENEID  CHR  START  END  STRAND  (tab-separated)
    # We'll use Ensembl gene_id (optionally with version removed) so it’s stable.
    seen = set()
    with gzip.open(gtf_gz, "rt") as f, open(out_loc, "w") as out:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs_s = parts
            if feature != "gene":
                continue

            # MAGMA expects autosomes as 1..22, and X/Y allowed; drop others
            chrom = chrom.replace("chr", "")
            if chrom in ["MT", "M"]:
                continue
            if not (chrom.isdigit() or chrom in ["X", "Y"]):
                continue

            attrs = parse_attrs(attrs_s)
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue
            gene_id = gene_id.split(".")[0]  # remove version
            if gene_id in seen:
                continue
            seen.add(gene_id)

            out.write(f"{gene_id}\t{chrom}\t{start}\t{end}\t{strand}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: make_gene_loc_from_gtf.py Homo_sapiens.GRCh38.115.gtf.gz out.gene.loc", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
