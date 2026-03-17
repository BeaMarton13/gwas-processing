#!/usr/bin/env python3
"""
Extract a gene ID → gene name / biotype lookup table from a gzipped Ensembl GTF.

Parses only 'gene' feature rows, strips Ensembl version suffixes from gene IDs
(e.g. ENSG00000123456.3 → ENSG00000123456), deduplicates on gene ID, and writes
a three-column TSV (gene_id, gene_name, gene_biotype).

Usage:
    gtf_gene_id_to_name.py <input.gtf.gz> <output.tsv>
"""
import gzip
import re
import sys

ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')

def parse_attrs(s):
    """Parse a GTF attribute string into a key/value dictionary.

    Extracts all `key "value"` pairs from the semicolon-separated attribute
    field (column 9) of a GTF record using a pre-compiled regex.

    Args:
        s: The raw attribute string from a GTF line, e.g.
           'gene_id "ENSG00000123456.3"; gene_name "BRCA2"; gene_biotype "protein_coding";'

    Returns:
        dict[str, str]: Mapping of attribute keys to their unquoted values.
    """
    return {m.group(1): m.group(2) for m in ATTR_RE.finditer(s)}

def main(gtf_gz, out_tsv):
    """Parse an Ensembl GTF and write a gene ID → name/biotype TSV.

    Streams through the gzipped GTF, processing only 'gene' feature rows.
    Skips comment lines and malformed records. For each gene entry, extracts
    the Ensembl gene ID (version suffix stripped), gene name, and biotype
    (falling back to 'gene_type' if 'gene_biotype' is absent). Duplicate
    gene IDs are silently skipped; the first occurrence is kept.

    Output columns (tab-separated, with header):
        gene_id      Ensembl gene ID without version suffix (e.g. ENSG00000123456)
        gene_name    HGNC gene symbol, or empty string if absent
        gene_biotype Biotype label (e.g. protein_coding, lncRNA), or empty string

    Args:
        gtf_gz:  Path to a gzipped Ensembl GTF file.
        out_tsv: Path for the output TSV file.
    """
    seen = set()
    with gzip.open(gtf_gz, "rt") as f, open(out_tsv, "w") as out:
        out.write("gene_id\tgene_name\tgene_biotype\n")
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            feature = parts[2]
            if feature != "gene":
                continue
            attrs = parse_attrs(parts[8])
            gid = attrs.get("gene_id")
            gname = attrs.get("gene_name")
            biotype = attrs.get("gene_biotype") or attrs.get("gene_type")
            if not gid:
                continue
            gid = gid.split(".")[0]  # strip version
            if gid in seen:
                continue
            seen.add(gid)
            out.write(f"{gid}\t{gname or ''}\t{biotype or ''}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: gtf_gene_id_to_name.py Homo_sapiens.GRCh38.115.gtf.gz out.tsv", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
