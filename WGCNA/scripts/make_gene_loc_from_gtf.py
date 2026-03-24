#!/usr/bin/env python3
"""
Build a MAGMA-format gene location file from a gzipped Ensembl GTF.

Parses only ‘gene’ feature rows, filters to canonical chromosomes (1–22, X, Y),
strips Ensembl version suffixes, deduplicates on gene ID, and writes a
tab-separated file with columns: GENEID, CHR, START, END, STRAND.

Usage:
    make_gene_loc_from_gtf.py <input.gtf.gz> <output.gene.loc>
"""
import gzip
import re
import sys

def parse_attrs(attr_str: str) -> dict:
    """Parse a GTF attribute string into a key/value dictionary.

    Extracts all `key "value"` pairs from the semicolon-separated attribute
    field (column 9) of a GTF record.

    Args:
        attr_str: The raw attribute string from a GTF line, e.g.
            ‘gene_id "ENSG00000123456.3"; gene_name "BRCA2"; gene_biotype "protein_coding";’

    Returns:
        dict[str, str]: Mapping of attribute keys to their unquoted values.
    """
    attrs = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_str):
        attrs[m.group(1)] = m.group(2)
    return attrs

def main(gtf_gz: str, out_loc: str):
    """Parse an Ensembl GTF and write a MAGMA gene location file.

    Streams through the gzipped GTF and processes only ‘gene’ feature rows.
    Skips comment lines, malformed records, mitochondrial contigs (MT/M), and
    any non-standard chromosomes (i.e. not 1–22, X, or Y). Strips Ensembl
    version suffixes from gene IDs and silently skips duplicate gene IDs.

    Output columns (tab-separated, no header):
        GENEID   Ensembl gene ID without version suffix
        CHR      Chromosome (e.g. 1, X)
        START    1-based start coordinate (from GTF)
        END      End coordinate (from GTF)
        STRAND   + or -

    Args:
        gtf_gz:  Path to a gzipped Ensembl GTF file.
        out_loc: Path for the output gene location file.
    """
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
