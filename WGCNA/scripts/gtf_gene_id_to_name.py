#!/usr/bin/env python3
import gzip
import re
import sys

ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')

def parse_attrs(s):
    return {m.group(1): m.group(2) for m in ATTR_RE.finditer(s)}

def main(gtf_gz, out_tsv):
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
