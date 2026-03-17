#!/usr/bin/env python3
import glob, os
import pandas as pd
import numpy as np
import pyranges as pr

GTF_LOC = "/magma/ref/ensembl115.GRCh38.gene.loc"
WINDOW = 10_000  # +/-10kb

def load_genes():
    genes = pd.read_csv(GTF_LOC, sep="\t", header=None,
                        names=["gene","chr","start","end","strand"])
    genes["Chromosome"] = genes["chr"].astype(str).str.replace("^chr","", regex=True)
    genes["Start"] = (genes["start"] - 1 - WINDOW).clip(lower=0)
    genes["End"]   = (genes["end"] + WINDOW)
    return pr.PyRanges(genes[["Chromosome","Start","End","gene"]])

GENES = load_genes()

def gene_scores_signed_z(gwas_tsv_gz: str) -> pd.Series:
    g = pd.read_csv(gwas_tsv_gz, sep="\t", compression="gzip", low_memory=False)

    g = g.dropna(subset=[
        "chromosome",
        "base_pair_location",
        "beta",
        "standard_error"
    ])

    g["Chromosome"] = g["chromosome"].astype(str).str.replace("^chr","", regex=True)
    g["Start"] = g["base_pair_location"].astype(int) - 1
    g["End"]   = g["Start"] + 1

    g["beta"] = pd.to_numeric(g["beta"], errors="coerce")
    g["standard_error"] = pd.to_numeric(g["standard_error"], errors="coerce")

    g = g[g["standard_error"] > 0]

    # Signed SNP Z-score
    g["Z"] = g["beta"] / g["standard_error"]

    snps = pr.PyRanges(g[["Chromosome","Start","End","Z"]])
    ov = snps.join(GENES)
    df = ov.df

    if df.empty:
        raise ValueError(f"No SNPs mapped to genes for {gwas_tsv_gz}")

    # Stouffer aggregation
    agg = df.groupby("gene")["Z"].agg(["sum","count"])
    gene_z = agg["sum"] / np.sqrt(agg["count"])

    return gene_z.rename("score")

def main():
    inputs = sorted(glob.glob("./gwas_raw/*.tsv.gz"))
    if len(inputs) < 2:
        raise SystemExit("Put your GWAS .tsv.gz files in ./gwas_raw/")

    mats = []
    for f in inputs:
        print("Processing:", f)
        name = os.path.basename(f).replace(".tsv.gz","")
        s = gene_scores_signed_z(f).rename(name)
        mats.append(s)
        print("Done:", f)

    mat = pd.concat(mats, axis=1, join="inner")

    # Column standardization (important for WGCNA)
    mat = (mat - mat.mean()) / mat.std()

    out = "wgcna/mat/gene_by_study.nomagma.signed_z.tsv"
    mat.to_csv(out, sep="\t")

    print("Wrote:", out)
    print("Shape genes x studies:", mat.shape)

if __name__ == "__main__":
    main()
