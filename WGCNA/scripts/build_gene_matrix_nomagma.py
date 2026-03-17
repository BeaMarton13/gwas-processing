#!/usr/bin/env python3
"""
Build a gene × study signal matrix from GWAS summary statistics without MAGMA.

For each GWAS study (*.tsv.gz in ./gwas_raw/), SNPs are mapped to genes using a
±WINDOW bp interval around each Ensembl gene body. Per-gene signal is aggregated
using Stouffer's signed Z-score method (Z = beta / SE, then summed and normalised
by sqrt(N_snps)). The resulting matrix is column-standardised and written to
wgcna/mat/gene_by_study.nomagma.signed_z.tsv for use as WGCNA input.
"""
import glob, os
import pandas as pd
import numpy as np
import pyranges as pr

GTF_LOC = "/magma/ref/ensembl115.GRCh38.gene.loc"
WINDOW = 10_000  # +/-10kb

def load_genes():
    """Load gene coordinates from the MAGMA gene location file and expand by WINDOW.

    Reads the tab-separated gene location file (columns: gene, chr, start, end, strand),
    strips any 'chr' prefix from chromosome names, and expands each gene interval by
    WINDOW bp upstream and downstream to capture nearby regulatory SNPs. Intervals are
    clipped at 0 to avoid negative start positions.

    Returns:
        pr.PyRanges: PyRanges object with columns Chromosome, Start, End, gene.
    """
    genes = pd.read_csv(GTF_LOC, sep="\t", header=None,
                        names=["gene","chr","start","end","strand"])
    genes["Chromosome"] = genes["chr"].astype(str).str.replace("^chr","", regex=True)
    genes["Start"] = (genes["start"] - 1 - WINDOW).clip(lower=0)
    genes["End"]   = (genes["end"] + WINDOW)
    return pr.PyRanges(genes[["Chromosome","Start","End","gene"]])

GENES = load_genes()

def gene_scores_signed_z(gwas_tsv_gz: str) -> pd.Series:
    """Compute a signed Stouffer Z-score per gene from a GWAS summary statistics file.

    Reads the gzipped TSV, drops rows with missing chromosome, position, beta, or
    standard error values, and computes a per-SNP Z-score (Z = beta / SE). SNPs are
    then mapped to genes using a genomic interval join against the pre-loaded GENES
    PyRanges (gene bodies ± WINDOW bp). For each gene, all overlapping SNP Z-scores
    are aggregated via Stouffer's method:

        gene_Z = sum(Z_snps) / sqrt(N_snps)

    Args:
        gwas_tsv_gz: Path to a gzipped GWAS summary statistics TSV file. Expected
            columns include 'chromosome', 'base_pair_location', 'beta', and
            'standard_error'.

    Returns:
        pd.Series: Series indexed by Ensembl gene ID with Stouffer Z-scores,
            named 'score'.

    Raises:
        ValueError: If no SNPs can be mapped to any gene after the interval join.
    """
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
    """Build and save the gene × study signed Z-score matrix.

    Discovers all GWAS summary statistics files matching ./gwas_raw/*.tsv.gz,
    computes a signed Stouffer Z-score vector per study via gene_scores_signed_z(),
    and joins them into a matrix on the inner set of genes (genes present in all
    studies). The matrix is column-standardised (zero mean, unit variance) before
    being written to wgcna/mat/gene_by_study.nomagma.signed_z.tsv.

    Raises:
        SystemExit: If fewer than two input files are found in ./gwas_raw/.
    """
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
