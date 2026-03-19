#!/usr/bin/env python3
"""
Select the top N genes by variance for WGCNA clustering.

Reads the full gene × study signed Z-score matrix, computes row-wise variance
across studies, selects the top N genes with highest variance, and writes a
filtered matrix for input to WGCNA clustering.

Usage:
    python select_top_genes.py

Inputs:
    wgcna/mat/gene_by_study.nomagma.signed_z.tsv

Outputs:
    wgcna/mat/gene_by_study.nomagma.top8k.z.tsv
"""
import pandas as pd
import numpy as np

# Configuration
IN_FILE = "wgcna/mat/gene_by_study.nomagma.signed_z.tsv"
OUT_FILE = "wgcna/mat/gene_by_study.nomagma.top8k.z.tsv"
TOP_N = 8000


def main():
    """Select top N genes by variance and write filtered matrix.

    Loads the full gene × study Z-score matrix, computes the variance across
    studies for each gene (row-wise), ranks genes by descending variance, and
    retains only the top TOP_N genes. The filtered matrix is written to OUT_FILE.

    Outputs:
        Writes a TSV file with TOP_N rows (genes) × K columns (studies).
        Prints summary statistics to stdout.
    """
    print(f"Loading matrix from {IN_FILE}...")
    mat = pd.read_csv(IN_FILE, sep="\t", index_col=0)

    print(f"Original matrix shape: {mat.shape} (genes × studies)")

    # Compute row-wise variance across studies
    gene_var = mat.var(axis=1)

    print(f"Variance range: [{gene_var.min():.4f}, {gene_var.max():.4f}]")

    # Select top N genes by variance
    top_genes = gene_var.nlargest(TOP_N).index
    mat_top = mat.loc[top_genes]

    print(f"Selected top {TOP_N} genes by variance")
    print(f"Filtered matrix shape: {mat_top.shape}")

    # Write output
    mat_top.to_csv(OUT_FILE, sep="\t")
    print(f"Wrote: {OUT_FILE}")

    # Summary statistics
    print(f"\nVariance threshold: >= {gene_var[top_genes].min():.4f}")
    print(f"Top gene variance: {gene_var[top_genes].max():.4f}")


if __name__ == "__main__":
    main()
