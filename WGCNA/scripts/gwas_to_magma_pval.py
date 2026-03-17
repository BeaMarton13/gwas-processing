#!/usr/bin/env python3
"""
Convert a raw GWAS summary statistics file to MAGMA SNP p-value format.

Reads a gzipped TSV from the NHGRI-EBI GWAS Catalog, validates required columns,
builds deterministic SNP IDs (CHR:POS:A1:A2), clamps p-values to the range
(1e-300, 1.0], drops invalid/duplicate rows, and writes a gzipped TSV with
columns SNP, P, N ready for MAGMA's --pval input.

Usage:
    gwas_to_magma_pval.py <input.tsv.gz> <output.magma.pval.tsv.gz>
"""
import gzip
import sys
import pandas as pd
import numpy as np

def read_tsv_gz(path: str) -> pd.DataFrame:
    """Read a gzipped tab-separated GWAS summary statistics file.

    Args:
        path: Path to a `.tsv.gz` file.

    Returns:
        pd.DataFrame: Full contents of the file with original column names.
    """
    return pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)

def clamp_p(p: pd.Series) -> pd.Series:
    """Coerce p-values to numeric and clamp to the open interval (1e-300, 1.0].

    MAGMA requires strictly positive p-values not exceeding 1. This function
    converts the series to float (setting non-numeric entries to NaN) and clips
    the range to avoid exact zeros or values above 1.

    Args:
        p: Series of raw p-value strings or floats.

    Returns:
        pd.Series: Numeric p-values clipped to [1e-300, 1.0].
    """
    # MAGMA expects 0<p<=1, avoid exact 0
    p = pd.to_numeric(p, errors="coerce")
    p = p.clip(lower=1e-300, upper=1.0)
    return p

def main(in_gz: str, out_gz: str):
    """Convert a GWAS summary stats file to MAGMA SNP p-value format.

    Reads the input file, validates that all required columns are present,
    standardises types, constructs a deterministic SNP ID of the form
    CHR:POS:EFFECT_ALLELE:OTHER_ALLELE, drops rows with missing or invalid
    values (N <= 0, non-finite coordinates), deduplicates on SNP ID, and
    writes the result as a gzipped TSV.

    Output columns (tab-separated, with header):
        SNP   Variant identifier: CHR:POS:A1:A2
        P     Clamped p-value in (1e-300, 1.0]
        N     Sample size (positive integer)

    Args:
        in_gz:  Path to a gzipped GWAS summary statistics TSV. Must contain
            columns: chromosome, base_pair_location, effect_allele,
            other_allele, p_value, n.
        out_gz: Path for the gzipped output TSV.

    Raises:
        SystemExit: If any required columns are absent from the input file.
    """
    df = read_tsv_gz(in_gz)

    required = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "p_value", "n"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns: {missing}. Found: {list(df.columns)}")

    # Standardize types
    df["chromosome"] = df["chromosome"].astype(str).str.replace("^chr", "", regex=True)
    df["base_pair_location"] = pd.to_numeric(df["base_pair_location"], errors="coerce")
    df["n"] = pd.to_numeric(df["n"], errors="coerce")
    df["P"] = clamp_p(df["p_value"])

    # Drop junk
    df = df.dropna(subset=["chromosome", "base_pair_location", "P", "n"])
    df = df[df["n"] > 0]

    # Build a deterministic SNP ID that you can also enforce in your LD reference
    # Format: CHR:POS:A1:A2  where A1=effect_allele, A2=other_allele
    # (Any nomenclature is fine as long as it matches the reference IDs.)  :contentReference[oaicite:7]{index=7}
    df["SNP"] = (
        df["chromosome"].astype(str)
        + ":"
        + df["base_pair_location"].astype(int).astype(str)
        + ":"
        + df["effect_allele"].astype(str)
        + ":"
        + df["other_allele"].astype(str)
    )

    out = df[["SNP", "P", "n"]].rename(columns={"n": "N"})
    out = out.drop_duplicates(subset=["SNP"])

    out.to_csv(out_gz, sep="\t", index=False, compression="gzip")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: gwas_to_magma_pval.py input.tsv.gz output.magma.pval.tsv.gz", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
