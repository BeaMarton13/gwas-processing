#!/usr/bin/env python3
import pandas as pd
import numpy as np

MAT  = "./wgcna/mat/gene_by_study.nomagma.top8k.z.tsv"
MOD  = "./wgcna/results/modules.simple.tsv"
OUT  = "./wgcna/results/module_summary.tsv"
HUBS = "./wgcna/results/module_hubs.tsv"

MIN_GENES = 5
TOP_HUBS = 10
SKIP_MODULE = 0  # usually "grey"/unassigned

def rowwise_corr_with_vector(A: np.ndarray, v: np.ndarray) -> np.ndarray:
    """
    Correlate each row of A (n_genes x n_traits) with vector v (n_traits,).
    Returns (n_genes,) correlations. Rows with zero variance -> NaN.
    """
    # center
    A0 = A - A.mean(axis=1, keepdims=True)
    v0 = v - v.mean()

    # denom
    A_ss = np.sum(A0 * A0, axis=1)
    v_ss = float(np.sum(v0 * v0))

    denom = np.sqrt(A_ss * v_ss)

    # dot
    num = A0 @ v0

    cor = np.full(A.shape[0], np.nan, dtype=float)
    ok = denom > 0
    cor[ok] = num[ok] / denom[ok]
    return cor

def main():
    X = pd.read_csv(MAT, sep="\t", index_col=0)
    M = pd.read_csv(MOD, sep="\t")

    # Ensure gene ids match as strings
    X.index = X.index.astype(str)
    M["gene"] = M["gene"].astype(str)

    study_cols = list(X.columns)

    summaries = []
    hubs_rows = []

    for mod, genes in M.groupby("module")["gene"]:
        if mod == SKIP_MODULE:
            continue

        g = [x for x in genes.tolist() if x in X.index]
        if len(g) < MIN_GENES:
            continue

        sub_df = X.loc[g]  # genes x traits
        sub = sub_df.to_numpy(dtype=float)

        # Drop genes with zero variance across traits (otherwise corr is undefined)
        row_sd = sub_df.std(axis=1, ddof=0).to_numpy()
        keep = row_sd > 0
        g = [gene for gene, k in zip(g, keep) if k]
        sub = sub[keep, :]

        if sub.shape[0] < MIN_GENES:
            continue

        # Module eigengene across studies = first right singular vector
        sub_centered = sub - sub.mean(axis=0, keepdims=True)
        _, _, Vt = np.linalg.svd(sub_centered, full_matrices=False)
        eig = Vt[0, :]  # length = n_traits

        # kME per gene = correlation with eigengene
        cor = rowwise_corr_with_vector(sub, eig)
        kME_abs = np.abs(cor)

        # Pick top hubs by |kME|
        order = np.argsort(-kME_abs)
        top = order[: min(TOP_HUBS, len(order))]

        for i in top:
            hubs_rows.append({
                "module": int(mod),
                "gene": g[i],
                "kME_abs": float(kME_abs[i])
            })

        # Summary row with eigengene values labeled by study names
        srow = {"module": int(mod), "n_genes": int(len(g))}
        for c, val in zip(study_cols, eig):
            srow[f"eig_{c}"] = float(val)
        summaries.append(srow)

    # Write outputs
    summary_df = pd.DataFrame(summaries).sort_values("n_genes", ascending=False)
    hubs_df = pd.DataFrame(hubs_rows).sort_values(["module", "kME_abs"], ascending=[True, False])

    summary_df.to_csv(OUT, sep="\t", index=False)
    hubs_df.to_csv(HUBS, sep="\t", index=False)

    # Print to screen
    print("\n=== Module summary (top 20 by size) ===")
    if len(summary_df) == 0:
        print("No modules passed filters.")
    else:
        print(summary_df.head(20).to_string(index=False))

    print("\n=== Top hubs per module (first 40 rows) ===")
    if len(hubs_df) == 0:
        print("No hubs computed.")
    else:
        print(hubs_df.head(40).to_string(index=False))

    print("\nWrote:", OUT)
    print("Wrote:", HUBS)

if __name__ == "__main__":
    main()
