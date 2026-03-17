
import pandas as pd
import numpy as np

# ----------------------------
# Settings
# ----------------------------
USE_SIGNED_Z = True      # True = signed Z matrix
Z_THRESHOLD = 5.0        # |Z| >= 3 considered strong
KME_THRESHOLD = 0.3      # optional hub cutoff (set None to disable)

# ----------------------------
# Load inputs
# ----------------------------
mods = pd.read_csv(
    "./wgcna/results/modules.simple.tsv",
    sep="\t"
)

hubs = pd.read_csv(
    "./wgcna/results/module_hubs.tsv",
    sep="\t"
)

names = pd.read_csv(
    "./annotation/ensembl115_gene_id_to_name.tsv",
    sep="\t"
)

gene_scores = pd.read_csv(
    "./wgcna/mat/gene_by_study.nomagma.signed_z.tsv",
    sep="\t",
    index_col=0
)

# Ensure consistent types
mods["gene"] = mods["gene"].astype(str)
hubs["gene"] = hubs["gene"].astype(str)
names["gene_id"] = names["gene_id"].astype(str)
gene_scores.index = gene_scores.index.astype(str)

# ----------------------------
# Compute strongest GWAS signal per gene
# ----------------------------
# For signed Z: use max absolute value
gwas_strength = gene_scores.abs().max(axis=1).rename("max_abs_z")

# ----------------------------
# Process each module
# ----------------------------
for mod in sorted(mods["module"].dropna().unique()):

    mdf = mods[mods["module"] == mod].copy()

    df = (
        mdf
        .merge(gwas_strength, left_on="gene", right_index=True, how="left")
        .merge(hubs, on="gene", how="left")
        .merge(
            names[["gene_id", "gene_name", "gene_biotype"]],
            left_on="gene",
            right_on="gene_id",
            how="left"
        )
    )

    # ----------------------------
    # Filtering
    # ----------------------------
    filt = (
        df["gene_name"].notna() &
        df["max_abs_z"].notna() &
        (df["max_abs_z"] >= Z_THRESHOLD)
    )

    # if KME_THRESHOLD is not None:
    #     filt &= df["kME_abs"].fillna(0) >= KME_THRESHOLD

    df_filt = df[filt].copy()

    # ----------------------------
    # Final formatting
    # ----------------------------
    out = df_filt[[
        "gene",
        "gene_name",
        "gene_biotype",
        "max_abs_z",
        "kME_abs"
    ]].copy()

    # Sort: strongest hubs first, then strongest GWAS
    out["kME_abs_sort"] = out["kME_abs"].fillna(-np.inf)
    out = (
        out
        # .sort_values(["kME_abs_sort", "max_abs_z"], ascending=[False, False])
        # .drop(columns="kME_abs_sort")
    )

    # ----------------------------
    # PRINT CLEANLY TO SCREEN
    # ----------------------------
    print("\n" + "=" * 70)
    print(f"Module {mod} — HIGH-CONFIDENCE genes (n = {len(out)})")
    print("=" * 70)

    if len(out) == 0:
        print("No genes passed filters.")
    else:
        print(out.to_string(index=False))

    # ----------------------------
    # Save file
    # ----------------------------
    outfile = f"./wgcna/results/module_{mod}_high_confidence_genes.tsv"
    out.to_csv(outfile, sep="\t", index=False)

print("\nDone.")
