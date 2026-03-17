#!/usr/bin/env python3
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

IN = "./wgcna/mat/gene_by_study.nomagma.top8k.z.tsv"
OUT = "./wgcna/results/modules.simple.tsv"

POWER = 6          # try 4,6,8
MIN_MODULE = 30    # minimum genes in a module
CORR_TYPE = "pearson"  # pearson is fine with 4 columns

def main():
    X = pd.read_csv(IN, sep="\t", index_col=0)  # genes x studies
    genes = X.index

    # correlation across studies (columns) => gene-gene correlation
    C = np.corrcoef(X.values)  # genes x genes
    C = np.nan_to_num(C, nan=0.0)

    # unsigned adjacency
    A = np.abs(C) ** POWER
    np.fill_diagonal(A, 0)

    # distance for clustering
    D = 1 - A
    D = np.clip(D, 0, 1)
    Z = linkage(squareform(D, checks=False), method="average")

    # choose a cut height (tune)
    # lower cut height => more modules, higher => fewer
    cut_height = 0.75
    labels = fcluster(Z, t=cut_height, criterion="distance")

    # enforce min module size
    lab_series = pd.Series(labels, index=genes)
    counts = lab_series.value_counts()
    small = counts[counts < MIN_MODULE].index
    lab_series.loc[lab_series.isin(small)] = 0  # 0 = "grey/unassigned"

    out = pd.DataFrame({"gene": genes, "module": lab_series.values})
    out.to_csv(OUT, sep="\t", index=False)
    print("Wrote:", OUT)
    print("Modules (excluding 0):", (out["module"]!=0).sum(), "genes assigned")

if __name__ == "__main__":
    main()
