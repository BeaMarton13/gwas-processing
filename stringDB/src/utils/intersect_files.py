import pandas as pd
import argparse
import pyranges as pr
import numpy as np
import re


# -----------------------------
# Utility to extract gene name
# -----------------------------
def extract_gene_name(attr):
    match = re.search(r'gene_name "([^"]+)"', attr)
    return match.group(1) if match else None

def select_min_pvalue_rows_multi(df, gene_col="Gene_0"):
    """
    For each gene, select the row(s) with the minimum p-value across all p_value_* columns.

    Parameters:
        df : pandas DataFrame
            Must contain gene_col and one or more p_value_* columns
        gene_col : str
            Column name containing gene identifiers

    Returns:
        pandas DataFrame with row(s) corresponding to the minimum p-value per gene
    """

    # Drop rows with NaN genes
    df_clean = df.dropna(subset=[gene_col])

    # Identify all p-value columns
    pval_cols = [c for c in df_clean.columns if c.startswith("p_value")]

    # Compute minimum p-value across all pval_cols for each row
    df_clean = df.dropna(subset=[gene_col]).copy()  # make an explicit copy

    # Compute row-wise minimum safely
    df_clean.loc[:, "_row_min_p"] = df_clean[pval_cols].min(axis=1, skipna=True)

    # Compute minimum p-value per gene
    min_per_gene = df_clean.groupby(gene_col)["_row_min_p"].min().reset_index()
    min_per_gene = min_per_gene.rename(columns={"_row_min_p": "_gene_min_p"})

    # Merge to keep only rows with minimum p-value for each gene
    df_min = df_clean.merge(min_per_gene, on=gene_col)
    df_min = df_min[df_min["_row_min_p"] == df_min["_gene_min_p"]].drop(
        columns=["_row_min_p", "_gene_min_p"]
    )

    return df_min.reset_index(drop=True)


# -----------------------------
# Read file in any supported format
# -----------------------------
def read_file(file_w_path):
    print(f'Reading file: {file_w_path}')
    try:
        if file_w_path.endswith('.csv'):
            return pd.read_csv(file_w_path)

        elif file_w_path.endswith('.xlsx'):
            return pd.read_excel(file_w_path)

        elif file_w_path.endswith('.bed'):
            try:
                df = pd.read_csv(
                    file_w_path,
                    sep="\t",
                    header=None,
                    usecols=[0, 1, 2, 3],
                    names=["Chromosome", "Start", "End", "p_value"],
                    dtype={"Chromosome": str, "Start": int, "End": int, "p_value": float}
                )
            except:
                # df = pd.read_csv(
                #     file_w_path,
                #     sep="\t",
                #     header=None,
                #     usecols=[0, 1, 2, 7, 9],
                #     names=["Chromosome", "Start", "End", "Feature", "Attributes"],
                # )
                # df["p_value"] = np.nan
                # df = df[df["Feature"] == "gene"]
                # df["Gene"] = df["Attributes"].apply(extract_gene_name)
                df = pd.read_csv(
                    file_w_path,
                    sep="\t",
                    header=None,
                    usecols=[0, 1, 2, 3, 5],
                    names=["Chromosome", "Start", "End", "GeneInfo", "Strand"],
                )

                df["Gene"] = df["GeneInfo"].str.split("|").str[0]
            print(file_w_path, "loaded with shape:", df.shape)
            print(df)
            return df

        else:
            raise ValueError("Unsupported file format")

    except Exception as e:
        print(f"Error reading file: {e}")
        return None


# ----------------------------------------------------------
# Merge files WITHOUT producing tuples, one gene per row
# ----------------------------------------------------------
def merge_files(dfs):
    gr_list = [pr.PyRanges(df) for df in dfs]

    # intersection regions only
    intersection = gr_list[0]
    for idx, gr in enumerate(gr_list[1:]):
        intersection = intersection.intersect(gr)
        if intersection.df.empty:
            print(f"⚠️ No overlapping regions found among the provided files at file {idx + 2}.")
            exit()
        
    base = intersection.df[["Chromosome", "Start", "End"]].drop_duplicates()

    out = base.copy()
    gene_candidates = ["Gene", "gene_name", "gene", "gene_id", "Gene_name"]

    for i, df in enumerate(dfs):
        gr = pr.PyRanges(df)
        joined = pr.PyRanges(base).join(gr)
        jdf = joined.df

        pcol = f"p_value_{i}"
        gcol = f"Gene_{i}"

        # p-value
        if "p_value" in jdf.columns:
            pv = jdf[["Chromosome", "Start", "End", "p_value"]].copy()
            pv = pv.dropna(subset=["p_value"])
            pv = pv.rename(columns={"p_value": pcol})
            out = out.merge(pv, on=["Chromosome", "Start", "End"], how="left")
        else:
            out[pcol] = np.nan

        # gene
        gene_col = next((c for c in gene_candidates if c in jdf.columns), None)
        if gene_col:
            g = jdf[["Chromosome", "Start", "End", gene_col]].dropna()
            g = g.rename(columns={gene_col: gcol})
            out = out.merge(g, on=["Chromosome", "Start", "End"], how="left")
        else:
            out[gcol] = np.nan

    # CLEAN STEP: split comma-separated gene strings into multiple rows
    for c in out.columns:
        if c.startswith("Gene"):
            out[c] = out[c].apply(
                lambda x: [g.strip() for g in x.split(",")] if isinstance(x, str) and "," in x else x
            )
            out = out.explode(c)

    return out.reset_index(drop=True)


# ----------------------------------------------------------
# MAIN
# ----------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Intersect files.")
    parser.add_argument("disease", type=str)
    parser.add_argument("files", nargs="+", type=str)
    args = parser.parse_args()

    dfs = [read_file(f) for f in args.files]
    res_df = merge_files(dfs)

    print(res_df.head(50))

    # -------------------------
    # Save unique gene list
    # -------------------------
    gene_col = "Gene_0"
    gene_list_path = f"./data/gene_names/{args.disease}.txt"

    res_df[gene_col].dropna().drop_duplicates().to_csv(
        gene_list_path, index=False, header=True
    )

    # -------------------------
    # Compute p-value stats
    # -------------------------
    pval_cols = [c for c in res_df.columns if c.startswith("p_value")]

    # ensure all numeric
    for col in pval_cols:
        res_df[col] = pd.to_numeric(res_df[col], errors="coerce")

    res_df = select_min_pvalue_rows_multi(res_df, gene_col=gene_col)

    avg_pvalues = res_df.groupby(gene_col)[pval_cols].mean().reset_index()
    avg_pvalues["mean_pvalue_all"] = avg_pvalues[pval_cols].mean(axis=1)
    avg_pvalues["min_pvalue_all"] = avg_pvalues[pval_cols].min(axis=1)
    avg_pvalues["median_pvalue_all"] = avg_pvalues[pval_cols].median(axis=1)

    out_file = f"./data/gene_pvalues/{args.disease}.csv"
    avg_pvalues.to_csv(out_file, index=False)

    print("Saved:", out_file)
