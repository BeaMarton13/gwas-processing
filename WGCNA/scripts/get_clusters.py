import pandas as pd

mods = pd.read_csv(
    "project/wgcna/results/modules.simple.tsv",
    sep="\t"
)

names = pd.read_csv(
    "project/annotation/ensembl115_gene_id_to_name.tsv",
    sep="\t"
)

selected_genes = ['NECTIN2', 'LNCOB1', 'TOMM40', 'APOE', 'APOC1' ]

modules = mods["module"].unique()
for mod in modules:

    print("==============================================")

    df = (
        mods[mods["module"] == mod]
        .merge(names, left_on="gene", right_on="gene_id", how="left")
        [["gene", "gene_name", "gene_biotype"]]
        .dropna(subset=["gene_name"])
        .sort_values("gene_name", na_position="last")
    )
    if any(g in df["gene_name"].values for g in selected_genes):
        print('Genes found in module:', mod, [g for g in selected_genes if g in df["gene_name"].values])
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ', selected_genes)

    print(f"Module {mod}: n = {len(df)} genes\n")
    print(df.to_string(index=False))