import pandas as pd

mods = pd.read_csv(
    "./wgcna/results/modules.simple.tsv",
    sep="\t"
)

names = pd.read_csv(
    "./annotation/ensembl115_gene_id_to_name.tsv",
    sep="\t"
)

df = (
    mods[mods["module"] == 3]
    .merge(names, left_on="gene", right_on="gene_id", how="left")
    [["gene", "gene_name", "gene_biotype"]]
    .sort_values("gene_name", na_position="last")
)

print(f"Module 3: n = {len(df)} genes\n")
print(df.to_string(index=False))