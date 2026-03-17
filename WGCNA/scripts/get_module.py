"""
Print all genes assigned to WGCNA module 3 with their names and biotypes.

Loads the module assignment table (modules.simple.tsv) and the Ensembl gene
ID → name/biotype annotation, filters to module 3, and prints the full gene
list sorted alphabetically by gene name.

Intended as a quick drill-down inspection script; output is written to stdout.

Inputs (relative to project root):
    wgcna/results/modules.simple.tsv
    annotation/ensembl115_gene_id_to_name.tsv
"""
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