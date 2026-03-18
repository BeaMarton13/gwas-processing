Clone the repository. 

<details>
<summary>GWAS-to-Gene-Network Pipeline (stringDB)</summary>


A pipeline for integrating multiple GWAS summary statistics into protein interaction networks using STRING and OmniPath databases, with community detection and network analysis.

The biological focus is the four dementia subtypes (Alzheimer's, vascular, other, and unspecified — GCST90473236, GCST90473240, GCST90473241, GCST90473242), with key locus genes APOE, APOC1, TOMM40, NECTIN2, LNCOB1.

---

## Overview

```
GWAS Catalog (data.xlsx)
    ↓ Look up GCST IDs by disease keyword
Download *.tsv.gz (raw GWAS summary stats)
    ↓ Filter by p-value
processed/*.json → *.bed (SNP positions)
    ↓ bedtools intersect with genome annotation
genes_with_pvalues.bed / bed_comparison/*.bed
    ↓ Extract gene names
gene_names/*.txt
    ↓ Query STRING
gene_networks/*.tsv → *.gml
    ↓ Community detection (InfoMap / Voronoi)
Cluster statistics + Plotly visualizations
```

---

## Project Structure

```
stringDB/
├── main.py                         # Entry point: download + p-value filter
├── process_disease.sh              # Orchestration script (download → BED)
├── compare_bed_files.sh            # bedtools intersect wrapper (two files)
├── intersect_beds.sh               # bedtools intersect wrapper (N files)
├── genes_with_pvalues.bed          # Combined SNP-GTF intersection BED
├── *_diseases_w_id.txt             # Auto-generated GCST ID manifests per disease
│
├── src/utils/
│   ├── downloader.py               # EBI FTP downloader
│   ├── excel_handler.py            # GWAS Catalog Excel query
│   ├── process_file.py             # Chunked TSV reader + GTF reader
│   ├── convert_json_to_bed.py      # JSON → BED conversion
│   ├── plotting.py                 # Keyword bubble chart
│   ├── intersect_files.py          # pyranges-based genomic intersection + gene p-values
│   ├── process_pathways.py         # ShinYGO GO enrichment analysis
│   ├── handle_network.py           # Network clustering, centrality, visualization
│   ├── export_network_from_tsv.py  # STRING TSV → igraph GML
│   └── voronoi_c.so             # Compiled C library for Voronoi community detection
│
├── data/
│   ├── <disease>/                  # Per-disease GWAS data
│   │   ├── *.tsv.gz                # Raw GWAS summary stats (~3.2 GB/study)
│   │   └── processed*/             # Filtered JSON + BED files per p-value threshold
│   ├── homo_sapiens/               # Ensembl GRCh38 v115 GTF + BED annotation
│   ├── gene_names/                 # Per-disease gene lists (CSV header + one gene per line)
│   ├── gene_pvalues/               # Per-gene aggregated p-value statistics (CSV)
│
├── bed_comparison/                 # bedtools intersection outputs
├── gene_names/                     # Root-level gene name extractions
├── gene_networks/                  # STRING TSVs and GML network files
```

---

## Pipeline Steps

### Step 1 — Download GWAS Summary Statistics

```bash
python3 main.py <disease> <processed_folder> <p_value_limit> <whole_word_match>
# e.g.: python3 main.py dementia processed_0_01 0.01 True
```

Or use the orchestration script (runs Steps 1–2 together):

```bash
./process_disease.sh <disease> <processed_folder> <p_value_limit> <whole_word_match>
```

`main.py` queries `data.xlsx` (GWAS Catalog spreadsheet) for matching GCST IDs using the
`"Reported trait"` column, downloads each `*.tsv.gz` from the EBI FTP server
(`src/utils/downloader.py`), and writes a `<disease>_diseases_w_id.txt` manifest.

It then reads each GWAS file in chunks of 1000 rows and saves rows where
`p_value < P_VALUE_LIMIT` to `data/<disease>/<processed_folder>/<GCST_ID>.json`.

**EBI FTP URL pattern:**
```
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST<bucket>/GCST<id>/GCST<id>.tsv.gz
```

---

### Step 2 — Convert Filtered SNPs to BED Format

```bash
python3 src/utils/convert_json_to_bed.py <disease> <processed_folder>
```

Reads each `.json` in `data/<disease>/<processed_folder>/` and writes a 5-column BED file
(0-based coordinates):

```
chromosome  start  end  p_value  .
```

---

### Step 3 — Convert Genome Annotation to BED (one-time setup)

```bash
gtf2bed < data/homo_sapiens/homo_sapiens.gtf > data/homo_sapiens/homo_sapiens.bed
```

Converts the Ensembl GRCh38 v115 GTF into BED format for genomic intersections.
Requires the `bedops` toolkit.

---

### Step 4 — Intersect SNPs with Gene Annotation

**Python route (recommended, uses `pyranges`):**
```bash
python3 src/utils/intersect_files.py <disease> \
    data/<disease>/processed/GCST1.bed \
    data/<disease>/processed/GCST2.bed \
    [...]
```

Generates:
- `data/gene_names/<disease>.txt` — unique gene symbols found in the intersection
- `data/gene_pvalues/<disease>.csv` — per-gene min/mean/median p-values across studies

**Shell route (two files):**
```bash
./compare_bed_files.sh data/<disease>/processed/<study>.bed \
                       data/homo_sapiens/homo_sapiens.bed \
                       bed_comparison/output.bed
```

**Shell route (N files):**
```bash
./intersect_beds.sh bed_comparison/output.bed \
    data/<disease>/processed/GCST1.bed \
    data/<disease>/processed/GCST2.bed [...]
```

Both shell scripts wrap `bedtools intersect -wa -wb`, producing wide BED records pairing
each SNP with its overlapping gene annotation.

---

### Step 5 — Extract Gene Names (shell route)

When using the shell intersection route, extract gene names from the BED output:

```bash
cat bed_comparison/common_alzheimer.bed \
  | cut -f10 | uniq \
  | cut -d";" -f3 | grep "gene_name" \
  | sed -E 's/.*gene_name "([^"]+)".*/\1/' \
  > gene_names/alzheimer.txt
```

The column index (`-f10` vs `-f9`) may differ depending on the input BED type — see
`notes.txt` for variants.

---

### Step 6 — Build Protein Interaction Network

Query STRING or OmniPath for the gene list, save the result as a TSV, then export to GML.

**STRING (downloaded TSV):**
```bash
python3 src/utils/export_network_from_tsv.py
```
Reads a STRING TSV (`#node1 node2 ... combined_score`), builds a directed `igraph` graph
weighted by `combined_score`, runs InfoMap and Voronoi community detection, and exports
a `.gml` file to `gene_networks/`.


---

### Step 7 — Community Detection and Visualization

```bash
python3 src/utils/handle_network.py
```

Loads a `.gml` network and a gene p-value CSV. Runs two community detection methods:

- **InfoMap** — information-flow-based (`igraph.Graph.community_infomap()`)
- **Voronoi** — custom algorithm in `src/utils/voronoi_c.so` (C library via `ctypes`)

For each cluster, computes an **energy score**: `centrality × −log10(p_value)` per gene,
then aggregates to cluster-level mean/median/min statistics.

Generates an interactive Plotly network visualization:
- Node color = cluster's p-value statistic
- Node border = cluster identity
- Key genes (default: NECTIN2, TOMM40, APOE, APOC1) highlighted in red with larger markers

---


## Key Results (Dementia)

**Intersection genes** (common across all 4 dementia subtypes at threshold 0.00001):
APOE, APOC1, NECTIN2, TOMM40, LNCOB1, plus CEACAM/IGSF/PVR locus genes on chr19.

**Community detection** (InfoMap on `dementia001.gml`): 29 clusters.

---

## Dependencies

- **Python ≥ 3.9** with:
  - `pandas`, `numpy` — data handling
  - `pyranges` — genomic interval intersection
  - `igraph` — network construction and community detection
  - `tqdm` — download progress bars
  - `matplotlib`, `plotly`, `adjust_text` — visualization
  - `gtfparse` — GTF file parsing
- **bedops** — `gtf2bed` for GTF-to-BED conversion
- **bedtools** — `bedtools intersect` / `bedtools subtract`

---

## Data Sources

| Data | Source |
|---|---|
| GWAS summary statistics | [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/) FTP |
| Genome annotation | Ensembl GRCh38 v115 GTF |
| Protein interactions | [STRING](https://string-db.org/) — TSV download for gene list |
| Reference paper | [Nature (2025)](https://doi.org/10.1038/s41586-025-09272-9) |


</details>


On MacOS (Silicone) environment
```bash
cd WGCNA
conda create -n wgcna                                          
conda activate wgcna
conda config --env --set subdir osx-64
conda env update -f environment.yml
conda install bioconda::plink2
```

On other enviromnets
```bash
cd WGCNA
conda env create -f environment.yml -n wgcna
conda install bioconda::plink2
```

<details>
<summary>GWAS-to-WGCNA Integrative Genomics Pipeline</summary>

A pipeline for integrating multiple GWAS summary statistics into gene co-expression modules using a WGCNA-style soft-threshold network approach, without requiring MAGMA gene analysis.

## Overview

This pipeline maps SNP-level GWAS signals to genes via genomic windows, aggregates them into a gene × study signal matrix using Stouffer's signed Z-score method, and clusters genes into co-signal modules using hierarchical network analysis. 

Four GWAS studies from the NHGRI-EBI GWAS Catalog are used as input (GCST90473236, GCST90473240, GCST90473241, GCST90473242), all mapped to the human GRCh38 genome with Ensembl v115 annotations.

---

## Project Structure

```
project/
├── annotation/                        # Gene ID → name/biotype lookup
├── gtf/                               # Ensembl v115 GTF annotation
├── gwas_clean/                        # MAGMA-format cleaned GWAS files
├── gwas_raw/                          # Raw GWAS summary statistics
├── ld_vcf/                            # 1000 Genomes GRCh38 phased VCFs (chr1–22)
├── magma/
│   ├── annot/                         # MAGMA annotation output (unused)
│   ├── genes/                         # MAGMA gene-level output (unused)
│   ├── genesets/                      # MAGMA gene-set output (unused)
│   └── ref/                           # LD reference panel (PLINK BED + pgen)
├── scripts/                           # All analysis scripts
└── wgcna/
    ├── mat/                           # Gene × study signal matrices
    ├── modules/                       # (alternative output dir, empty)
    └── results/                       # Module assignments, hubs, and filtered gene lists
```

---

## Pipeline Steps

### Step 1 — Download LD Reference Data

```bash
bash scripts/download_1000G.sh
```

Downloads 1000 Genomes Project phased biallelic SNV+INDEL VCF files (GRCh38, 2019-03-12 release) for chromosomes 1–22 into `ld_vcf/`. Uses resumable downloads (`wget -c`) and skips already-present files.

---

### Step 2 — Build LD Reference Panel

```bash
bash scripts/generate_chrs.sh
```

Converts each per-chromosome VCF into PLINK1 BED format using `plink2`. Filters to biallelic, ACGT-only SNPs. Output: `magma/ref/plink_chr/chr{1..22}.{bed,bim,fam}`.

---

### Step 3 — Build Gene Annotations

Create a `gtf` folder and download the [annotation data available for human](https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz) and place it into the previously created `gtf` folder.

```bash
mkdir annotation
python scripts/gtf_gene_id_to_name.py \
  gtf/Homo_sapiens.GRCh38.115.gtf.gz \
  annotation/ensembl115_gene_id_to_name.tsv
python scripts/make_gene_loc_from_gtf.py \
  gtf/Homo_sapiens.GRCh38.115.gtf.gz \ 
  magma/ref/ensembl115.GRCh38.gene.loc 
```

- `gtf_gene_id_to_name.py` — Parses `gtf/Homo_sapiens.GRCh38.115.gtf.gz` and outputs `annotation/ensembl115_gene_id_to_name.tsv` (Ensembl ID → gene symbol + biotype).
- `make_gene_loc_from_gtf.py` — Produces `magma/ref/ensembl115.GRCh38.gene.loc` in MAGMA format (GENEID, CHR, START, END, STRAND) for chromosomes 1–22, X, Y.

---

### Step 4 — Preprocess GWAS Summary Statistics

```bash
pip install pandas
mkdir gwas_clean
bash scripts/run_gwas_to_magma.sh
```

Loops over all `gwas_raw/*.tsv.gz` files and calls `scripts/gwas_to_magma_pval.py` on each. Outputs MAGMA-ready SNP files to `gwas_clean/` with columns `SNP | P | N`, where SNP IDs are `CHR:POS:A1:A2`.

Key operations in `gwas_to_magma_pval.py`:
- Strips `chr` prefixes from chromosome fields
- Clamps p-values to [1e-300, 1.0]
- Deduplicates on SNP ID

> **Note:** The `magma/annot/`, `magma/genes/`, and `magma/genesets/` directories are empty — the MAGMA gene analysis step has not been run. The MAGMA-free route below is used instead.

---

### Step 5 — Build Gene × Study Signal Matrix

```bash
python scripts/build_gene_matrix_nomagma.py
```

Maps SNPs to genes using `pyranges` with a ±10 kb window around each Ensembl gene body. For each gene in each study, computes a **signed Stouffer Z-score**:

```
Z_snp  = beta / SE
Z_gene = sum(Z_snps) / sqrt(N_snps)
```

The resulting matrix is column-standardised. Output:
- `wgcna/mat/gene_by_study.nomagma.signed_z.tsv` — 74,940 genes × 4 studies
- `wgcna/mat/gene_by_study.nomagma.minp_log10.tsv` — alternative min(−log10(p)) scoring

Then filter to the top 8,000 genes by variance for input to WGCNA:
```bash
python scripts/filter_genes.py
```
Output: `wgcna/mat/gene_by_study.nomagma.top8k.z.tsv`

---

### Step 6 — WGCNA Clustering

```bash
python scripts/wgcna_simple.py
```

Implements WGCNA-style network analysis in pure Python using `scipy`:

1. Computes all-vs-all Pearson correlations across 4 GWAS studies
2. Builds unsigned soft-threshold adjacency: `A = |r|^6` (power = 6)
3. Converts to a distance matrix: `D = 1 − A`
4. Performs hierarchical clustering (average linkage)
5. Cuts the dendrogram at height 0.75
6. Enforces minimum module size of 30 genes; smaller modules → module 0 ("grey")

Output: `wgcna/results/modules.simple.tsv` (gene → module integer).

**Resulting modules:**

| Module | Genes | Key signal |
|--------|-------|-----------|
| 3      | 6,116 | Uniformly negative across all 4 studies; APOE locus (APOE, TOMM40, APOC1), MAPT/17q21 locus |
| 6      | 715   | APOE-flanking locus; NECTIN2, LNCOB1 |
| 4      | 642   | Driven predominantly by study GCST90473240 |
| 5      | 236   | Mixed signal; CEACAM16-AS1, CCSER1, ANK2 |
| 2      | 204   | TYW1, MAPT-IT1 |
| 1      | 87    | RSRC1 |

---

### Step 7 — Module Summary and Hub Genes

```bash
python scripts/module_summary.py
```

For each module:
- Computes the **module eigengene** (first right singular vector from SVD)
- Computes **kME** (module membership) = Pearson correlation of each gene with the eigengene
- Identifies top 10 hub genes by |kME|

Outputs:
- `wgcna/results/module_summary.tsv` — module sizes + eigengene loadings per study
- `wgcna/results/module_hubs.tsv` — top 10 hub genes per module

---

### Step 8 — Inspect and Filter Results

```bash
python scripts/get_clusters.py   # Print all module memberships; search for candidate genes
python scripts/get_module.py     # Print module 3 genes in detail
python scripts/filter_genes.py   # Filter to high-confidence candidates
```

`filter_genes.py` retains genes that:
- Have a known gene symbol (after annotation join)
- Have a maximum |signed Z-score| ≥ 5.0 across all 4 studies

Outputs: `wgcna/results/module_{1..6}_high_confidence_genes.tsv`

**High-confidence hits (|Z| ≥ 5.0):**

| Module | Count | Top genes |
|--------|-------|-----------|
| 3 | ~94 | LINC02210-CRHR1 (Z=13.2), TOMM40 (Z=13.1), MAPT (Z=13.1), KANSL1 (Z=11.9), APOC1 (Z=10.3) |
| 6 | 13 | LNCOB1 (Z=9.9), NECTIN2 (Z=8.4), CADM2 (Z=7.5), MRTFB (Z=7.5) |
| 4 | 8 | MYRIP (Z=7.8), FMO1-AS1 (Z=6.7) |
| 5 | 4 | CEACAM16-AS1 (Z=7.8) |
| 2 | 2 | TYW1, MAPT-IT1 |
| 1 | 1 | RSRC1 |

A prioritised ranking for module 3 (by kME + GWAS signal) is in `wgcna/results/module_3_prioritised_genes.tsv`.

---

## Key Parameters

| Parameter | Value | Location |
|-----------|-------|----------|
| Gene window (SNP mapping) | ±10 kb | `build_gene_matrix_nomagma.py` |
| Soft-threshold power | 6 | `wgcna_simple.py` |
| Dendrogram cut height | 0.75 | `wgcna_simple.py` |
| Minimum module size | 30 genes | `wgcna_simple.py` |
| Genes used for WGCNA | top 8,000 by variance | `filter_genes.py` |
| High-confidence Z threshold | 5.0 | `filter_genes.py` |

---

## Dependencies

- `plink2` — reference panel construction
- `python` ≥ 3.9 with:
  - `pandas`, `numpy`, `scipy` — core data processing and clustering
  - `pyranges` — genomic interval overlaps (SNP → gene mapping)

---

## Data Sources

| Data | Source |
|------|--------|
| GWAS summary statistics | [NHGRI-EBI GWAS Catalog – Summary Statistics (FTP repository)](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/) — accessions GCST90473236, GCST90473240, GCST90473241, GCST90473242 |
| LD reference panel | 1000 Genomes Project, GRCh38 phased release 2019-03-12 (IGSR/EBI FTP) |
| Gene annotation | Ensembl v115, `Homo_sapiens.GRCh38.115.gtf.gz` |

</details>

