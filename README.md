# PSC-BurdenMap: PCA + k-means Clustering of Gene-Level Variant Burden

## How it works (5 lines)
1) Load a sample × gene variant-burden matrix (PSC + controls).
2) Standardize features (and optionally log-transform burdens).
3) Apply PCA to reduce dimensionality and denoise the burden space.
4) Run k-means clustering in PCA space and pick k using the silhouette score.
5) Export plots + tables to interpret clusters and gene drivers (PCA loadings).
 ---
## Project question
Do individuals (PSC patients vs controls) show natural structure in gene-level variant-burden space, and which genes drive the dominant variance axes?
---
## Input data
- `data/burden_matrix.csv`: rows are individuals, columns are genes, values are burden (0/1, counts, or weighted scores).
- `data/labels.csv` (optional): `sample_id,label` where label ∈ {PSC, Control}.
---
## Reto structure
```text
PSC-BurdenMap/
├─ README.md
├─ requirements.txt
├─ src/
│  └─ run_pca_kmeans.py
├─ data/
│  ├─ burden_matrix.csv
│  └─ labels.csv                
└─ artifacts/
```
---
## Data file formats
data/burden_matrix.csv

first column: sample_id

all other columns: genes (e.g., ENSG or symbols)

values: numeric burden (0/1, counts, or weighted)
---
## Run
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

python src/run_pca_kmeans.py \
  --burden data/burden_matrix.csv \
  --labels data/labels.csv \
  --outdir artifacts \
  --log1p \
  --pca_var 0.90 \
  --kmin 2 --kmax 8
---
