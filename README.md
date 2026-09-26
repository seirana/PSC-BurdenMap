# PSC-BurdenMap

PSC-BurdenMap is a reproducible unsupervised analysis pipeline for exploring sample-level structure in gene-level variant-burden matrices.

It combines:

- strict burden-matrix validation;
- optional log1p transformation;
- variance filtering;
- feature standardization;
- PCA for dimensionality reduction;
- stability-aware k-means model selection;
- optional descriptive comparison with PSC/control labels;
- provenance, QC, tests, CI, and Docker.

The analysis is exploratory. PCA loadings identify genes that contribute to variance in the analyzed burden matrix; they do not establish causality. k-means clusters are mathematical partitions of the reduced feature space and are not automatically clinical subtypes.

## Analysis question

Given a sample × gene burden matrix:

1. is there reproducible low-dimensional structure among samples;
2. how stable is k-means clustering across random initializations;
3. which genes contribute most strongly to the leading PCs;
4. if labels are available, how closely do the unsupervised clusters align with those observed labels?

Labels are never used to fit PCA or k-means.

## Repository structure

```text
PSC-BurdenMap/
├── src/
│   ├── psc_burdenmap/
│   │   ├── core.py
│   │   └── cli.py
│   └── run_pca_kmeans.py
├── tests/
├── data/
├── .github/workflows/ci.yml
├── Dockerfile
├── METHODS.md
├── pyproject.toml
├── requirements.txt
└── README.md
```

The historical `src/run_pca_kmeans.py` entry point is retained as a thin compatibility wrapper.

## Installation

```bash
git clone https://github.com/seirana/PSC-BurdenMap.git
cd PSC-BurdenMap

python -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -e .
```

For development:

```bash
python -m pip install -e ".[dev]"
```

## Input data

### Burden matrix

Default:

```text
data/burden_matrix.csv
```

Format:

```text
sample_id,GeneA,GeneB,GeneC,...
S001,0,1,0,...
S002,2,0,1,...
```

Requirements:

- `sample_id` must be present, non-empty, and unique;
- at least two gene columns are required;
- gene values must be numeric and finite;
- missing or non-numeric values are rejected by default.

If zero imputation is scientifically justified, enable it explicitly with `--impute-zero`. The number of imputed values is recorded in QC and metrics outputs.

### Optional labels

Default:

```text
data/labels.csv
```

Format:

```text
sample_id,label
S001,PSC
S002,Control
```

A missing labels file is treated as an unlabeled analysis. When labels are present, they are merged only after PCA/k-means fitting.

See [data/README.md](data/README.md).

## Run

```bash
psc-burdenmap \
  --burden data/burden_matrix.csv \
  --labels data/labels.csv \
  --outdir artifacts \
  --log1p \
  --pca-var 0.90 \
  --kmin 2 \
  --kmax 8 \
  --stability-runs 10 \
  --seed 42
```

The historical command remains valid:

```bash
python src/run_pca_kmeans.py ...
```

The older underscore forms `--pca_var` and `--var_thresh` are also accepted.

## Model-selection behavior

For each feasible value of `k`, the pipeline repeats k-means with several seeds.

It reports:

- mean and standard deviation of silhouette score;
- mean and standard deviation of inertia;
- mean and minimum pairwise Adjusted Rand Index (ARI) across repeated clusterings.

The selected `k` maximizes mean silhouette, then stability ARI, with smaller `k` used as the final tie-breaker.

If the requested `kmax` exceeds the silhouette-feasible maximum `n_samples - 1`, the evaluated range is clipped and recorded in `run_metrics.json`.

## Outputs

Typical outputs:

```text
artifacts/
├── pca_embedding.csv
├── pca_explained_variance.csv
├── pca_variance.png
├── pca_scatter_by_cluster.png
├── pca_scatter_by_label.png        # only when labels are available
├── k_selection.csv
├── sample_clusters.csv
├── cluster_summary.csv
├── top_gene_loadings.csv
├── data_qc.json
├── run_metrics.json
└── run_metadata.json
```

`run_metadata.json` records the command arguments, Python/platform information, core package versions, Git commit when available, and SHA-256 checksums for the input files.

## Label comparison

When labels are available, the pipeline reports descriptive:

- Adjusted Rand Index;
- Normalized Mutual Information.

These metrics describe agreement between an unsupervised clustering and observed labels. They are not classification accuracy, because labels were not used to train the model.

## Testing

```bash
python -m pytest
python -m ruff check src tests
```

Tests cover:

- duplicate sample rejection;
- strict handling of missing/non-numeric burden values;
- explicit zero-imputation behavior;
- log1p input validation;
- variance-filter edge cases;
- PCA component selection;
- feasible k-range handling;
- deterministic k-selection;
- cluster/label agreement metrics;
- end-to-end artifact creation.

GitHub Actions runs the maintained package on Python 3.10, 3.11, and 3.12 and builds the Docker image.

## Docker

Build:

```bash
docker build -t psc-burdenmap .
```

Run:

```bash
mkdir -p artifacts

docker run --rm \
  -v "$PWD/data:/app/data:ro" \
  -v "$PWD/artifacts:/app/artifacts" \
  psc-burdenmap \
  --burden /app/data/burden_matrix.csv \
  --labels /app/data/labels.csv \
  --outdir /app/artifacts
```

## Interpretation

Appropriate interpretation:

> The analyzed cohort contains an unsupervised cluster structure with the reported silhouette/stability values, and the listed genes have large PCA loadings in this dataset.

Not appropriate:

> The clusters are validated PSC subtypes, or high-loading genes are causal PSC genes.

See [METHODS.md](METHODS.md) for statistical and interpretation details.

## License

No explicit license file is currently included. Repository visibility alone does not grant reuse rights.
