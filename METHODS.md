# Methods and interpretation

## Preprocessing

The burden matrix is validated before analysis.

The maintained pipeline does not silently:

- drop duplicate sample IDs;
- coerce invalid text to zero;
- continue with non-finite values;
- continue after variance filtering leaves fewer than two features.

Zero imputation is available only through an explicit `--impute-zero` option and is recorded in the output metadata.

If `--log1p` is enabled, all burden values must be non-negative.

## Standardization

Each retained gene is standardized with `StandardScaler` before PCA.

This gives every retained gene unit variance in the PCA input and prevents large-scale features from dominating solely because of scale.

## PCA

A full PCA is first fit to estimate the explained-variance curve.

The number of retained PCs is the smallest number reaching `--pca-var`, subject to retaining at least two components and the rank available from the sample/feature matrix.

PCA is unsupervised. Sample labels are not used in fitting.

## Gene loadings

For the first five retained PCs, the pipeline reports genes with the largest absolute loadings.

A large loading means a gene contributes strongly to a variance direction in this cohort after preprocessing and standardization.

It does not imply:

- disease causality;
- statistical association with PSC;
- biological direction of effect;
- replication in another cohort.

## k-means selection

Silhouette score depends on both within-cluster compactness and separation from other clusters.

A single k-means run can be sensitive to initialization, so PSC-BurdenMap repeats each candidate `k` with multiple random seeds and reports both silhouette and ARI-based stability.

The selected `k` is an exploratory model-selection choice, not proof that the cohort contains a fixed biological number of subtypes.

## Cluster stability

For a fixed `k`, pairwise Adjusted Rand Index is computed between repeated k-means solutions.

High ARI means the optimization repeatedly finds similar partitions on the same PCA representation.

It does not measure cohort-resampling stability. A future extension could add bootstrap/subsampling stability to assess sensitivity to which individuals are included.

## Optional labels

Observed labels are merged after unsupervised fitting.

ARI and Normalized Mutual Information summarize cluster/label agreement without treating either label set as ordered.

These are descriptive concordance statistics, not supervised predictive performance.

## Important limitations

- no external cohort validation;
- no formal significance test for cluster count;
- no correction for ancestry, batch effects, sequencing depth, or other covariates unless they were addressed before constructing the burden matrix;
- no supervised PSC-vs-control classifier;
- no causal inference;
- PCA and k-means can be sensitive to preprocessing choices;
- rare-variant burden definitions strongly affect the geometry being analyzed.

Any publication-scale analysis should preserve the exact burden-generation procedure, inclusion/exclusion criteria, variant filtering rules, gene annotation version, and cohort metadata separately from this clustering code.
