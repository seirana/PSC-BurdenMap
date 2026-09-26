# Data inputs

PSC-BurdenMap expects local cohort-derived inputs. Real individual-level genetic data should not be committed to a public repository.

## burden_matrix.csv

Required columns:

- `sample_id`;
- at least two gene-burden columns.

Each row represents one individual. Each gene column contains a numeric burden value such as a count, binary indicator, or weighted score.

The pipeline does not assume that different burden definitions are equivalent. Record how burden values were constructed.

## labels.csv

Optional columns:

- `sample_id`;
- `label`.

Labels are used only after unsupervised PCA/k-means fitting for plotting and descriptive agreement statistics.

## Reproducibility metadata to preserve

For a real study, record:

- cohort inclusion/exclusion criteria;
- sequencing platform and calling pipeline;
- genome build;
- variant QC thresholds;
- annotation tool/database versions;
- allele-frequency source and cutoff;
- burden aggregation rule;
- gene identifier version;
- ancestry/covariate handling;
- date/version of the exported burden matrix.

## Privacy

Individual-level genomic and phenotype data can be sensitive. Keep controlled-access or identifiable data outside the public repository and mount them locally when running the pipeline.
