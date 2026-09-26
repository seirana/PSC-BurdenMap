from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)


@dataclass(frozen=True)
class AnalysisConfig:
    log1p: bool = False
    impute_zero: bool = False
    var_thresh: float = 0.0
    pca_var: float = 0.90
    kmin: int = 2
    kmax: int = 8
    seed: int = 42
    stability_runs: int = 10
    n_init: int = 20

    def __post_init__(self) -> None:
        if self.var_thresh < 0:
            raise ValueError("var_thresh must be non-negative")
        if not 0.0 < self.pca_var <= 1.0:
            raise ValueError("pca_var must be in (0, 1]")
        if self.kmin < 2:
            raise ValueError("kmin must be at least 2")
        if self.kmax < self.kmin:
            raise ValueError("kmax must be greater than or equal to kmin")
        if self.stability_runs <= 0:
            raise ValueError("stability_runs must be greater than 0")
        if self.n_init <= 0:
            raise ValueError("n_init must be greater than 0")


@dataclass(frozen=True)
class PreparedBurden:
    sample_ids: list[str]
    features: pd.DataFrame
    n_genes_input: int
    n_genes_after_filter: int
    n_missing_imputed: int
    n_genes_removed_by_variance: int


def _validate_sample_ids(series: pd.Series) -> list[str]:
    if series.isna().any():
        raise ValueError("sample_id values must not be missing")

    sample_ids = series.astype(str).str.strip()
    if (sample_ids == "").any():
        raise ValueError("sample_id values must not be empty")
    if sample_ids.duplicated().any():
        duplicates = sorted(sample_ids[sample_ids.duplicated()].unique())
        raise ValueError(
            "sample_id values must be unique; duplicates include: "
            + ", ".join(duplicates[:10])
        )
    return sample_ids.tolist()


def load_burden_matrix(
    path: str | Path,
    *,
    impute_zero: bool = False,
) -> pd.DataFrame:
    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(f"Burden matrix not found: {source}")

    frame = pd.read_csv(source)
    if "sample_id" not in frame.columns:
        raise ValueError("Burden matrix must contain a 'sample_id' column.")
    _validate_sample_ids(frame["sample_id"])

    gene_columns = [column for column in frame.columns if column != "sample_id"]
    if len(gene_columns) < 2:
        raise ValueError("Burden matrix must contain at least two gene columns.")

    numeric = frame[gene_columns].apply(pd.to_numeric, errors="coerce")
    missing_count = int(numeric.isna().sum().sum())
    if missing_count and not impute_zero:
        columns = numeric.columns[numeric.isna().any()].astype(str).tolist()
        raise ValueError(
            "Burden matrix contains missing or non-numeric values. "
            "Use --impute-zero only when zero imputation is scientifically justified. "
            f"Affected columns include: {', '.join(columns[:10])}"
        )
    if missing_count:
        numeric = numeric.fillna(0.0)

    values = numeric.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Burden matrix contains non-finite values.")

    output = numeric.copy()
    output.insert(0, "sample_id", frame["sample_id"].astype(str).str.strip())
    output.attrs["n_missing_imputed"] = missing_count
    return output


def load_labels(path: str | Path | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame(columns=["sample_id", "label"])

    source = Path(path)
    if not source.exists():
        return pd.DataFrame(columns=["sample_id", "label"])

    labels = pd.read_csv(source)
    required = {"sample_id", "label"}
    if not required.issubset(labels.columns):
        raise ValueError("Labels file must contain columns: sample_id,label")

    labels = labels[["sample_id", "label"]].copy()
    labels["sample_id"] = labels["sample_id"].astype(str).str.strip()
    if labels["sample_id"].duplicated().any():
        raise ValueError("Labels file contains duplicate sample_id values.")
    if labels["sample_id"].eq("").any():
        raise ValueError("Labels file contains empty sample_id values.")

    labels["label"] = labels["label"].astype("string")
    return labels


def prepare_burden(
    burden_df: pd.DataFrame,
    config: AnalysisConfig,
) -> PreparedBurden:
    sample_ids = _validate_sample_ids(burden_df["sample_id"])
    features = burden_df.drop(columns=["sample_id"]).copy()

    if len(sample_ids) < 3:
        raise ValueError(
            "At least three samples are required for PCA plus silhouette-based clustering."
        )

    if config.log1p:
        if (features.to_numpy(dtype=float) < 0).any():
            raise ValueError("log1p requires non-negative burden values.")
        features = np.log1p(features)

    variances = features.var(axis=0, ddof=0)
    keep = variances > config.var_thresh
    filtered = features.loc[:, keep]

    if filtered.shape[1] < 2:
        raise ValueError(
            "Variance filtering left fewer than two genes; lower --var-thresh."
        )

    max_components = min(filtered.shape[0], filtered.shape[1])
    if max_components < 2:
        raise ValueError("At least two PCA components must be estimable.")

    return PreparedBurden(
        sample_ids=sample_ids,
        features=filtered,
        n_genes_input=int(features.shape[1]),
        n_genes_after_filter=int(filtered.shape[1]),
        n_missing_imputed=int(burden_df.attrs.get("n_missing_imputed", 0)),
        n_genes_removed_by_variance=int((~keep).sum()),
    )


def choose_n_components(
    explained_var_ratio: Sequence[float],
    target: float,
    *,
    max_components: int | None = None,
) -> int:
    if not 0.0 < target <= 1.0:
        raise ValueError("target must be in (0, 1]")

    ratios = np.asarray(explained_var_ratio, dtype=float)
    if ratios.ndim != 1 or ratios.size == 0:
        raise ValueError("explained_var_ratio must be a non-empty 1D sequence")
    if not np.isfinite(ratios).all() or (ratios < 0).any():
        raise ValueError("explained_var_ratio must contain finite non-negative values")

    available = int(ratios.size)
    if max_components is not None:
        available = min(available, int(max_components))
    if available < 2:
        raise ValueError("At least two PCA components are required.")

    cumulative = np.cumsum(ratios[:available])
    selected = int(np.searchsorted(cumulative, target, side="left") + 1)
    return min(max(2, selected), available)


def valid_k_values(
    n_samples: int,
    kmin: int,
    kmax: int,
) -> list[int]:
    if n_samples < 3:
        raise ValueError("At least three samples are required for silhouette scoring.")
    if kmin < 2:
        raise ValueError("kmin must be at least 2")
    if kmax < kmin:
        raise ValueError("kmax must be greater than or equal to kmin")

    upper = min(kmax, n_samples - 1)
    values = list(range(kmin, upper + 1))
    if not values:
        raise ValueError(
            "No valid k values remain; silhouette requires 2 <= k <= n_samples - 1."
        )
    return values


def _pairwise_ari(clusterings: list[np.ndarray]) -> list[float]:
    if len(clusterings) <= 1:
        return [1.0]

    return [
        float(adjusted_rand_score(left, right))
        for left, right in combinations(clusterings, 2)
    ]


def evaluate_kmeans_grid(
    pcs: np.ndarray,
    *,
    kmin: int,
    kmax: int,
    seed: int,
    stability_runs: int = 10,
    n_init: int = 20,
) -> tuple[pd.DataFrame, int]:
    array = np.asarray(pcs, dtype=float)
    if array.ndim != 2:
        raise ValueError("pcs must be a 2D array")
    if not np.isfinite(array).all():
        raise ValueError("pcs must contain only finite values")
    if stability_runs <= 0:
        raise ValueError("stability_runs must be greater than 0")
    if n_init <= 0:
        raise ValueError("n_init must be greater than 0")

    rows: list[dict[str, float | int]] = []
    for k in valid_k_values(array.shape[0], kmin, kmax):
        silhouettes: list[float] = []
        inertias: list[float] = []
        clusterings: list[np.ndarray] = []

        for run_index in range(stability_runs):
            model = KMeans(
                n_clusters=k,
                random_state=int(seed + run_index),
                n_init=n_init,
            )
            clusters = model.fit_predict(array)
            silhouettes.append(float(silhouette_score(array, clusters)))
            inertias.append(float(model.inertia_))
            clusterings.append(clusters)

        stability = _pairwise_ari(clusterings)
        rows.append(
            {
                "k": int(k),
                "silhouette_mean": float(np.mean(silhouettes)),
                "silhouette_std": float(np.std(silhouettes, ddof=0)),
                "inertia_mean": float(np.mean(inertias)),
                "inertia_std": float(np.std(inertias, ddof=0)),
                "stability_ari_mean": float(np.mean(stability)),
                "stability_ari_min": float(np.min(stability)),
                "n_stability_runs": int(stability_runs),
            }
        )

    results = pd.DataFrame(rows).sort_values("k").reset_index(drop=True)
    ranked = results.sort_values(
        ["silhouette_mean", "stability_ari_mean", "k"],
        ascending=[False, False, True],
        kind="stable",
    )
    best_k = int(ranked.iloc[0]["k"])
    return results, best_k


def top_loadings(
    components: np.ndarray,
    feature_names: Sequence[str],
    *,
    top_n: int = 25,
    max_pcs: int = 5,
) -> pd.DataFrame:
    if top_n <= 0:
        raise ValueError("top_n must be greater than 0")
    if max_pcs <= 0:
        raise ValueError("max_pcs must be greater than 0")

    matrix = np.asarray(components, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("components must be a 2D array")
    if matrix.shape[1] != len(feature_names):
        raise ValueError("feature_names length must match component width")

    rows: list[dict[str, float | int | str]] = []
    for pc_index in range(min(max_pcs, matrix.shape[0])):
        component = matrix[pc_index]
        order = np.argsort(np.abs(component))[::-1][:top_n]
        for rank, feature_index in enumerate(order, start=1):
            loading = float(component[feature_index])
            rows.append(
                {
                    "PC": f"PC{pc_index + 1}",
                    "rank": int(rank),
                    "gene": str(feature_names[feature_index]),
                    "loading": loading,
                    "abs_loading": abs(loading),
                }
            )

    return pd.DataFrame(rows)


def cluster_label_alignment(
    clusters: Sequence[int],
    labels: Sequence[object],
) -> dict[str, float | int | None]:
    cluster_array = np.asarray(clusters)
    label_series = pd.Series(labels, dtype="string")
    complete = label_series.notna()

    n_labeled = int(complete.sum())
    n_label_classes = int(label_series[complete].nunique())
    if n_labeled < 2 or n_label_classes < 2:
        return {
            "n_labeled_samples": n_labeled,
            "n_label_classes": n_label_classes,
            "adjusted_rand_index": None,
            "normalized_mutual_information": None,
        }

    labels_complete = label_series[complete].astype(str).to_numpy()
    clusters_complete = cluster_array[complete.to_numpy()]
    return {
        "n_labeled_samples": n_labeled,
        "n_label_classes": n_label_classes,
        "adjusted_rand_index": float(
            adjusted_rand_score(labels_complete, clusters_complete)
        ),
        "normalized_mutual_information": float(
            normalized_mutual_info_score(labels_complete, clusters_complete)
        ),
    }
