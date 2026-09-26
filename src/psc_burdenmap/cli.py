from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .core import (
    AnalysisConfig,
    choose_n_components,
    cluster_label_alignment,
    evaluate_kmeans_grid,
    load_burden_matrix,
    load_labels,
    prepare_burden,
    top_loadings,
)


def _save_json(payload: dict[str, object], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_commit() -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    value = result.stdout.strip()
    return value or None


def _save_scree_plot(explained_ratio: np.ndarray, outpath: Path) -> None:
    cumulative = np.cumsum(explained_ratio)
    plt.figure()
    plt.plot(np.arange(1, len(cumulative) + 1), cumulative)
    plt.xlabel("Number of Principal Components")
    plt.ylabel("Cumulative Explained Variance")
    plt.title("PCA Explained Variance")
    plt.ylim(0.0, 1.02)
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()


def _save_scatter(
    frame: pd.DataFrame,
    outpath: Path,
    *,
    color_col: str | None,
    title: str,
) -> None:
    plt.figure()

    if color_col is not None and color_col in frame.columns:
        values = frame[color_col].astype("string").fillna("NA")
        codes, names = pd.factorize(values, sort=True)
        scatter = plt.scatter(
            frame["PC1"],
            frame["PC2"],
            c=codes,
        )
        handles, _ = scatter.legend_elements()
        if len(handles) == len(names):
            plt.legend(
                handles,
                [str(name) for name in names],
                title=color_col,
                loc="best",
            )
    else:
        plt.scatter(frame["PC1"], frame["PC2"])

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(title)
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()


def _cluster_summary(frame: pd.DataFrame) -> pd.DataFrame:
    totals = (
        frame.groupby("cluster", dropna=False)["sample_id"]
        .count()
        .rename("n_total")
        .reset_index()
    )

    if "label" not in frame.columns or frame["label"].isna().all():
        return totals

    by_label = (
        frame.dropna(subset=["label"])
        .groupby(["cluster", "label"], dropna=False)["sample_id"]
        .count()
        .rename("n")
        .reset_index()
    )
    summary = by_label.merge(totals, on="cluster", how="left")
    summary["fraction_within_cluster"] = summary["n"] / summary["n_total"]
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "PSC-BurdenMap: validated PCA and stability-aware k-means "
            "for gene-level burden matrices."
        )
    )
    parser.add_argument(
        "--burden",
        type=Path,
        default=Path("data/burden_matrix.csv"),
    )
    parser.add_argument(
        "--labels",
        type=Path,
        default=Path("data/labels.csv"),
        help="Optional sample_id,label CSV. Missing file is treated as no labels.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("artifacts"),
    )
    parser.add_argument(
        "--log1p",
        action="store_true",
        help="Apply log1p; requires non-negative burden values.",
    )
    parser.add_argument(
        "--impute-zero",
        action="store_true",
        help=(
            "Explicitly replace missing/non-numeric burden entries with zero. "
            "Without this flag such values are rejected."
        ),
    )
    parser.add_argument(
        "--var-thresh",
        type=float,
        default=0.0,
        help="Drop genes with variance <= this threshold after optional log1p.",
    )
    parser.add_argument(
        "--pca-var",
        type=float,
        default=0.90,
        help="Keep enough PCs to reach this cumulative explained variance.",
    )
    parser.add_argument("--kmin", type=int, default=2)
    parser.add_argument("--kmax", type=int, default=8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--stability-runs",
        type=int,
        default=10,
        help="Repeat k-means with different seeds for each candidate k.",
    )
    parser.add_argument(
        "--n-init",
        type=int,
        default=20,
        help="KMeans n_init for each stability run.",
    )
    parser.add_argument(
        "--top-loadings",
        type=int,
        default=25,
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    config = AnalysisConfig(
        log1p=args.log1p,
        impute_zero=args.impute_zero,
        var_thresh=args.var_thresh,
        pca_var=args.pca_var,
        kmin=args.kmin,
        kmax=args.kmax,
        seed=args.seed,
        stability_runs=args.stability_runs,
        n_init=args.n_init,
    )

    args.outdir.mkdir(parents=True, exist_ok=True)

    burden = load_burden_matrix(
        args.burden,
        impute_zero=config.impute_zero,
    )
    labels = load_labels(args.labels)
    prepared = prepare_burden(burden, config)

    scaler = StandardScaler()
    standardized = scaler.fit_transform(
        prepared.features.to_numpy(dtype=float)
    )

    pca_full = PCA(svd_solver="full")
    full_scores = pca_full.fit_transform(standardized)
    max_components = min(
        standardized.shape[0],
        standardized.shape[1],
    )
    n_components = choose_n_components(
        pca_full.explained_variance_ratio_,
        config.pca_var,
        max_components=max_components,
    )

    pca = PCA(
        n_components=n_components,
        svd_solver="full",
    )
    pcs = pca.fit_transform(standardized)

    k_selection, best_k = evaluate_kmeans_grid(
        pcs,
        kmin=config.kmin,
        kmax=config.kmax,
        seed=config.seed,
        stability_runs=config.stability_runs,
        n_init=config.n_init,
    )
    final_kmeans = KMeans(
        n_clusters=best_k,
        random_state=config.seed,
        n_init=config.n_init,
    )
    clusters = final_kmeans.fit_predict(pcs)

    pc_columns = [
        f"PC{index + 1}"
        for index in range(pcs.shape[1])
    ]
    embedding = pd.DataFrame(
        pcs,
        columns=pc_columns,
    )
    embedding.insert(0, "sample_id", prepared.sample_ids)

    n_labels_extra = 0
    if not labels.empty:
        burden_ids = set(prepared.sample_ids)
        n_labels_extra = int(
            (~labels["sample_id"].isin(burden_ids)).sum()
        )
        embedding = embedding.merge(
            labels,
            on="sample_id",
            how="left",
            validate="one_to_one",
        )

    embedding["cluster"] = clusters.astype(int)

    embedding.to_csv(
        args.outdir / "pca_embedding.csv",
        index=False,
    )
    embedding[
        [
            column
            for column in embedding.columns
            if column == "sample_id"
            or column == "label"
            or column == "cluster"
        ]
    ].to_csv(
        args.outdir / "sample_clusters.csv",
        index=False,
    )
    k_selection.to_csv(
        args.outdir / "k_selection.csv",
        index=False,
    )
    _cluster_summary(embedding).to_csv(
        args.outdir / "cluster_summary.csv",
        index=False,
    )

    loading_frame = top_loadings(
        pca.components_,
        prepared.features.columns.astype(str).tolist(),
        top_n=args.top_loadings,
    )
    loading_frame.to_csv(
        args.outdir / "top_gene_loadings.csv",
        index=False,
    )

    cumulative = np.cumsum(
        pca_full.explained_variance_ratio_
    )
    variance_frame = pd.DataFrame(
        {
            "component": np.arange(
                1,
                len(cumulative) + 1,
            ),
            "explained_variance_ratio": (
                pca_full.explained_variance_ratio_
            ),
            "cumulative_explained_variance": cumulative,
        }
    )
    variance_frame.to_csv(
        args.outdir / "pca_explained_variance.csv",
        index=False,
    )

    _save_scree_plot(
        pca_full.explained_variance_ratio_,
        args.outdir / "pca_variance.png",
    )
    _save_scatter(
        embedding,
        args.outdir / "pca_scatter_by_cluster.png",
        color_col="cluster",
        title=f"PCA + k-means (k={best_k})",
    )
    if "label" in embedding.columns and embedding["label"].notna().any():
        _save_scatter(
            embedding,
            args.outdir / "pca_scatter_by_label.png",
            color_col="label",
            title="PCA of Gene Burden by Observed Label",
        )

    best_row = k_selection.loc[
        k_selection["k"] == best_k
    ].iloc[0]

    label_metrics: dict[str, float | int | None]
    if "label" in embedding.columns:
        label_metrics = cluster_label_alignment(
            clusters,
            embedding["label"].tolist(),
        )
    else:
        label_metrics = {
            "n_labeled_samples": 0,
            "n_label_classes": 0,
            "adjusted_rand_index": None,
            "normalized_mutual_information": None,
        }

    metrics: dict[str, object] = {
        "n_samples": int(len(prepared.sample_ids)),
        "n_genes_input": prepared.n_genes_input,
        "n_genes_after_var_filter": prepared.n_genes_after_filter,
        "n_genes_removed_by_variance": prepared.n_genes_removed_by_variance,
        "n_missing_values_imputed_to_zero": prepared.n_missing_imputed,
        "log1p": config.log1p,
        "var_thresh": config.var_thresh,
        "pca_cumvar_target": config.pca_var,
        "pca_n_components_kept": int(n_components),
        "pca_cumulative_variance_kept": float(
            pca.explained_variance_ratio_.sum()
        ),
        "requested_k_range": [config.kmin, config.kmax],
        "evaluated_k_range": [
            int(k_selection["k"].min()),
            int(k_selection["k"].max()),
        ],
        "best_k": best_k,
        "best_silhouette_mean": float(best_row["silhouette_mean"]),
        "best_silhouette_std": float(best_row["silhouette_std"]),
        "best_stability_ari_mean": float(best_row["stability_ari_mean"]),
        "best_stability_ari_min": float(best_row["stability_ari_min"]),
        "cluster_label_alignment": label_metrics,
    }
    _save_json(
        metrics,
        args.outdir / "run_metrics.json",
    )

    qc = {
        "burden_file": str(args.burden.resolve()),
        "labels_file": (
            str(args.labels.resolve())
            if args.labels.exists()
            else None
        ),
        "n_samples": int(len(prepared.sample_ids)),
        "n_labels_rows": int(len(labels)),
        "n_labels_not_in_burden_matrix": n_labels_extra,
        "n_burden_samples_without_label": int(
            embedding["label"].isna().sum()
            if "label" in embedding.columns
            else len(embedding)
        ),
        "n_missing_values_imputed_to_zero": prepared.n_missing_imputed,
        "impute_zero_enabled": config.impute_zero,
    }
    _save_json(
        qc,
        args.outdir / "data_qc.json",
    )

    metadata = {
        "command": "psc-burdenmap",
        "arguments": {
            key: (
                str(value)
                if isinstance(value, Path)
                else value
            )
            for key, value in vars(args).items()
        },
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scikit_learn": sklearn.__version__,
        "git_commit": _git_commit(),
        "input_sha256": {
            "burden": _sha256(args.burden),
            "labels": (
                _sha256(args.labels)
                if args.labels.exists()
                else None
            ),
        },
    }
    _save_json(
        metadata,
        args.outdir / "run_metadata.json",
    )

    print(f"Done. Artifacts in: {args.outdir.resolve()}")
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
