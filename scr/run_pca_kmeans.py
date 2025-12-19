#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: seirana
"""

import argparse
import os
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def load_burden_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "sample_id" not in df.columns:
        raise ValueError("burden_matrix.csv must contain a 'sample_id' column.")
    df = df.drop_duplicates("sample_id")
    return df


def load_labels(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["sample_id", "label"])
    lab = pd.read_csv(path)
    if not {"sample_id", "label"}.issubset(lab.columns):
        raise ValueError("labels.csv must contain columns: sample_id,label")
    lab = lab.drop_duplicates("sample_id")
    return lab


def choose_n_components(explained_var_ratio: np.ndarray, target: float) -> int:
    cum = np.cumsum(explained_var_ratio)
    n = int(np.searchsorted(cum, target) + 1)
    return max(2, n)


def save_scree_plot(pca: PCA, outpath: Path) -> None:
    evr = pca.explained_variance_ratio_
    plt.figure()
    plt.plot(np.arange(1, len(evr) + 1), np.cumsum(evr))
    plt.xlabel("Number of Principal Components")
    plt.ylabel("Cumulative Explained Variance")
    plt.title("PCA Explained Variance (Cumulative)")
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()


def save_scatter(pc_df: pd.DataFrame, outpath: Path, color_col: str, title: str) -> None:
    plt.figure()
    if color_col in pc_df.columns:
        # Map categories to integers for plotting without choosing colors
        cats = pc_df[color_col].astype(str).fillna("NA")
        cat_codes, cat_names = pd.factorize(cats)
        plt.scatter(pc_df["PC1"], pc_df["PC2"], c=cat_codes)
        # legend-like text (simple, no colors specified)
        for i, name in enumerate(cat_names[:12]):
            plt.text(
                0.02,
                0.98 - i * 0.05,
                f"{i}: {name}",
                transform=plt.gca().transAxes,
                verticalalignment="top",
            )
        plt.ylabel("PC2")
        plt.xlabel("PC1")
        plt.title(title + f" (color={color_col})")
    else:
        plt.scatter(pc_df["PC1"], pc_df["PC2"])
        plt.ylabel("PC2")
        plt.xlabel("PC1")
        plt.title(title)
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()


def top_loadings(pca: PCA, feature_names: list[str], top_n: int = 25) -> pd.DataFrame:
    # For PC1 and PC2 (extend if you want)
    rows = []
    for pc_idx in range(min(5, pca.components_.shape[0])):  # top 5 PCs by default
        comp = pca.components_[pc_idx]
        order = np.argsort(np.abs(comp))[::-1][:top_n]
        for rank, j in enumerate(order, 1):
            rows.append(
                {
                    "PC": f"PC{pc_idx+1}",
                    "rank": rank,
                    "gene": feature_names[j],
                    "loading": float(comp[j]),
                    "abs_loading": float(abs(comp[j])),
                }
            )
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description="PSC-BurdenMap: PCA + k-means on gene-level burden")
    ap.add_argument("--burden", type=str, default="data/burden_matrix.csv")
    ap.add_argument("--labels", type=str, default="data/labels.csv")
    ap.add_argument("--outdir", type=str, default="artifacts")
    ap.add_argument("--log1p", action="store_true", help="Apply log1p transform to burden values")
    ap.add_argument("--var_thresh", type=float, default=0.0, help="Drop genes with variance <= threshold after preprocessing")
    ap.add_argument("--pca_var", type=float, default=0.90, help="Keep enough PCs to reach this cumulative explained variance")
    ap.add_argument("--kmin", type=int, default=2)
    ap.add_argument("--kmax", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    burden_df = load_burden_matrix(Path(args.burden))
    labels_df = load_labels(Path(args.labels))

    sample_ids = burden_df["sample_id"].astype(str).tolist()
    X = burden_df.drop(columns=["sample_id"]).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    feature_names = X.columns.astype(str).tolist()

    if args.log1p:
        X = np.log1p(X)

    # variance filter (optional)
    variances = X.var(axis=0)
    keep = variances > args.var_thresh
    X = X.loc[:, keep]
    feature_names = list(np.array(feature_names)[keep.values])

    # standardize
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X.values)

    # PCA fit (full), then choose n components by target explained variance
    pca_full = PCA(random_state=args.seed)
    pcs_full = pca_full.fit_transform(Xs)
    n_keep = choose_n_components(pca_full.explained_variance_ratio_, args.pca_var)

    # Refit PCA with chosen number for clean outputs
    pca = PCA(n_components=n_keep, random_state=args.seed)
    pcs = pca.fit_transform(Xs)

    # Save PCA artifacts
    save_scree_plot(pca_full, outdir / "pca_variance.png")

    pc_df = pd.DataFrame(pcs[:, :2], columns=["PC1", "PC2"])
    pc_df.insert(0, "sample_id", sample_ids)

    if not labels_df.empty:
        pc_df = pc_df.merge(labels_df, on="sample_id", how="left")

    pc_df.to_csv(outdir / "pca_embedding.csv", index=False)

    save_scatter(
        pc_df,
        outdir / "pca_scatter_by_label.png",
        color_col="label",
        title="PCA of Gene Burden",
    )

    # k-means sweep on PCA space (use all kept PCs)
    k_rows = []
    best = None

    for k in range(args.kmin, args.kmax + 1):
        km = KMeans(n_clusters=k, random_state=args.seed, n_init=20)
        cl = km.fit_predict(pcs)
        sil = silhouette_score(pcs, cl) if k > 1 else np.nan
        k_rows.append({"k": k, "silhouette": float(sil), "inertia": float(km.inertia_)})
        if best is None or sil > best["silhouette"]:
            best = {"k": k, "silhouette": sil, "model": km, "clusters": cl}

    ksel = pd.DataFrame(k_rows).sort_values("k")
    ksel.to_csv(outdir / "k_selection.csv", index=False)

    # Final clustering at best k
    cluster = best["clusters"]
    pc_df["cluster"] = cluster.astype(int)
    pc_df.to_csv(outdir / "sample_clusters.csv", index=False)

    save_scatter(
        pc_df,
        outdir / "pca_scatter_by_cluster.png",
        color_col="cluster",
        title=f"PCA + k-means (best k={best['k']}, silhouette={best['silhouette']:.3f})",
    )

    # Cluster summary (with PSC/Control proportions if labels exist)
    if "label" in pc_df.columns:
        summ = (
            pc_df.groupby(["cluster", "label"])["sample_id"]
            .count()
            .rename("n")
            .reset_index()
        )
        summ_tot = pc_df.groupby("cluster")["sample_id"].count().rename("n_total").reset_index()
        summ = summ.merge(summ_tot, on="cluster", how="left")
        summ["fraction"] = summ["n"] / summ["n_total"]
        summ.to_csv(outdir / "cluster_summary.csv", index=False)
    else:
        summ_tot = pc_df.groupby("cluster")["sample_id"].count().rename("n_total").reset_index()
        summ_tot.to_csv(outdir / "cluster_summary.csv", index=False)

    # Top gene loadings
    load_df = top_loadings(pca, feature_names, top_n=25)
    load_df.to_csv(outdir / "top_gene_loadings.csv", index=False)

    # Metrics JSON for the run
    metrics = {
        "n_samples": int(len(sample_ids)),
        "n_genes_input": int(burden_df.shape[1] - 1),
        "n_genes_after_var_filter": int(len(feature_names)),
        "log1p": bool(args.log1p),
        "var_thresh": float(args.var_thresh),
        "pca_cumvar_target": float(args.pca_var),
        "pca_n_components_kept": int(n_keep),
        "best_k": int(best["k"]),
        "best_silhouette": float(best["silhouette"]),
    }
    with open(outdir / "run_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)

    print("Done. Artifacts in:", str(outdir.resolve()))
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
