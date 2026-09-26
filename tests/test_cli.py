import json
import sys

import pandas as pd

from psc_burdenmap.cli import main


def test_cli_writes_reproducible_artifacts(tmp_path, monkeypatch):
    burden_path = tmp_path / "burden.csv"
    labels_path = tmp_path / "labels.csv"
    output_dir = tmp_path / "artifacts"

    pd.DataFrame(
        {
            "sample_id": [
                "S1",
                "S2",
                "S3",
                "S4",
                "S5",
                "S6",
                "S7",
                "S8",
            ],
            "G1": [0, 0, 1, 1, 8, 9, 8, 9],
            "G2": [0, 1, 0, 1, 9, 8, 9, 8],
            "G3": [1, 1, 0, 0, 7, 7, 8, 8],
            "G4": [0, 1, 1, 0, 8, 8, 7, 7],
        }
    ).to_csv(burden_path, index=False)

    pd.DataFrame(
        {
            "sample_id": [
                "S1",
                "S2",
                "S3",
                "S4",
                "S5",
                "S6",
                "S7",
                "S8",
            ],
            "label": [
                "Control",
                "Control",
                "Control",
                "Control",
                "PSC",
                "PSC",
                "PSC",
                "PSC",
            ],
        }
    ).to_csv(labels_path, index=False)

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "psc-burdenmap",
            "--burden",
            str(burden_path),
            "--labels",
            str(labels_path),
            "--outdir",
            str(output_dir),
            "--kmin",
            "2",
            "--kmax",
            "3",
            "--stability-runs",
            "3",
            "--n-init",
            "5",
        ],
    )

    assert main() == 0

    expected = {
        "pca_embedding.csv",
        "sample_clusters.csv",
        "k_selection.csv",
        "cluster_summary.csv",
        "top_gene_loadings.csv",
        "pca_explained_variance.csv",
        "pca_variance.png",
        "pca_scatter_by_cluster.png",
        "pca_scatter_by_label.png",
        "run_metrics.json",
        "data_qc.json",
        "run_metadata.json",
    }
    assert expected.issubset(
        {path.name for path in output_dir.iterdir()}
    )

    metrics = json.loads(
        (output_dir / "run_metrics.json").read_text(
            encoding="utf-8"
        )
    )
    metadata = json.loads(
        (output_dir / "run_metadata.json").read_text(
            encoding="utf-8"
        )
    )

    assert metrics["n_samples"] == 8
    assert metrics["best_k"] in {2, 3}
    assert metadata["input_sha256"]["burden"]
    assert metadata["input_sha256"]["labels"]
