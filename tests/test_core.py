import numpy as np
import pandas as pd
import pytest

from psc_burdenmap.core import (
    AnalysisConfig,
    choose_n_components,
    cluster_label_alignment,
    evaluate_kmeans_grid,
    load_burden_matrix,
    prepare_burden,
    valid_k_values,
)


def test_duplicate_sample_ids_are_rejected(tmp_path):
    path = tmp_path / "burden.csv"
    pd.DataFrame(
        {
            "sample_id": ["S1", "S1", "S2"],
            "G1": [0, 1, 2],
            "G2": [1, 0, 1],
        }
    ).to_csv(path, index=False)

    with pytest.raises(ValueError, match="unique"):
        load_burden_matrix(path)


def test_missing_or_nonnumeric_values_require_explicit_imputation(tmp_path):
    path = tmp_path / "burden.csv"
    pd.DataFrame(
        {
            "sample_id": ["S1", "S2", "S3"],
            "G1": [0, "bad", 2],
            "G2": [1, 0, 1],
        }
    ).to_csv(path, index=False)

    with pytest.raises(ValueError, match="impute-zero"):
        load_burden_matrix(path)

    loaded = load_burden_matrix(path, impute_zero=True)
    assert loaded.loc[1, "G1"] == 0.0
    assert loaded.attrs["n_missing_imputed"] == 1


def test_log1p_rejects_negative_values():
    frame = pd.DataFrame(
        {
            "sample_id": ["S1", "S2", "S3"],
            "G1": [0.0, -1.0, 2.0],
            "G2": [1.0, 2.0, 3.0],
        }
    )

    with pytest.raises(ValueError, match="non-negative"):
        prepare_burden(
            frame,
            AnalysisConfig(log1p=True),
        )


def test_variance_filter_cannot_remove_almost_every_gene():
    frame = pd.DataFrame(
        {
            "sample_id": ["S1", "S2", "S3"],
            "G1": [1.0, 1.0, 1.0],
            "G2": [2.0, 2.0, 2.0],
        }
    )

    with pytest.raises(ValueError, match="fewer than two genes"):
        prepare_burden(
            frame,
            AnalysisConfig(var_thresh=0.0),
        )


def test_choose_n_components_is_bounded_and_keeps_two():
    assert choose_n_components(
        [0.8, 0.15, 0.05],
        0.75,
    ) == 2

    assert choose_n_components(
        [0.4, 0.3, 0.2, 0.1],
        0.95,
        max_components=3,
    ) == 3


def test_k_values_clip_to_silhouette_feasible_range():
    assert valid_k_values(
        n_samples=5,
        kmin=2,
        kmax=10,
    ) == [2, 3, 4]


def test_kmeans_grid_is_reproducible():
    rng = np.random.default_rng(5)
    first_cluster = rng.normal(
        loc=-3.0,
        scale=0.2,
        size=(12, 2),
    )
    second_cluster = rng.normal(
        loc=3.0,
        scale=0.2,
        size=(12, 2),
    )
    pcs = np.vstack([first_cluster, second_cluster])

    first, best_first = evaluate_kmeans_grid(
        pcs,
        kmin=2,
        kmax=4,
        seed=42,
        stability_runs=4,
        n_init=10,
    )
    second, best_second = evaluate_kmeans_grid(
        pcs,
        kmin=2,
        kmax=4,
        seed=42,
        stability_runs=4,
        n_init=10,
    )

    pd.testing.assert_frame_equal(first, second)
    assert best_first == best_second == 2


def test_cluster_label_alignment_reports_descriptive_agreement():
    metrics = cluster_label_alignment(
        [0, 0, 1, 1],
        ["PSC", "PSC", "Control", "Control"],
    )

    assert metrics["n_labeled_samples"] == 4
    assert metrics["adjusted_rand_index"] == pytest.approx(1.0)
    assert metrics["normalized_mutual_information"] == pytest.approx(1.0)
