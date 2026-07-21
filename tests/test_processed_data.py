from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
PROCESSED = REPO_ROOT / "data" / "processed"


def test_processed_data_invariants() -> None:
    blocks = pd.read_parquet(PROCESSED / "ma_block_points_proj_2023.parquet")
    clr = pd.read_parquet(PROCESSED / "ma_tract_clr_2023.parquet")
    strata = pd.read_parquet(PROCESSED / "ma_tract_strata_2023.parquet")

    assert blocks.shape == (119_181, 7)
    assert clr.shape == (1_598, 19)
    assert strata.shape == (1_598, 5)
    assert int(blocks["jobs"].sum()) == 3_203_251
    assert blocks["sector"].nunique() == 19
    assert blocks[["x_m", "y_m"]].apply(np.isfinite).all().all()

    density = blocks[["tract", "density_q"]].drop_duplicates()
    assert len(density) == 1_598
    assert density["density_q"].value_counts().to_dict() == {
        "Q1_low": 400,
        "Q2": 399,
        "Q3": 399,
        "Q4_high": 400,
    }

    assert clr.index.name == "tract"
    assert clr.index.is_unique
    assert np.abs(clr.sum(axis=1)).max() < 1e-8
    assert set(clr.index) == set(strata["tract"])
    assert set(clr.index) == set(density["tract"])
    assert strata["county_fips"].nunique() == 14
    assert not strata.isna().any().any()
