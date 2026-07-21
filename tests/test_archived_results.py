from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"


def test_mfunction_results() -> None:
    free = pd.read_csv(TABLES / "rq1_results.csv")
    county = pd.read_csv(TABLES / "rq1_county_restricted_results.csv")
    assert len(free) == 18
    assert set(free.loc[free["sig_BH"], "sector"]) == {
        "CNS01", "CNS05", "CNS06", "CNS07", "CNS12", "CNS19"
    }
    assert set(county.loc[county["sig_restr"], "sector"]) == {"CNS05", "CNS12"}
    assert int((county["flipped"] == "Y").sum()) == 4


def test_compositional_results() -> None:
    omnibus = pd.read_csv(TABLES / "rq2_omnibus_results.csv")
    free = omnibus[omnibus["permutation"] == "free"].set_index("test")
    assert np.isclose(free.loc["PERMANOVA", "statistic"], 23.268, atol=0.001)
    assert np.isclose(free.loc["PERMDISP", "statistic"], 55.764, atol=0.001)

    topological = pd.read_csv(TABLES / "rq2_topotest_free.csv")
    assert len(topological) == 6
    assert int(topological["sig_BH"].sum()) == 4
    extreme = topological[topological["pair"] == "Q1_low vs Q4_high"].iloc[0]
    assert np.isclose(extreme["D_obs"], 669.8175, atol=1e-4)


def test_ball_mapper_and_reconciliation_results() -> None:
    nodes = pd.read_csv(TABLES / "fig_5_4_nodes.csv")
    edges = pd.read_csv(TABLES / "fig_5_4_edges.csv")
    regional = pd.read_csv(TABLES / "regional_agreement.csv").iloc[0]
    sector = pd.read_csv(TABLES / "rq3_sector_spearman.csv").iloc[0]
    assert len(nodes) == 158
    assert len(edges) == 3_205
    assert int(nodes["is_isolated"].sum()) == 9
    assert np.isclose(regional["ARI"], 0.0274, atol=1e-4)
    assert np.isclose(regional["NMI"], 0.0366, atol=1e-4)
    assert np.isclose(sector["rho"], 0.1499486, atol=1e-7)
    assert np.isclose(sector["between_stratum_ss"], 2714.882, atol=0.001)
