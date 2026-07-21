#!/usr/bin/env python3
"""Pairwise Euler-characteristic-curve tests for CLR point-cloud shape."""

from __future__ import annotations

import argparse
import os
from itertools import combinations
from pathlib import Path

import gudhi
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.spatial.distance import pdist, squareform
from statsmodels.stats.multitest import multipletests


ORDER = ["Q1_low", "Q2", "Q3", "Q4_high"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode", choices=("free", "county", "both"), default="both"
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=int(os.environ.get("ERP_PERMUTATIONS", "9999")),
    )
    parser.add_argument(
        "--jobs", type=int, default=int(os.environ.get("ERP_JOBS", "4"))
    )
    parser.add_argument("--seed", type=int, default=14196142)
    parser.add_argument("--grid-steps", type=int, default=200)
    parser.add_argument("--distance-quantile", type=float, default=0.10)
    parser.add_argument("--max-dimension", type=int, default=2)
    return parser.parse_args()


def load_analysis_frame(repo_root: Path) -> tuple[pd.DataFrame, list[str]]:
    processed = repo_root / "data" / "processed"
    clr = pd.read_parquet(processed / "ma_tract_clr_2023.parquet").reset_index()
    if "tract" not in clr.columns:
        clr = clr.rename(columns={clr.columns[0]: "tract"})
    density = (
        pd.read_parquet(processed / "ma_block_points_proj_2023.parquet")
        [["tract", "density_q"]]
        .drop_duplicates("tract")
    )
    county = pd.read_parquet(processed / "ma_tract_strata_2023.parquet")[
        ["tract", "county_fips"]
    ]
    dat = clr.merge(density, on="tract", validate="one_to_one").merge(
        county, on="tract", validate="one_to_one"
    )
    clr_columns = [column for column in dat if column.startswith("clr_CNS")]
    if len(dat) != 1_598 or len(clr_columns) != 19:
        raise ValueError("Expected 1,598 tracts and 19 CLR coordinates")
    if dat[["density_q", "county_fips"]].isna().any().any():
        raise ValueError("Density or county labels are missing")
    if np.abs(dat[clr_columns].sum(axis=1)).max() >= 1e-8:
        raise ValueError("CLR rows do not sum to zero")
    return dat, clr_columns


def euler_curve(points: np.ndarray, grid: np.ndarray, max_dimension: int) -> np.ndarray:
    distance_matrix = squareform(pdist(points))
    complex_ = gudhi.RipsComplex(
        distance_matrix=distance_matrix, max_edge_length=float(grid[-1])
    )
    tree = complex_.create_simplex_tree(max_dimension=max_dimension)
    simplices = [(len(simplex) - 1, filtration) for simplex, filtration in tree.get_filtration()]
    dimensions = np.fromiter((item[0] for item in simplices), dtype=int)
    filtrations = np.fromiter((item[1] for item in simplices), dtype=float)
    signs = (-1.0) ** dimensions
    return np.array([signs[filtrations <= radius].sum() for radius in grid])


def shape_distance(
    x: np.ndarray, y: np.ndarray, grid: np.ndarray, max_dimension: int
) -> float:
    curve_x = euler_curve(x, grid, max_dimension) / len(x)
    curve_y = euler_curve(y, grid, max_dimension) / len(y)
    return float(np.max(np.abs(curve_x - curve_y)))


def common_grid(
    x: np.ndarray, y: np.ndarray, quantile: float, steps: int
) -> np.ndarray:
    pooled = np.vstack([x, y])
    upper = float(np.quantile(pdist(pooled), quantile))
    return np.linspace(0.0, upper, steps)


def one_free_permutation(
    child_seed: np.random.SeedSequence,
    pooled: np.ndarray,
    first_size: int,
    grid: np.ndarray,
    max_dimension: int,
) -> float:
    rng = np.random.default_rng(child_seed)
    order = rng.permutation(len(pooled))
    return shape_distance(
        pooled[order[:first_size]], pooled[order[first_size:]], grid, max_dimension
    )


def county_mask(
    rng: np.random.Generator,
    county_codes: np.ndarray,
    first_counts: dict[object, int],
) -> np.ndarray:
    mask = np.zeros(len(county_codes), dtype=bool)
    for county in np.unique(county_codes):
        positions = np.flatnonzero(county_codes == county)
        selected = rng.permutation(positions)[: first_counts[county]]
        mask[selected] = True
    return mask


def one_county_permutation(
    child_seed: np.random.SeedSequence,
    pooled: np.ndarray,
    county_codes: np.ndarray,
    first_counts: dict[object, int],
    grid: np.ndarray,
    max_dimension: int,
) -> float:
    rng = np.random.default_rng(child_seed)
    mask = county_mask(rng, county_codes, first_counts)
    return shape_distance(pooled[mask], pooled[~mask], grid, max_dimension)


def run_pair(
    dat: pd.DataFrame,
    clr_columns: list[str],
    first: str,
    second: str,
    mode: str,
    permutations: int,
    jobs: int,
    seed: int,
    quantile: float,
    steps: int,
    max_dimension: int,
) -> tuple[dict[str, object], np.ndarray]:
    first_rows = dat[dat["density_q"] == first]
    second_rows = dat[dat["density_q"] == second]
    x = first_rows[clr_columns].to_numpy(float)
    y = second_rows[clr_columns].to_numpy(float)
    pooled = np.vstack([x, y])
    grid = common_grid(x, y, quantile, steps)
    observed = shape_distance(x, y, grid, max_dimension)
    child_seeds = np.random.SeedSequence(seed).spawn(permutations)

    if mode == "free":
        null = Parallel(n_jobs=jobs)(
            delayed(one_free_permutation)(
                child, pooled, len(x), grid, max_dimension
            )
            for child in child_seeds
        )
    else:
        counties = np.concatenate(
            [first_rows["county_fips"].to_numpy(), second_rows["county_fips"].to_numpy()]
        )
        is_first = np.concatenate(
            [np.ones(len(x), dtype=bool), np.zeros(len(y), dtype=bool)]
        )
        first_counts = {
            county: int(is_first[counties == county].sum())
            for county in np.unique(counties)
        }
        null = Parallel(n_jobs=jobs)(
            delayed(one_county_permutation)(
                child, pooled, counties, first_counts, grid, max_dimension
            )
            for child in child_seeds
        )

    null_array = np.asarray(null)
    p_value = (1 + np.sum(null_array >= observed)) / (permutations + 1)
    return (
        {
            "pair": f"{first} vs {second}",
            "n": len(x) + len(y),
            "D_obs": observed,
            "null_median": float(np.median(null_array)),
            "p_raw": p_value,
            "permutation": mode,
            "permutations": permutations,
            "seed": seed,
        },
        null_array,
    )


def run_mode(
    dat: pd.DataFrame,
    clr_columns: list[str],
    mode: str,
    args: argparse.Namespace,
    table_dir: Path,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for first, second in combinations(ORDER, 2):
        row, null = run_pair(
            dat,
            clr_columns,
            first,
            second,
            mode,
            args.permutations,
            args.jobs,
            args.seed,
            args.distance_quantile,
            args.grid_steps,
            args.max_dimension,
        )
        rows.append(row)
        print(
            f"{mode:6s} {row['pair']}: D={row['D_obs']:.4f}, "
            f"p={row['p_raw']:.4f}"
        )
        if (first, second) == ("Q1_low", "Q4_high"):
            pd.DataFrame({"D_null": null}).to_csv(
                table_dir / f"topotest_q1_q4_null_{mode}.csv", index=False
            )

    result = pd.DataFrame(rows)
    result["p_BH"] = multipletests(result["p_raw"], method="fdr_bh")[1]
    result["sig_BH"] = result["p_BH"] < 0.05
    result.to_csv(table_dir / f"rq2_topotest_{mode}.csv", index=False)
    return result


def main() -> None:
    args = parse_args()
    if args.permutations < 1 or args.jobs < 1:
        raise ValueError("permutations and jobs must be positive")
    repo_root = Path(__file__).resolve().parents[1]
    table_dir = repo_root / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    dat, clr_columns = load_analysis_frame(repo_root)
    modes = ("free", "county") if args.mode == "both" else (args.mode,)
    for mode in modes:
        result = run_mode(dat, clr_columns, mode, args, table_dir)
        print(result.to_string(index=False))


if __name__ == "__main__":
    main()
