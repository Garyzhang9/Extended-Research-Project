#!/usr/bin/env python3
"""Power validation for the pairwise Euler-characteristic-curve test.

Reruns the retained PR8 simulation design against the EXACT analysis
statistic, imported directly from 05_topotest.py: the supremum difference
of normalised dimension-2 Vietoris-Rips ECC curves on a per-pair 200-point
grid over [0, q10 of the pair-pooled Aitchison distances]. The earlier PR8
record validated an L2 statistic on a fixed q90 grid and is superseded for
dissertation reporting.

Design: Gaussian generative model in the rank-18 zero-sum subspace
calibrated to the pooled CLR matrix; effect arms are a location shift along
the leading covariance eigenvector, covariance scale inflation, and a
two-component mixture separated along that eigenvector. Stage 1 is a
calibrated screen across all effect levels; Stage 2 is the exact
label-permutation test (999 permutations x 200 replicates) at the two
moderate shape effects (inflate c = 1.25, mixture s = 1.0).

Run a pipeline check first: python scripts/08_power_simulation.py --smoke
(its numbers are NOT a record). Full Stage 2 takes roughly 1.5-3 h on
8 cores because the q10 edge cutoff keeps the Rips complexes small.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

_SPEC = importlib.util.spec_from_file_location(
    "topotest05", Path(__file__).resolve().parent / "05_topotest.py"
)
_TOPOTEST = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_TOPOTEST)
shape_distance = _TOPOTEST.shape_distance
common_grid = _TOPOTEST.common_grid


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--jobs", type=int, default=int(os.environ.get("ERP_JOBS", "4"))
    )
    parser.add_argument("--seed", type=int, default=14196142)
    parser.add_argument("--n-cloud", type=int, default=400)
    parser.add_argument("--n-null", type=int, default=2000)
    parser.add_argument("--n-reps", type=int, default=500)
    parser.add_argument("--stage2-reps", type=int, default=200)
    parser.add_argument("--stage2-perms", type=int, default=999)
    parser.add_argument("--alpha", type=float, default=0.05)
    # Statistic parameters: keep identical to the 05_topotest.py defaults.
    parser.add_argument("--grid-steps", type=int, default=200)
    parser.add_argument("--distance-quantile", type=float, default=0.10)
    parser.add_argument("--max-dimension", type=int, default=2)
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def load_clr(repo_root: Path) -> np.ndarray:
    clr = pd.read_parquet(
        repo_root / "data" / "processed" / "ma_tract_clr_2023.parquet"
    )
    matrix = clr.select_dtypes("number").to_numpy(float)
    if matrix.shape != (1_598, 19):
        raise ValueError(f"Expected 1,598 x 19 CLR matrix, got {matrix.shape}")
    if np.abs(matrix.sum(axis=1)).max() >= 1e-6:
        raise ValueError("CLR rows do not sum to zero")
    return matrix


class Generator:
    """Gaussian model in the rank-18 zero-sum subspace, calibrated to pooled CLR."""

    def __init__(self, clr: np.ndarray, seed) -> None:
        self.rng = np.random.default_rng(seed)
        self.mu = clr.mean(axis=0)
        eigenvalues, eigenvectors = np.linalg.eigh(np.cov(clr, rowvar=False))
        self.eigenvalues = np.clip(eigenvalues, 0.0, None)
        self.eigenvectors = eigenvectors
        self.sigma1 = float(np.sqrt(self.eigenvalues[-1]))
        self.v1 = eigenvectors[:, -1]

    def _gauss(self, n: int, mu: np.ndarray, scale: float = 1.0) -> np.ndarray:
        z = self.rng.standard_normal((n, len(mu)))
        points = mu + (z * np.sqrt(scale * self.eigenvalues)) @ self.eigenvectors.T
        return points - points.mean(axis=1, keepdims=True)

    def null(self, n: int) -> np.ndarray:
        return self._gauss(n, self.mu)

    def location(self, delta: float, n: int) -> np.ndarray:
        return self._gauss(n, self.mu + delta * self.sigma1 * self.v1)

    def inflate(self, scale: float, n: int) -> np.ndarray:
        return self._gauss(n, self.mu, scale=scale)

    def mixture(self, separation: float, n: int) -> np.ndarray:
        half = n // 2
        low = self._gauss(half, self.mu - 0.5 * separation * self.sigma1 * self.v1)
        high = self._gauss(n - half, self.mu + 0.5 * separation * self.sigma1 * self.v1)
        return np.vstack([low, high])


def pair_statistic(
    x: np.ndarray, y: np.ndarray, quantile: float, steps: int, max_dimension: int
) -> float:
    grid = common_grid(x, y, quantile, steps)
    return shape_distance(x, y, grid, max_dimension)


def one_exact_replicate(
    child_seed: np.random.SeedSequence,
    clr: np.ndarray,
    arm: str,
    parameter: float,
    settings: dict,
) -> bool:
    generator_seed, permutation_seed = child_seed.spawn(2)
    generator = Generator(clr, generator_seed)
    rng = np.random.default_rng(permutation_seed)
    x = generator.null(settings["n_cloud"])
    y = getattr(generator, arm)(parameter, settings["n_cloud"])
    pooled = np.vstack([x, y])
    grid = common_grid(x, y, settings["quantile"], settings["steps"])
    observed = shape_distance(x, y, grid, settings["max_dimension"])
    exceed = 0
    for _ in range(settings["n_perm"]):
        order = rng.permutation(len(pooled))
        statistic = shape_distance(
            pooled[order[: settings["n_cloud"]]],
            pooled[order[settings["n_cloud"] :]],
            grid,
            settings["max_dimension"],
        )
        if statistic >= observed:
            exceed += 1
    p_value = (1 + exceed) / (settings["n_perm"] + 1)
    return p_value <= settings["alpha"]


def main() -> None:
    args = parse_args()
    if args.smoke:
        args.n_null, args.n_reps = 100, 50
        args.stage2_reps, args.stage2_perms = 5, 49
        print("SMOKE MODE: numbers below are NOT a record")
    if min(args.jobs, args.n_null, args.n_reps, args.stage2_reps,
           args.stage2_perms) < 1:
        raise ValueError("all replication counts and jobs must be positive")

    repo_root = Path(__file__).resolve().parents[1]
    table_dir = repo_root / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    clr = load_clr(repo_root)
    generator = Generator(clr, args.seed)
    settings = dict(
        n_cloud=args.n_cloud,
        quantile=args.distance_quantile,
        steps=args.grid_steps,
        max_dimension=args.max_dimension,
        n_perm=args.stage2_perms,
        alpha=args.alpha,
    )

    # Stage 1: calibrated screen. Clouds are drawn serially in the parent
    # process from one generator (reproducible); only the statistic is
    # evaluated in parallel.
    print(f"Stage 1: {args.n_null} null pairs for calibration ...")
    start = time.time()
    null_statistics = np.asarray(
        Parallel(n_jobs=args.jobs)(
            delayed(pair_statistic)(
                generator.null(args.n_cloud),
                generator.null(args.n_cloud),
                args.distance_quantile,
                args.grid_steps,
                args.max_dimension,
            )
            for _ in range(args.n_null)
        )
    )
    critical_value = float(np.quantile(null_statistics, 1 - args.alpha))
    print(
        f"  done in {time.time() - start:.0f}s; "
        f"sup-ECC critical value = {critical_value:.4f}"
    )

    levels = (
        [("location", delta) for delta in (0.5, 1.0, 2.0)]
        + [("inflate", scale) for scale in (1.10, 1.25, 1.50)]
        + [("mixture", separation) for separation in (0.5, 1.0, 1.5)]
    )
    screen_rows = []
    for arm, parameter in levels:
        statistics = np.asarray(
            Parallel(n_jobs=args.jobs)(
                delayed(pair_statistic)(
                    generator.null(args.n_cloud),
                    getattr(generator, arm)(parameter, args.n_cloud),
                    args.distance_quantile,
                    args.grid_steps,
                    args.max_dimension,
                )
                for _ in range(args.n_reps)
            )
        )
        power = float(np.mean(statistics > critical_value))
        screen_rows.append(
            dict(arm=arm, parameter=parameter, power_screen=power,
                 n_reps=args.n_reps)
        )
        print(f"  {arm:>8} {parameter:>5}: screen power {power:.3f}")
    screen = pd.DataFrame(screen_rows)
    screen.to_csv(table_dir / "power_simulation_screen.csv", index=False)

    # Stage 2: exact label-permutation power at the two moderate shape
    # effects, replicates parallelised with per-replicate spawned seeds.
    gate_levels = [("inflate", 1.25), ("mixture", 1.0)]
    child_seeds = np.random.SeedSequence(args.seed).spawn(
        args.stage2_reps * len(gate_levels)
    )
    exact_rows = []
    for index, (arm, parameter) in enumerate(gate_levels):
        start = time.time()
        seeds = child_seeds[
            index * args.stage2_reps : (index + 1) * args.stage2_reps
        ]
        rejections = Parallel(n_jobs=args.jobs)(
            delayed(one_exact_replicate)(seed, clr, arm, parameter, settings)
            for seed in seeds
        )
        power = float(np.mean(rejections))
        exact_rows.append(
            dict(arm=arm, parameter=parameter, power_exact=power,
                 n_reps=args.stage2_reps, n_perm=args.stage2_perms)
        )
        print(
            f"  {arm:>8} {parameter:>5}: EXACT power {power:.3f} "
            f"({time.time() - start:.0f}s)"
        )
    exact = pd.DataFrame(exact_rows)
    exact.to_csv(table_dir / "power_simulation_exact.csv", index=False)

    summary = dict(
        statistic=(
            "sup |normalised dim-2 Rips ECC difference| on a per-pair "
            f"{args.grid_steps}-point grid over [0, "
            f"q{args.distance_quantile:.2f} of pair-pooled distances], "
            "imported from 05_topotest.py"
        ),
        smoke_mode=args.smoke,
        seed=args.seed,
        n_cloud=args.n_cloud,
        alpha=args.alpha,
        stage1=dict(
            n_null=args.n_null, n_reps=args.n_reps,
            critical_value=critical_value,
        ),
        stage2=exact.to_dict("records"),
        note=(
            "Supersedes the PR8 L2/q90 power record for dissertation "
            "reporting; that run validated a different statistic."
        ),
    )
    (table_dir / "power_simulation_summary.json").write_text(
        json.dumps(summary, indent=2)
    )
    print(screen.to_string(index=False))
    print(exact.to_string(index=False))
    print(f"Outputs written to {table_dir}/")


if __name__ == "__main__":
    main()
