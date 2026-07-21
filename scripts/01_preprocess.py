#!/usr/bin/env python3
"""Build the ERP analysis datasets from public source files.

Pipeline
--------
1. Aggregate 2023 Massachusetts LODES WAC JT02 jobs from blocks to tracts.
2. Join TIGER/Line land area, filter the analysis sample, and assign density
   quartiles.
3. Apply the Jeffreys 0.5 pseudocount and centred log-ratio transform.
4. Attach county metadata used by county-restricted permutations.
5. Join LODES block coordinates and project them to EPSG:26986 for the
   distance-based M-function pipeline.

The three load-bearing outputs reproduce the files used by the paper:

* ma_block_points_proj_2023.parquet
* ma_tract_clr_2023.parquet
* ma_tract_strata_2023.parquet
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import tempfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow
import pyproj


YEAR = 2023
MIN_PRIVATE_JOBS = 50
STATE_FIPS = "25"
TARGET_CRS = "EPSG:26986"
NAICS_COLUMNS = [f"CNS{i:02d}" for i in range(1, 20)]
DENSITY_LEVELS = ["Q1_low", "Q2", "Q3", "Q4_high"]

EXPECTED_TRACTS = 1_598
EXPECTED_PRIVATE_JOBS = 3_203_251
EXPECTED_BLOCK_SECTOR_ROWS = 119_181
EXPECTED_DENSITY_COUNTS = {
    "Q1_low": 400,
    "Q2": 399,
    "Q3": 399,
    "Q4_high": 400,
}

OUTPUT_NAMES = (
    "ma_tract_jobs_2023.parquet",
    "ma_tract_shares_2023.parquet",
    "ma_tract_shares_smoothed_2023.parquet",
    "ma_tract_clr_2023.parquet",
    "ma_tract_strata_2023.parquet",
    "ma_block_points_proj_2023.parquet",
    "preprocessing_summary.json",
    "checksums.sha256",
)


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=repo_root / "data" / "raw",
        help="Directory containing the four public source files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo_root / "data" / "processed",
        help="Destination for deterministic parquet outputs.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace existing generated outputs after all checks pass.",
    )
    parser.add_argument(
        "--reference-block",
        type=Path,
        help="Optional reference ma_block_points_proj_2023.parquet to compare.",
    )
    parser.add_argument(
        "--reference-clr",
        type=Path,
        help="Optional reference ma_tract_clr_2023.parquet to compare.",
    )
    parser.add_argument(
        "--reference-strata",
        type=Path,
        help="Optional reference ma_tract_strata_2023.parquet to compare.",
    )
    return parser.parse_args()


def require_file(path: Path) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Required input does not exist: {path}")
    return path


def load_inputs(raw_dir: Path) -> dict[str, Path]:
    raw_dir = raw_dir.expanduser().resolve()
    return {
        "wac": require_file(raw_dir / "ma_wac_S000_JT02_2023.csv"),
        "xwalk": require_file(raw_dir / "ma_xwalk.csv.gz"),
        "tiger_zip": require_file(raw_dir / "tl_2023_25_tract.zip"),
        "county_list": require_file(raw_dir / "list1_2023.xlsx"),
    }


def load_wac(path: Path) -> pd.DataFrame:
    wac = pd.read_csv(path, dtype={"w_geocode": str})
    required = {"w_geocode", "C000", *NAICS_COLUMNS}
    missing = required.difference(wac.columns)
    if missing:
        raise ValueError(f"WAC is missing required columns: {sorted(missing)}")

    wac["w_geocode"] = wac["w_geocode"].str.zfill(15)
    if not wac["w_geocode"].str.fullmatch(r"\d{15}").all():
        raise ValueError("WAC w_geocode values must be 15-digit block GEOIDs")
    if wac["w_geocode"].duplicated().any():
        raise ValueError("WAC contains duplicate w_geocode rows")

    wac["tract"] = wac["w_geocode"].str[:11]
    return wac


def load_tiger(path: Path) -> gpd.GeoDataFrame:
    inner = "tl_2023_25_tract/tl_2023_25_tract.shp"
    tiger_uri = f"zip://{path.as_posix()}!{inner}"
    tiger = gpd.read_file(tiger_uri).rename(columns={"GEOID": "tract"})
    required = {"tract", "ALAND"}
    missing = required.difference(tiger.columns)
    if missing:
        raise ValueError(f"TIGER file is missing required columns: {sorted(missing)}")
    tiger["tract"] = tiger["tract"].astype(str).str.zfill(11)
    tiger["aland_km2"] = tiger["ALAND"] / 1_000_000.0
    return tiger


def make_tract_data(
    wac: pd.DataFrame, tiger: gpd.GeoDataFrame
) -> tuple[pd.DataFrame, list[float], dict[str, int]]:
    tract = wac.groupby("tract", as_index=False)[["C000", *NAICS_COLUMNS]].sum()

    max_sector_residual = (tract[NAICS_COLUMNS].sum(axis=1) - tract["C000"]).abs().max()
    if max_sector_residual != 0:
        raise ValueError(
            "CNS01-CNS19 do not sum to C000; CNS20 is not structural zero "
            f"(maximum residual {max_sector_residual})"
        )

    tract = tract.merge(
        tiger[["tract", "aland_km2"]], on="tract", how="left", validate="one_to_one"
    )
    if tract["aland_km2"].isna().any():
        missing = tract.loc[tract["aland_km2"].isna(), "tract"].tolist()
        raise ValueError(f"Tracts unmatched to TIGER: {missing[:10]}")

    below_jobs = tract["C000"] < MIN_PRIVATE_JOBS
    non_positive_land = tract["aland_km2"] <= 0
    exclusion_counts = {
        "jobs_below_50": int(below_jobs.sum()),
        "non_positive_land_area": int(non_positive_land.sum()),
        "both_conditions": int((below_jobs & non_positive_land).sum()),
        "excluded_unique_tracts": int((below_jobs | non_positive_land).sum()),
    }
    tract = tract.loc[
        (tract["C000"] >= MIN_PRIVATE_JOBS) & (tract["aland_km2"] > 0)
    ].reset_index(drop=True)
    tract["emp_density"] = tract["C000"] / tract["aland_km2"]
    tract["density_q"], bins = pd.qcut(
        tract["emp_density"], 4, labels=DENSITY_LEVELS, retbins=True
    )

    density_counts = {
        str(level): int(count)
        for level, count in tract["density_q"].value_counts().sort_index().items()
    }
    if len(tract) != EXPECTED_TRACTS:
        raise AssertionError(f"Expected {EXPECTED_TRACTS} tracts, found {len(tract)}")
    if int(tract["C000"].sum()) != EXPECTED_PRIVATE_JOBS:
        raise AssertionError(
            f"Expected {EXPECTED_PRIVATE_JOBS} jobs, found {int(tract['C000'].sum())}"
        )
    if density_counts != EXPECTED_DENSITY_COUNTS:
        raise AssertionError(
            f"Unexpected density quartile counts: {density_counts}"
        )

    return tract, [float(value) for value in bins], exclusion_counts


def make_compositions(
    tract: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts = tract.set_index("tract")[NAICS_COLUMNS]
    raw_shares = counts.div(counts.sum(axis=1), axis=0)

    adjusted = counts.astype(float) + 0.5
    smoothed_shares = adjusted.div(adjusted.sum(axis=1), axis=0)
    logged = np.log(smoothed_shares)
    clr = logged.sub(logged.mean(axis=1), axis=0)
    clr.columns = [f"clr_{column}" for column in NAICS_COLUMNS]

    if not np.allclose(raw_shares.sum(axis=1), 1.0, atol=1e-12, rtol=0):
        raise AssertionError("Raw sector shares do not close to one")
    max_clr_sum = float(clr.sum(axis=1).abs().max())
    if max_clr_sum >= 1e-8:
        raise AssertionError(f"CLR rows do not sum to zero: max abs={max_clr_sum}")

    return raw_shares, smoothed_shares, clr


def make_county_lookup(path: Path) -> pd.DataFrame:
    columns = {
        "FIPS State Code": str,
        "FIPS County Code": str,
    }
    source = pd.read_excel(path, sheet_name="List 1", header=2, dtype=columns)
    source = source.loc[source["FIPS State Code"] == STATE_FIPS].copy()
    source["county_fips"] = (
        source["FIPS State Code"] + source["FIPS County Code"]
    )

    def cbsa_stratum(row: pd.Series) -> str:
        if "Boston-Cambridge-Newton" in str(row["CBSA Title"]):
            return "boston_msa"
        if (
            row["Metropolitan/Micropolitan Statistical Area"]
            == "Metropolitan Statistical Area"
        ):
            return "other_msa"
        return "non_metro"

    source["stratum"] = source.apply(cbsa_stratum, axis=1)
    lookup = source[
        ["county_fips", "County/County Equivalent", "CBSA Title", "stratum"]
    ].rename(
        columns={
            "County/County Equivalent": "county_name",
            "CBSA Title": "cbsa_name",
        }
    )
    if len(lookup) != 14 or lookup["county_fips"].duplicated().any():
        raise AssertionError("Expected one lookup row for each of 14 MA counties")
    return lookup


def attach_counties(
    tract: pd.DataFrame, county_lookup: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tract_jobs = tract.copy()
    tract_jobs["county_fips"] = tract_jobs["tract"].str[:5]
    tract_jobs = tract_jobs.merge(
        county_lookup, on="county_fips", how="left", validate="many_to_one"
    )
    if tract_jobs[["county_name", "stratum"]].isna().any().any():
        raise AssertionError("Some retained tracts are missing county metadata")

    # Preserve the exact schema used by the validated downstream analysis.
    strata = tract_jobs[
        ["tract", "county_fips", "county_name", "stratum", "C000"]
    ].copy()
    return tract_jobs, strata


def make_projected_block_points(
    wac: pd.DataFrame, xwalk_path: Path, tract: pd.DataFrame
) -> pd.DataFrame:
    xwalk = pd.read_csv(
        xwalk_path,
        dtype={"tabblk2020": str},
        usecols=["tabblk2020", "blklatdd", "blklondd"],
    )
    xwalk["tabblk2020"] = xwalk["tabblk2020"].str.zfill(15)
    if xwalk["tabblk2020"].duplicated().any():
        raise ValueError("LODES crosswalk contains duplicate tabblk2020 rows")

    blocks = wac.merge(
        xwalk,
        left_on="w_geocode",
        right_on="tabblk2020",
        how="left",
        validate="one_to_one",
    )
    if blocks[["blklatdd", "blklondd"]].isna().any().any():
        raise AssertionError("Some WAC blocks are missing crosswalk coordinates")

    density_lookup = tract.set_index("tract")["density_q"]
    blocks["density_q"] = blocks["tract"].map(density_lookup)
    blocks = blocks.loc[blocks["density_q"].notna()].copy()

    projected = gpd.GeoDataFrame(
        blocks,
        geometry=gpd.points_from_xy(blocks["blklondd"], blocks["blklatdd"]),
        crs="EPSG:4326",
    ).to_crs(TARGET_CRS)
    blocks["x_m"] = projected.geometry.x.to_numpy()
    blocks["y_m"] = projected.geometry.y.to_numpy()

    long = blocks.melt(
        id_vars=["w_geocode", "x_m", "y_m", "tract", "density_q"],
        value_vars=NAICS_COLUMNS,
        var_name="sector",
        value_name="jobs",
    )
    long = long.loc[long["jobs"] > 0].reset_index(drop=True)

    if len(long) != EXPECTED_BLOCK_SECTOR_ROWS:
        raise AssertionError(
            f"Expected {EXPECTED_BLOCK_SECTOR_ROWS} block-sector rows, found {len(long)}"
        )
    if int(long["jobs"].sum()) != EXPECTED_PRIVATE_JOBS:
        raise AssertionError("Block-sector weights do not sum to retained private jobs")
    if not np.isfinite(long[["x_m", "y_m"]].to_numpy()).all():
        raise AssertionError("Projected block coordinates contain non-finite values")

    return long


def compare_reference(actual: pd.DataFrame, reference: Path, label: str) -> None:
    reference_frame = pd.read_parquet(require_file(reference))
    pd.testing.assert_frame_equal(
        actual,
        reference_frame,
        check_exact=True,
        check_categorical=True,
        check_like=False,
    )
    print(f"Reference check passed: {label}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ensure_targets_are_safe(output_dir: Path, force: bool) -> None:
    existing = [output_dir / name for name in OUTPUT_NAMES if (output_dir / name).exists()]
    if existing and not force:
        formatted = "\n".join(f"  - {path}" for path in existing)
        raise FileExistsError(
            "Generated outputs already exist. Use --force to replace them:\n" + formatted
        )


def write_outputs(
    output_dir: Path,
    frames: dict[str, pd.DataFrame],
    summary: dict[str, object],
    force: bool,
) -> None:
    output_dir = output_dir.expanduser().resolve()
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    ensure_targets_are_safe(output_dir, force)

    with tempfile.TemporaryDirectory(
        prefix="erp-preprocess-", dir=output_dir.parent
    ) as temp_name:
        temp_dir = Path(temp_name)
        for filename, frame in frames.items():
            frame.to_parquet(temp_dir / filename)

        summary_path = temp_dir / "preprocessing_summary.json"
        summary_path.write_text(
            json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

        checksum_lines = []
        for path in sorted(temp_dir.iterdir(), key=lambda item: item.name):
            if path.name != "checksums.sha256":
                checksum_lines.append(f"{sha256(path)}  {path.name}")
        (temp_dir / "checksums.sha256").write_text(
            "\n".join(checksum_lines) + "\n", encoding="utf-8"
        )

        output_dir.mkdir(parents=True, exist_ok=True)
        for path in temp_dir.iterdir():
            destination = output_dir / path.name
            if destination.exists() and force:
                destination.unlink()
            shutil.move(str(path), str(destination))


def main() -> None:
    args = parse_args()
    inputs = load_inputs(args.raw_dir)

    wac = load_wac(inputs["wac"])
    tiger = load_tiger(inputs["tiger_zip"])
    tract, density_bins, exclusion_counts = make_tract_data(wac, tiger)
    raw_shares, smoothed_shares, clr = make_compositions(tract)
    county_lookup = make_county_lookup(inputs["county_list"])
    tract_jobs, strata = attach_counties(tract, county_lookup)
    block_points = make_projected_block_points(wac, inputs["xwalk"], tract)

    if args.reference_block:
        compare_reference(block_points, args.reference_block, "block points")
    if args.reference_clr:
        compare_reference(clr, args.reference_clr, "CLR")
    if args.reference_strata:
        compare_reference(strata, args.reference_strata, "strata")

    frames = {
        "ma_tract_jobs_2023.parquet": tract_jobs,
        "ma_tract_shares_2023.parquet": raw_shares,
        "ma_tract_shares_smoothed_2023.parquet": smoothed_shares,
        "ma_tract_clr_2023.parquet": clr,
        "ma_tract_strata_2023.parquet": strata,
        "ma_block_points_proj_2023.parquet": block_points,
    }
    summary = {
        "year": YEAR,
        "job_type": "JT02 (all private jobs)",
        "minimum_private_jobs": MIN_PRIVATE_JOBS,
        "target_crs": TARGET_CRS,
        "tracts": len(tract),
        "private_jobs": int(tract["C000"].sum()),
        "block_sector_rows": len(block_points),
        "density_quartile_counts": EXPECTED_DENSITY_COUNTS,
        "density_cutpoints_jobs_per_km2": density_bins[1:-1],
        "excluded": exclusion_counts,
        "max_abs_clr_row_sum": float(clr.sum(axis=1).abs().max()),
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "geopandas": gpd.__version__,
            "pyarrow": pyarrow.__version__,
            "pyproj": pyproj.__version__,
        },
    }
    write_outputs(args.output_dir, frames, summary, args.force)

    print(f"Wrote {len(frames)} parquet files to {args.output_dir.resolve()}")
    print(
        "Density cut-points (jobs/km²): "
        + ", ".join(f"{value:.6f}" for value in density_bins[1:-1])
    )
    print("All preprocessing invariants passed.")


if __name__ == "__main__":
    main()
