# Data provenance and files

All inputs are public 2023 United States government statistical or geographic
files. The repository keeps the exact source files used by the analysis so the
preprocessing stage can run without manual downloads.

## Raw inputs

| File | Role |
|---|---|
| `ma_wac_S000_JT02_2023.csv` | Massachusetts LODES WAC, all workers, private primary jobs |
| `ma_xwalk.csv.gz` | LODES block crosswalk and block coordinates |
| `tl_2023_25_tract.zip` | 2023 TIGER/Line Massachusetts census tract geometry and land area |
| `list1_2023.xlsx` | 2023 county/CBSA reference table used for county names and strata |

The WAC and crosswalk files are from the U.S. Census Bureau LEHD Origin-
Destination Employment Statistics release. The tract archive is from the U.S.
Census Bureau TIGER/Line 2023 release.

## Processed outputs

`scripts/01_preprocess.py` creates:

- `ma_tract_jobs_2023.parquet`: tract job and land-area totals;
- `ma_tract_shares_2023.parquet`: raw 19-sector shares;
- `ma_tract_shares_smoothed_2023.parquet`: shares after adding 0.5 to each cell;
- `ma_tract_clr_2023.parquet`: 19 centred log-ratio coordinates;
- `ma_tract_strata_2023.parquet`: density and county lookup;
- `ma_block_points_proj_2023.parquet`: weighted block-sector points in
  EPSG:26986, whose distance unit is metres;
- `preprocessing_summary.json` and `checksums.sha256`: integrity metadata.

Tracts with fewer than 50 private primary jobs or non-positive land area are
excluded before density quartiles are assigned.
