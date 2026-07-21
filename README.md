# Extended Research Project — reproducible analysis

This repository reproduces the data preparation, statistical analyses, and
figures for *Do Distance-Based and Compositional–Topological Measures of
Industrial Co-location Agree? Evidence from Massachusetts LODES*.

The repository is self-contained: it includes the four public source files, a
validated processed-data snapshot, executable analysis scripts, final result
tables, and final figures. Every script resolves paths from the repository
root; no user-specific absolute paths are used.

## Analysis inventory

| Stage | Script | Main output |
|---|---|---|
| Data preparation | `scripts/01_preprocess.py` | Six parquet files and integrity summary |
| M-function, free permutations | `scripts/02_mfunction_free.R` | Sector-level D statistics and null distribution |
| M-function, county-restricted | `scripts/03_mfunction_county.R` | Restricted p-values and comparison table |
| PERMANOVA and PERMDISP | `scripts/04_permanova_permdisp.R` | Omnibus and pairwise tables |
| Euler-curve TopoTest | `scripts/05_topotest.py` | Free and county-restricted pairwise tables |
| Ball Mapper | `scripts/06_ball_mapper.R` | Graph, communities, isolated balls, ARI and NMI |
| Reconciliation | `scripts/07_reconciliation.R` | Sector-level Spearman comparison |

Additional checks live in `scripts/robustness/`; publication figure and table
builders live in `scripts/figures/`.

## Environment

Python package versions are recorded in `environment.yml`:

```bash
conda env create -f environment.yml
conda activate erp-repro
```

R 4.3.1 and the exact package versions used for the reference run are listed
in `R-packages.txt`. Verify an R installation with:

```bash
Rscript scripts/check_r_environment.R
```

## Quick validation

From the repository root:

```bash
python scripts/01_preprocess.py --force
pytest -q
python -m py_compile scripts/*.py scripts/figures/*.py
Rscript -e 'for (f in list.files("scripts", "[.]R$", recursive=TRUE, full.names=TRUE)) parse(f)'
```

The preprocessing stage checks the 1,598-tract sample, 3,203,251 jobs,
400/399/399/400 quartile sizes, 119,181 block-sector rows, finite projected
coordinates, and CLR closure.

## Full analysis

The computational settings used for the archived results are:

- M-function: 9,999 permutations, radius 1,000 m, seed 20260612, four workers;
- PERMANOVA/PERMDISP: 9,999 permutations, seed 14196142;
- TopoTest: 9,999 permutations per pair, seed 14196142, four workers;
- Ball Mapper: epsilon 7.0 and seed 1;
- regional null: 1,999 permutations; landmark bootstrap: 2,000 runs.

Run stages individually in numeric order. The two longest stages are the
M-function permutations and the TopoTest. For a fast execution check, set a
smaller repetition count without changing the archived results:

```bash
ERP_PERMUTATIONS=19 ERP_CORES=2 Rscript scripts/02_mfunction_free.R
ERP_PERMUTATIONS=19 ERP_JOBS=2 python scripts/05_topotest.py --mode free
```

Generated outputs go to `results/tables/`, `results/figures/`, and
`results/null_distributions/`. Temporary restart files go to the ignored
`results/checkpoints/` directory.

## Repository contents

- `data/raw/`: public input files used by preprocessing;
- `data/processed/`: deterministic processed datasets;
- `results/tables/`: machine-readable final results and figure backing data;
- `results/figures/`: final PNG figures and tables;
- `results/null_distributions/`: archived long-run M-function null arrays;
- `tests/`: processed-data invariants.

See `data/README.md` for provenance and `REPRODUCIBILITY.md` for the execution
map and expected checkpoints.
