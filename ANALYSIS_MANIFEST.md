# Analysis manifest

This manifest records how the analyses previously stored in separate desktop
folders are represented in the consolidated repository.

| Former folder | Consolidated files |
|---|---|
| `Mfunction` | `scripts/02_mfunction_free.R`, `scripts/03_mfunction_county.R`, M-function figure and robustness scripts, final tables and null arrays |
| `PERMANOVA` | `scripts/04_permanova_permdisp.R`, `scripts/robustness/01_permdisp_leaveout.R`, final location/dispersion tables |
| `Topotest` | `scripts/05_topotest.py`, free and county-restricted result tables, null-distribution figure builder |
| `TDABM` | `scripts/06_ball_mapper.R`, `scripts/07_reconciliation.R`, Ball Mapper uncertainty and epsilon checks |
| `Figures_Output` | `scripts/figures/`, final PNG files, and machine-readable backing tables |

Only final analytical material is included. Session images, shell history,
checkpoint chunks, duplicated implementations, truncated notebooks, document-
editing utilities, draft notes, and backup manuscripts are excluded.

The complete TopoTest implementation is the executable Python script in this
repository. It contains both permutation schemes, independent parallel random
streams, pairwise correction, and null-table export; it supersedes the partial
desktop notebook.
