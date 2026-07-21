#!/usr/bin/env python3
"""Table 3.1: sample and employment-density stratification."""

from pathlib import Path
import textwrap

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
PROCESSED_DIR = ROOT / "data" / "processed"
TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
STRATA = ["Q1_low", "Q2", "Q3", "Q4_high"]


def render_table(df: pd.DataFrame, title: str, out_png: Path, out_pdf: Path) -> None:
    columns = ["\n".join(textwrap.wrap(str(column), 16)) for column in df.columns]
    figure, axis = plt.subplots(figsize=(max(9, 1.55 * len(columns)), 3.7))
    axis.axis("off")
    table = axis.table(
        cellText=df.values, colLabels=columns, cellLoc="center", loc="center"
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.auto_set_column_width(col=list(range(len(columns))))
    table.scale(1, 2.0)
    for (row, _), cell in table.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if row == 0:
            cell.set_facecolor("#2c3e50")
            cell.set_text_props(color="white", weight="bold")
        else:
            cell.set_facecolor("#f2f2f2" if row % 2 == 0 else "white")
    axis.set_title(title, fontsize=13, weight="bold", pad=18)
    figure.tight_layout()
    figure.savefig(out_png, dpi=300, bbox_inches="tight")
    figure.savefig(out_pdf, bbox_inches="tight")
    plt.close(figure)


tracts = pd.read_parquet(PROCESSED_DIR / "ma_tract_jobs_2023.parquet")
tracts["density_q"] = pd.Categorical(
    tracts["density_q"], categories=STRATA, ordered=True
)
assert len(tracts) == 1_598
assert int(tracts["C000"].sum()) == 3_203_251
assert tracts["density_q"].value_counts().sort_index().tolist() == [400, 399, 399, 400]

cutpoints = (
    tracts.groupby("density_q", observed=True)["emp_density"]
    .max()
    .reindex(STRATA[:-1])
    .to_numpy()
)
assert all(abs(value - target) < 0.01 for value, target in zip(
    cutpoints, (93.718079, 329.871241, 997.766017)
))

rows = []
for stratum in STRATA:
    subset = tracts[tracts["density_q"] == stratum]
    main_counties = ", ".join(
        f"{name.replace(' County', '')} ({count})"
        for name, count in subset["county_name"].value_counts().head(3).items()
    )
    rows.append({
        "Density stratum": stratum,
        "Tracts": len(subset),
        "Density range (jobs/km2)": (
            f"{subset['emp_density'].min():.2f}-{subset['emp_density'].max():.2f}"
        ),
        "Median density (jobs/km2)": round(subset["emp_density"].median(), 2),
        "Total jobs": int(subset["C000"].sum()),
        "Share of state jobs (%)": round(100 * subset["C000"].sum() / 3_203_251, 1),
        "Main counties (tract count)": main_counties,
    })

result = pd.DataFrame(rows)
result.to_csv(TABLE_DIR / "table_3_1_sample_description.csv", index=False)
render_table(
    result,
    "Table 3.1. Sample and employment-density stratification description",
    FIGURE_DIR / "table_3_1_sample_description.png",
    FIGURE_DIR / "table_3_1_sample_description.pdf",
)
print(result.to_string(index=False))
