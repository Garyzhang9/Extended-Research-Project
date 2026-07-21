"""Table A.1: county by employment-density quartile cross-tabulation."""
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import textwrap

ROOT = Path(__file__).resolve().parents[2]
PROCESSED_DIR = ROOT / "data" / "processed"
TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def render_table(df, title, out_png, out_pdf, fontsize=9.5, header_wrap=16, index_label=None):
    cols = ([index_label] if index_label else []) + list(df.columns)
    cell_text = df.reset_index().values if index_label else df.values
    wrapped_cols = ["\n".join(textwrap.wrap(str(c), header_wrap)) for c in cols]
    n_rows, n_cols = len(df), len(cols)
    fig_w = max(9, 1.4 * n_cols)
    fig_h = 0.5 * (n_rows + 1) + 0.9
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    tbl = ax.table(cellText=cell_text, colLabels=wrapped_cols, cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)
    tbl.auto_set_column_width(col=list(range(n_cols)))
    tbl.scale(1, 1.9)
    for (row, col), cell in tbl.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if row == 0:
            cell.set_facecolor("#2c3e50")
            cell.set_text_props(color="white", weight="bold")
        else:
            cell.set_facecolor("#f2f2f2" if row % 2 == 0 else "white")
    ax.set_title(title, fontsize=13, weight="bold", pad=18)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
STRATA = ["Q1_low", "Q2", "Q3", "Q4_high"]

blk = pd.read_parquet(PROCESSED_DIR / "ma_block_points_proj_2023.parquet")
dq = blk[["tract", "density_q"]].drop_duplicates()
dq["tract"] = dq["tract"].astype(str)

strata_tab = pd.read_parquet(PROCESSED_DIR / "ma_tract_strata_2023.parquet")
strata_tab["tract"] = strata_tab["tract"].astype(str)

dat = dq.merge(strata_tab[["tract", "county_fips", "county_name"]], on="tract", how="inner")
assert len(dat) == 1598

ct = pd.crosstab(dat["county_name"], dat["density_q"])[STRATA]
ct["Total"] = ct.sum(axis=1)
ct["n_quartiles_present"] = (ct[STRATA] > 0).sum(axis=1)
ct = ct.sort_values("Total", ascending=False)

print(ct.to_string())
assert len(ct) == 14
assert (ct["n_quartiles_present"] >= 2).all()

suffolk = ct.loc["Suffolk County"]
barnstable = ct.loc["Barnstable County"]
print(f"\nSuffolk: {suffolk['Q4_high']} / {suffolk['Total']} in Q4_high (target 137/225)")
print(f"Barnstable: {barnstable['Q1_low']} / {barnstable['Total']} in Q1_low (target 39/56)")
assert (suffolk["Q4_high"], suffolk["Total"]) == (137, 225)
assert (barnstable["Q1_low"], barnstable["Total"]) == (39, 56)

ct.to_csv(TABLE_DIR / "table_A1_county_by_quartile.csv")
print(f"\nwrote table_A1_county_by_quartile.csv ({len(ct)} rows)")

ct_display = ct.rename(columns={"n_quartiles_present": "Quartiles present"})
render_table(ct_display, "Table A.1. Tract counts by county and employment-density quartile",
             FIGURE_DIR / "table_A1_county_by_quartile.png",
             FIGURE_DIR / "table_A1_county_by_quartile.pdf",
             index_label="County")
print("wrote table_A1_county_by_quartile.png/.pdf")
