"""Table A.2 and Figure A.3: profiles of nine isolated tracts."""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
PROCESSED_DIR = ROOT / "data" / "processed"
TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

SECTOR_LABELS = {
    "CNS01": "Agriculture", "CNS02": "Mining", "CNS03": "Utilities",
    "CNS04": "Construction", "CNS05": "Manufacturing", "CNS06": "Wholesale",
    "CNS07": "Retail", "CNS08": "Transport/Warehouse", "CNS09": "Information",
    "CNS10": "Finance/Insurance", "CNS11": "Real Estate", "CNS12": "Prof/Sci/Tech",
    "CNS13": "Management", "CNS14": "Admin/Waste", "CNS15": "Education",
    "CNS16": "Health/Social", "CNS17": "Arts/Rec", "CNS18": "Accommodation/Food",
    "CNS19": "Other Services",
}

# ---------------------------------------------------------------------
# 1) Invert CLR -> compositional shares (exact identity, no new data)
# ---------------------------------------------------------------------
clr = pd.read_parquet(PROCESSED_DIR / "ma_tract_clr_2023.parquet")
clr.index = clr.index.astype(str)
clr_cols = [c for c in clr.columns if c.startswith("clr_CNS")]

E = np.exp(clr[clr_cols].values)
SH = E / E.sum(axis=1, keepdims=True)
shares = pd.DataFrame(SH, index=clr.index, columns=[c.replace("clr_", "") for c in clr_cols])
assert np.abs(shares.sum(axis=1) - 1).max() < 1e-9

mean_share = shares.mean()
hhi_all = (shares ** 2).sum(axis=1)

# ---------------------------------------------------------------------
# 2) The nine isolated tracts (from Figure 5.4's original-package run)
# ---------------------------------------------------------------------
iso = pd.read_csv(TABLE_DIR / "fig_5_4_isolated_balls.csv", dtype={"tract": str})
strata_tab = pd.read_parquet(PROCESSED_DIR / "ma_tract_strata_2023.parquet")
strata_tab["tract"] = strata_tab["tract"].astype(str)

iso = iso.merge(strata_tab[["tract", "county_name", "C000"]], on="tract", how="left")
iso["hhi"] = hhi_all.reindex(iso["tract"]).values
iso["hhi_percentile"] = 100 * hhi_all.rank(pct=True).reindex(iso["tract"]).values
# GEOID = 2-digit state + 3-digit county + 6-digit tract number (implied
# decimal after the first 4 digits, e.g. "981800" = tract 9818.00), so the
# Census 9800-series special-use range is the first 4 digits of that
# 6-digit tract code, not the last 4 digits of the GEOID.
tract_code_6 = iso["tract"].str[-6:]
iso["is_9800_series"] = tract_code_6.str[:4].between("9800", "9899")

top3 = []
for t in iso["tract"]:
    s = shares.loc[t].sort_values(ascending=False)
    top3.append("; ".join(f"{SECTOR_LABELS[sec]} {100*v:.1f}%" for sec, v in s.head(3).items()))
iso["top3_sectors"] = top3

iso_out = iso[["tract", "county_name", "density_q", "C000", "hhi", "hhi_percentile",
               "is_9800_series", "top3_sectors"]].sort_values("hhi", ascending=False)
iso_out.to_csv(TABLE_DIR / "table_A2_isolated_tracts.csv", index=False)
print(iso_out.to_string(index=False))


def render_table(df, title, out_png, out_pdf, fontsize=8.5, header_wrap=14):
    import textwrap
    wrapped_cols = ["\n".join(textwrap.wrap(str(c), header_wrap)) for c in df.columns]
    n_rows, n_cols = df.shape
    fig_w = max(9, 1.7 * n_cols)
    fig_h = 0.55 * (n_rows + 1) + 0.9
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    tbl = ax.table(cellText=df.values, colLabels=wrapped_cols, cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)
    tbl.auto_set_column_width(col=list(range(n_cols)))
    tbl.scale(1, 2.1)
    for (row, col), cell in tbl.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if row == 0:
            cell.set_facecolor("#2c3e50")
            cell.set_text_props(color="white", weight="bold")
        else:
            cell.set_facecolor("#f2f2f2" if row % 2 == 0 else "white")
            cell.set_text_props(wrap=True)
    ax.set_title(title, fontsize=13, weight="bold", pad=18)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


iso_fmt = iso_out.copy()
iso_fmt["hhi"] = iso_fmt["hhi"].map(lambda v: f"{v:.3f}")
iso_fmt["hhi_percentile"] = iso_fmt["hhi_percentile"].map(lambda v: f"{v:.1f}")
iso_fmt = iso_fmt.rename(columns={
    "tract": "GEOID", "county_name": "County", "density_q": "Density q",
    "C000": "Jobs (C000)", "hhi": "HHI", "hhi_percentile": "HHI pctile",
    "is_9800_series": "9800-series", "top3_sectors": "Top 3 sectors (share)",
})
render_table(iso_fmt, "Table A.2. The nine Ball Mapper isolated tracts: full attributes",
             FIGURE_DIR / "table_A2_isolated_tracts.png",
             FIGURE_DIR / "table_A2_isolated_tracts.pdf")
print("wrote table_A2_isolated_tracts.png/.pdf")

# spot-check the paper's specific callouts (ERP Sec 5.4)
def share_of(tract, sector):
    return shares.loc[tract, sector]

checks = [
    ("25025981800", "CNS15", 0.946),  # Education
    ("25025981300", "CNS08", 0.743),  # Transport/Warehousing
    ("25025081001", "CNS16", 0.693),  # Health
    ("25025981900", "CNS03", 0.656),  # Utilities
]
print("\nverification against ERP Sec 5.4 callouts:")
for tract, sector, target in checks:
    v = share_of(tract, sector)
    print(f"  {tract} {SECTOR_LABELS[sector]}: {v*100:.1f}% (target {target*100:.1f}%)")
    assert abs(v - target) < 0.01
c000_check = strata_tab.set_index("tract").loc["25025081001", "C000"]
print(f"  25025081001 C000 = {c000_check:,} (target 48,020)")
assert c000_check == 48020
n_9800 = int(iso["is_9800_series"].sum())
print(f"  9800-series tracts among the nine: {n_9800} (target 4)")
assert n_9800 == 4

# ---------------------------------------------------------------------
# 3) Figure A.3 -- 9 (tracts) x 19 (sectors) share heatmap
# ---------------------------------------------------------------------
all_sectors = [f"CNS{i:02d}" for i in range(1, 20)]
mat = shares.loc[iso_out["tract"], all_sectors].values

fig, ax = plt.subplots(figsize=(11, 5.5))
im = ax.imshow(mat, cmap="YlOrRd", aspect="auto", vmin=0, vmax=1)
ax.set_xticks(range(len(all_sectors)))
ax.set_xticklabels([f"{s}\n{SECTOR_LABELS[s]}" for s in all_sectors], rotation=90, fontsize=7)
row_labels = [f"{t}  ({dq}, HHI={h:.2f}{'*' if is9800 else ''})"
              for t, dq, h, is9800 in zip(iso_out["tract"], iso_out["density_q"],
                                           iso_out["hhi"], iso_out["is_9800_series"])]
ax.set_yticks(range(len(iso_out)))
ax.set_yticklabels(row_labels, fontsize=8)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if mat[i, j] > 0.10:
            ax.text(j, i, f"{100*mat[i,j]:.0f}", ha="center", va="center", fontsize=6,
                    color="white" if mat[i, j] > 0.5 else "black")
ax.set_title("Sector employment shares for the nine Ball Mapper isolated tracts\n"
              "(* = Census 9800-series special-use tract)")
cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label("employment share")
fig.tight_layout()
fig.savefig(FIGURE_DIR / "fig_A3_isolated_tract_profiles.png", dpi=300)
fig.savefig(FIGURE_DIR / "fig_A3_isolated_tract_profiles.pdf")
plt.close(fig)
print("\nsaved fig_A3_isolated_tract_profiles.png/.pdf and table_A2_isolated_tracts.csv")
