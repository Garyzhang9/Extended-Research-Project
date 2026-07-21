"""Figure 3.1: Massachusetts employment-density quartile map."""
from pathlib import Path

import geopandas as gpd
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[2]
PROCESSED_DIR = ROOT / "data" / "processed"
RAW_DIR = ROOT / "data" / "raw"
TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
STRATA = ["Q1_low", "Q2", "Q3", "Q4_high"]
STRATA_COLOR = {"Q1_low": "#2c7bb6", "Q2": "#abd9e9", "Q3": "#fdae61", "Q4_high": "#d7191c"}
CUTPOINTS = (93.72, 329.87, 997.77)

# ---------------------------------------------------------------------
# 1) Load tract boundaries (GEOID kept as string, never numeric)
# ---------------------------------------------------------------------
tiger_zip = RAW_DIR / "tl_2023_25_tract.zip"
gdf = gpd.read_file(f"zip://{tiger_zip}!tl_2023_25_tract/tl_2023_25_tract.shp")
gdf["GEOID"] = gdf["GEOID"].astype(str).str.zfill(11)
assert gdf["GEOID"].str.len().eq(11).all()
print(f"loaded {len(gdf)} MA tract polygons")

# ---------------------------------------------------------------------
# 2) Load density_q labels and join
# ---------------------------------------------------------------------
blk = pd.read_parquet(PROCESSED_DIR / "ma_block_points_proj_2023.parquet")
dq = blk[["tract", "density_q"]].drop_duplicates()
dq["tract"] = dq["tract"].astype(str).str.zfill(11)
assert len(dq) == 1598

merged = gdf.merge(dq, left_on="GEOID", right_on="tract", how="inner")
assert len(merged) == 1598, f"expected 1598 matched tracts, got {len(merged)}"
merged["density_q"] = pd.Categorical(merged["density_q"], categories=STRATA, ordered=True)
print(f"joined {len(merged)} tracts (of {len(gdf)} shapefile tracts; "
      f"{len(gdf) - len(merged)} excluded: <50 private jobs or zero land area)")

# ---------------------------------------------------------------------
# 3) County outlines dissolved directly from the tract shapefile
# ---------------------------------------------------------------------
counties = gdf.dissolve(by=["STATEFP", "COUNTYFP"])
assert len(counties) == 14

# ---------------------------------------------------------------------
# 4) Verification: Suffolk (25025) Q4_high, Barnstable (25001) Q1_low
# ---------------------------------------------------------------------
merged["county_fips"] = merged["GEOID"].str[:5]
suffolk = merged[merged["county_fips"] == "25025"]
barnstable = merged[merged["county_fips"] == "25001"]
n_suffolk_q4 = (suffolk["density_q"] == "Q4_high").sum()
n_barnstable_q1 = (barnstable["density_q"] == "Q1_low").sum()
print(f"Suffolk: {n_suffolk_q4} / {len(suffolk)} tracts in Q4_high (target 137 / 225)")
print(f"Barnstable: {n_barnstable_q1} / {len(barnstable)} tracts in Q1_low (target 39 / 56)")
assert (n_suffolk_q4, len(suffolk)) == (137, 225)
assert (n_barnstable_q1, len(barnstable)) == (39, 56)

merged.to_file(TABLE_DIR / "fig_3_1_tract_density_quartiles.geojson", driver="GeoJSON")
merged.drop(columns="geometry").to_csv(TABLE_DIR / "fig_3_1_tract_density_quartiles.csv", index=False)

# ---------------------------------------------------------------------
# 5) Render map (300dpi PNG + vector PDF)
# ---------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))
merged.plot(color=merged["density_q"].map(STRATA_COLOR), linewidth=0, ax=ax)
counties.boundary.plot(ax=ax, color="#4d4d4d", linewidth=0.6, alpha=0.8)

legend_handles = [Patch(facecolor=STRATA_COLOR[s], label=s) for s in STRATA]
ax.legend(handles=legend_handles, title="Employment-density\nquartile", loc="lower right", fontsize=9)
ax.set_axis_off()
ax.set_title("Massachusetts employment-density quartiles (2023 LODES WAC, private jobs)")
fig.text(0.5, 0.01,
         f"Cut-points: {CUTPOINTS[0]}, {CUTPOINTS[1]}, {CUTPOINTS[2]} jobs km$^{{-2}}$ "
         "(private jobs per sq. km land area). County boundaries in grey.",
         ha="center", fontsize=8)

fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(FIGURE_DIR / "fig_3_1_density_map.png", dpi=300)
fig.savefig(FIGURE_DIR / "fig_3_1_density_map.pdf")
plt.close(fig)
print("saved fig_3_1_density_map.png / .pdf")
