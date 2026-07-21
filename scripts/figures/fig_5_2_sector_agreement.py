"""Figure 5.2 and Table 5.2: sector-level measurement comparison."""
import subprocess
import io
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

np.random.seed(14196142)  # no stochastic step in this script; set for consistency with the pipeline's seed

ROOT = Path(__file__).resolve().parents[2]
PROCESSED_DIR = ROOT / "data" / "processed"
TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

SECTOR_NAMES = {
    "CNS01": "Agriculture, Forestry, Fishing and Hunting",
    "CNS03": "Utilities",
    "CNS04": "Construction",
    "CNS05": "Manufacturing",
    "CNS06": "Wholesale Trade",
    "CNS07": "Retail Trade",
    "CNS08": "Transportation and Warehousing",
    "CNS09": "Information",
    "CNS10": "Finance and Insurance",
    "CNS11": "Real Estate and Rental and Leasing",
    "CNS12": "Professional, Scientific and Technical Services",
    "CNS13": "Management of Companies and Enterprises",
    "CNS14": "Administrative and Support and Waste Management",
    "CNS15": "Educational Services",
    "CNS16": "Health Care and Social Assistance",
    "CNS17": "Arts, Entertainment and Recreation",
    "CNS18": "Accommodation and Food Services",
    "CNS19": "Other Services (except Public Administration)",
}

# ---------------------------------------------------------------------
# 1) CLR coordinates + density_q labels (tract-level, raw parquet)
# ---------------------------------------------------------------------
clr = pd.read_parquet(PROCESSED_DIR / "ma_tract_clr_2023.parquet")
clr.index = clr.index.astype(str)
clr.index.name = "tract"

blk = pd.read_parquet(PROCESSED_DIR / "ma_block_points_proj_2023.parquet")
dq = blk[["tract", "density_q"]].drop_duplicates()
dq["tract"] = dq["tract"].astype(str)

clr = clr.join(dq.set_index("tract"), how="inner")
assert len(clr) == 1598, f"expected 1598 tracts, got {len(clr)}"
assert clr["density_q"].isna().sum() == 0

clr_cols = [c for c in clr.columns if c.startswith("clr_CNS")]
assert len(clr_cols) == 19, f"expected 19 CLR coordinates, got {len(clr_cols)}"

# ---------------------------------------------------------------------
# 2) Exact between-stratum sum-of-squares decomposition (no permutation)
# ---------------------------------------------------------------------
grand_mean = clr[clr_cols].mean()
ss_between = pd.Series(
    {cn: sum(len(sub) * (sub[cn].mean() - grand_mean[cn]) ** 2
             for _, sub in clr.groupby("density_q", observed=True))
     for cn in clr_cols}
)
total_ss = ss_between.sum()
print(f"total between-stratum SS (all 19 CLR coords) = {total_ss:.4f}  "
      f"(must equal PERMANOVA density_q SumOfSqs = 2714.882)")
assert abs(total_ss - 2714.882) < 0.01

contrib = ss_between / total_ss  # composition_contribution, sums to 1 over 19 coords
assert abs(contrib.sum() - 1.0) < 1e-9

# ---------------------------------------------------------------------
# 3) Merge with the M-function Ds pipeline outputs (18 sectors; CNS02 has
#    no Ds because it fails the >=30-occupied-blocks rule and drops out
#    automatically via the inner join; CNS20 never enters clr_cols)
# ---------------------------------------------------------------------
rq1_free = pd.read_csv(TABLE_DIR / "rq1_results.csv")
rq1_restr = pd.read_csv(TABLE_DIR / "rq1_county_restricted_results.csv")

# ---------------------------------------------------------------------
# 3b) Recompute D_obs at full double precision (dbmss::Mhat, r*=1000m,
#     identical method to RQ1_Mfunction.R) so that sectors which happen to
#     round to the same 4dp value in rq1_results.csv (CNS18/CNS11/CNS15,
#     all "0.0079") are not falsely treated as exact ties.
# ---------------------------------------------------------------------
R_SNIPPET = r'''
suppressMessages({library(arrow); library(dbmss)})
pts <- read_parquet("__BLOCK_FILE__")
rq1 <- read.csv("__RESULT_FILE__")
Family1 <- rq1$sector
strata_levels <- c("Q1_low","Q2","Q3","Q4_high")
W <- spatstat.geom::owin(xrange=range(pts$x_m), yrange=range(pts$y_m))
P_by_q <- lapply(strata_levels, function(q){
  keep <- pts$density_q==q
  wmppp(data.frame(X=pts$x_m[keep],Y=pts$y_m[keep],PointType=pts$sector[keep],PointWeight=pts$jobs[keep]), window=W)
})
names(P_by_q) <- strata_levels
varpop <- function(v) mean((v-mean(v))^2)
D_precise <- sapply(Family1, function(s){
  Ms <- sapply(strata_levels, function(q) Mhat(P_by_q[[q]], ReferenceType=s, NeighborType=s, r=c(0,1000))$M[2])
  varpop(Ms)
})
write.csv(data.frame(sector=Family1, D_obs=D_precise), stdout(), row.names=FALSE)
'''
R_SNIPPET = (R_SNIPPET
             .replace("__BLOCK_FILE__", (PROCESSED_DIR / "ma_block_points_proj_2023.parquet").as_posix())
             .replace("__RESULT_FILE__", (TABLE_DIR / "rq1_results.csv").as_posix()))
res = subprocess.run(["Rscript", "-e", R_SNIPPET], capture_output=True, text=True, check=True)
d_precise = pd.read_csv(io.StringIO(res.stdout))
assert len(d_precise) == 18
# sanity: full-precision values must round to the same 4dp as the pipeline's own file
_chk = d_precise.merge(rq1_free[["sector", "D_obs"]], on="sector", suffixes=("_precise", "_pipeline"))
assert (_chk["D_obs_precise"].round(4) == _chk["D_obs_pipeline"]).all()
print("full-precision D_obs recomputation matches rq1_results.csv to 4dp for all 18 sectors")
for s in ["CNS18", "CNS11", "CNS15"]:
    print(f"  {s}: D_obs = {d_precise.loc[d_precise.sector == s, 'D_obs'].iloc[0]:.7f}")

comp = pd.DataFrame({
    "sector": [c.replace("clr_", "") for c in clr_cols],
    "composition_contribution": contrib.values,
})

table = (
    d_precise[["sector", "D_obs"]]
    .merge(rq1_free[["sector", "p_raw", "p_BH"]], on="sector", how="left")
    .merge(rq1_restr[["sector", "p_restr", "p_restr_BH"]], on="sector", how="left")
    .merge(comp, on="sector", how="inner")
    .rename(columns={"p_restr": "p_restricted", "p_restr_BH": "p_restricted_BH"})
)
assert len(table) == 18, f"expected 18 sectors, got {len(table)}"
assert "CNS02" not in table["sector"].values
assert "CNS20" not in table["sector"].values
contrib_by_sector = contrib.rename(index=lambda c: c.replace("clr_", ""))
print(f"composition_contribution sums to {table['composition_contribution'].sum():.4f} over the "
      f"18 sectors with a Ds statistic (full 19-coordinate contrib sums to {contrib.sum():.4f}; "
      f"CNS02 -- excluded here because it has no Ds -- carries the remaining "
      f"{contrib_by_sector.get('CNS02', 0):.4f})")

table["industry_name"] = table["sector"].map(SECTOR_NAMES)
table["D_rank"] = table["D_obs"].rank(ascending=False, method="min").astype(int)
table["composition_rank"] = table["composition_contribution"].rank(ascending=False, method="min").astype(int)
table = table[["sector", "industry_name", "D_obs", "p_raw", "p_BH",
               "p_restricted", "p_restricted_BH", "composition_contribution",
               "D_rank", "composition_rank"]].sort_values("D_rank")

table.to_csv(TABLE_DIR / "table_5_2_sector_results.csv", index=False)
print(f"wrote table_5_2_sector_results.csv  ({len(table)} rows)")


def render_table(df, title, out_png, out_pdf, fontsize=8.5, header_wrap=13):
    import textwrap
    wrapped_cols = ["\n".join(textwrap.wrap(str(c), header_wrap)) for c in df.columns]
    n_rows, n_cols = df.shape
    fig_w = max(9, 1.55 * n_cols)
    fig_h = 0.42 * (n_rows + 1) + 0.9
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    tbl = ax.table(cellText=df.values, colLabels=wrapped_cols, cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)
    tbl.auto_set_column_width(col=list(range(n_cols)))
    tbl.scale(1, 1.7)
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


table_fmt = table.copy()
# 6 significant figures throughout, so near-tied small values (CNS18/CNS11/
# CNS15, all ~0.0079) are shown as genuinely distinct rather than collapsing
# back to the same-looking "0.0079" the pipeline's 4dp storage produced.
table_fmt["D_obs"] = table_fmt["D_obs"].map(lambda v: f"{v:.6g}")
for c in ["p_raw", "p_BH", "p_restricted", "p_restricted_BH"]:
    table_fmt[c] = table_fmt[c].map(lambda v: f"{v:.4f}")
table_fmt["composition_contribution"] = table_fmt["composition_contribution"].map(lambda v: f"{v:.4f}")
table_fmt = table_fmt.rename(columns={
    "sector": "Sector", "industry_name": "Industry", "D_obs": "D_obs",
    "p_raw": "p_raw", "p_BH": "p_BH", "p_restricted": "p_restr",
    "p_restricted_BH": "p_restr_BH", "composition_contribution": "Comp. contrib.",
    "D_rank": "D rank", "composition_rank": "Comp. rank",
})
render_table(table_fmt, "Table 5.2. Eighteen-sector distance-based and compositional results",
             FIGURE_DIR / "table_5_2_sector_results.png",
             FIGURE_DIR / "table_5_2_sector_results.pdf")
print("wrote table_5_2_sector_results.png/.pdf")

# ---------------------------------------------------------------------
# 4) Spearman rank correlation (sector-level agreement test, ERP Sec 4.5)
# ---------------------------------------------------------------------
rho, p_two_sided = stats.spearmanr(table["D_obs"], table["composition_contribution"])
p_one_sided = p_two_sided / 2 if rho > 0 else 1 - p_two_sided / 2
print(f"Spearman rho = {rho:.4f}  one-tailed p = {p_one_sided:.4f}  "
      f"(full-precision D_obs, no artificial ties; paper's own 4dp-rounded "
      f"reading gives rho=0.1499, p=0.2763 -- both well below the rho>=0.50 criterion)")
assert abs(rho - 0.1517) < 1e-3
assert abs(p_one_sided - 0.2740) < 1e-3

# ---------------------------------------------------------------------
# 5) Figure 5.2 -- rank scatter with 45-degree agreement line
# ---------------------------------------------------------------------
highlight = ["CNS01", "CNS03", "CNS04", "CNS05", "CNS12"]  # Agriculture, Utilities,
                                                            # Construction, Manufacturing,
                                                            # Professional Services
fig, ax = plt.subplots(figsize=(6.5, 6.5))

sig = table["p_BH"] < 0.05
ax.scatter(table.loc[~sig, "D_rank"], table.loc[~sig, "composition_rank"],
           s=70, facecolors="white", edgecolors="black", linewidths=1.2,
           label="not significant (free-permutation BH)", zorder=3)
ax.scatter(table.loc[sig, "D_rank"], table.loc[sig, "composition_rank"],
           s=70, facecolors="black", edgecolors="black", linewidths=1.2,
           label="significant (free-permutation BH)", zorder=3)

lims = [0.5, 18.5]
ax.plot(lims, lims, linestyle="--", color="grey", linewidth=1.3, zorder=1,
        label="45° agreement line")

for _, row in table[table["sector"].isin(highlight)].iterrows():
    ax.annotate(f"{row['sector']}\n{row['industry_name'].split(',')[0].split(' and ')[0]}",
                (row["D_rank"], row["composition_rank"]),
                textcoords="offset points", xytext=(7, 6), fontsize=8, zorder=4)

ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel("$D_s$ rank (distance-based concentration, 1 = strongest)")
ax.set_ylabel("Composition-contribution rank (1 = largest CLR between-stratum share)")
ax.set_title("Sector-level agreement: distance-based vs.\ncompositional measures of industrial co-location")
ax.set_aspect("equal")
ax.legend(loc="lower right", fontsize=8, framealpha=0.9)

textstr = (f"Spearman $\\rho$ = {rho:.3f}\n"
           f"one-tailed $p$ = {p_one_sided:.3f}\n"
           f"agreement criterion: $\\rho \\geq 0.50$")
ax.text(0.03, 0.97, textstr, transform=ax.transAxes, fontsize=9,
        va="top", ha="left", bbox=dict(boxstyle="round", facecolor="white", edgecolor="grey"))

fig.tight_layout()
fig.savefig(FIGURE_DIR / "fig_5_2_sector_agreement.png", dpi=300)
fig.savefig(FIGURE_DIR / "fig_5_2_sector_agreement.pdf")
plt.close(fig)
print("saved fig_5_2_sector_agreement.png / .pdf")
