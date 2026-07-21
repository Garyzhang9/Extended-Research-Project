"""Figure A.4: free and county-restricted TopoTest null distributions."""
import os
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from joblib import Parallel, delayed
import gudhi
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
PROCESSED_DIR = ROOT / "data" / "processed"
TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
SEED = 14196142
B = int(os.environ.get("ERP_PERMUTATIONS", "9999"))
N_JOBS = int(os.environ.get("ERP_JOBS", "4"))

# ---------------------------------------------------------------------
# 1) Load CLR point clouds + density_q + county_fips (same 3 raw files
#    and same merge logic as the notebook)
# ---------------------------------------------------------------------
clr = pd.read_parquet(PROCESSED_DIR / "ma_tract_clr_2023.parquet").reset_index()
clr_cols = [c for c in clr.columns if c.startswith("clr_CNS")]

dq = (pd.read_parquet(PROCESSED_DIR / "ma_block_points_proj_2023.parquet")[["tract", "density_q"]]
      .drop_duplicates(subset=["tract"]))
cty = pd.read_parquet(PROCESSED_DIR / "ma_tract_strata_2023.parquet")[["tract", "county_fips"]]

clr["tract"] = clr["tract"].astype(str)
dq["tract"] = dq["tract"].astype(str)
cty["tract"] = cty["tract"].astype(str)

dat = clr.merge(dq, on="tract", how="left").merge(cty, on="tract", how="left")
assert len(dat) == 1598
assert dat["density_q"].notna().all() and dat["county_fips"].notna().all()
assert np.abs(dat[clr_cols].sum(axis=1)).max() < 1e-8

ORDER = ["Q1_low", "Q2", "Q3", "Q4_high"]
clouds = {q: dat.loc[dat["density_q"] == q, clr_cols].to_numpy(dtype=float) for q in ORDER}

# ---------------------------------------------------------------------
# 2) ECC / two-sample D / permutation machinery (verbatim from the notebook)
# ---------------------------------------------------------------------
def ecc_on_grid(X, grid, max_dim=2):
    Dm = pd.DataFrame(X).values
    from scipy.spatial.distance import squareform
    Dm = squareform(pdist(X))
    rips = gudhi.RipsComplex(distance_matrix=Dm, max_edge_length=float(grid[-1]))
    st = rips.create_simplex_tree(max_dimension=max_dim)
    simplices = [(len(s) - 1, f) for s, f in st.get_filtration()]
    dims = np.array([d for d, _ in simplices], dtype=int)
    filt = np.array([f for _, f in simplices], dtype=float)
    sign = (-1.0) ** dims
    return np.array([sign[filt <= r].sum() for r in grid])

def topo_D(X, Y, grid, max_dim=2):
    cx = ecc_on_grid(X, grid, max_dim) / len(X)
    cy = ecc_on_grid(Y, grid, max_dim) / len(Y)
    return np.max(np.abs(cx - cy)), cx, cy

def common_grid(X, Y, q=0.10, n_steps=200):
    Z = np.vstack([X, Y])
    dz = pdist(Z)
    return np.linspace(0.0, float(np.quantile(dz, q)), n_steps)

def _perm_one_D(seed_p, Z, m, N, grid, max_dim):
    rng = np.random.default_rng(seed_p)
    idx = rng.permutation(N)
    Xp, Yp = Z[idx[:m]], Z[idx[m:]]
    return topo_D(Xp, Yp, grid, max_dim)[0]

def topo_twosample_perm_par(X, Y, n_perm=9999, q=0.10, n_steps=200, max_dim=2, seed=SEED, n_jobs=N_JOBS):
    Z = np.vstack([X, Y]); m, N = len(X), len(X) + len(Y)
    grid = np.linspace(0.0, float(np.quantile(pdist(Z), q)), n_steps)
    D_obs = topo_D(X, Y, grid, max_dim)[0]
    child_seeds = np.random.SeedSequence(seed).spawn(n_perm)
    D_null = Parallel(n_jobs=n_jobs)(delayed(_perm_one_D)(s, Z, m, N, grid, max_dim) for s in child_seeds)
    D_null = np.asarray(D_null)
    pval = (np.sum(D_null >= D_obs) + 1) / (n_perm + 1)
    return D_obs, pval, D_null

def county_resplit(rng, county_codes, m_per_county):
    mask = np.zeros(len(county_codes), dtype=bool)
    for c in np.unique(county_codes):
        pos = np.where(county_codes == c)[0]
        k = m_per_county[c]
        chosen = rng.permutation(pos)[:k]
        mask[chosen] = True
    return mask

def _perm_one_D_restr(seed_p, Z, cc, m_per_county, grid, max_dim):
    rng = np.random.default_rng(seed_p)
    mask = county_resplit(rng, cc, m_per_county)
    return topo_D(Z[mask], Z[~mask], grid, max_dim)[0]

def topo_twosample_restr(X, Y, county_X, county_Y, n_perm=9999, q=0.10, n_steps=200, max_dim=2, seed=SEED, n_jobs=N_JOBS):
    Z = np.vstack([X, Y])
    cc = np.concatenate([county_X, county_Y])
    is_X = np.concatenate([np.ones(len(X), bool), np.zeros(len(Y), bool)])
    m_per_county = {c: int(is_X[cc == c].sum()) for c in np.unique(cc)}
    grid = np.linspace(0.0, float(np.quantile(pdist(Z), q)), n_steps)
    D_obs = topo_D(X, Y, grid, max_dim)[0]
    child = np.random.SeedSequence(seed).spawn(n_perm)
    D_null = Parallel(n_jobs=n_jobs)(delayed(_perm_one_D_restr)(s, Z, cc, m_per_county, grid, max_dim) for s in child)
    D_null = np.asarray(D_null)
    pval = (np.sum(D_null >= D_obs) + 1) / (n_perm + 1)
    return D_obs, pval, D_null

# ---------------------------------------------------------------------
# 3) Run both permutation schemes for Q1_low vs Q4_high
# ---------------------------------------------------------------------
a, b = "Q1_low", "Q4_high"
X, Y = clouds[a], clouds[b]
cX = dat.loc[dat["density_q"] == a, "county_fips"].to_numpy()
cY = dat.loc[dat["density_q"] == b, "county_fips"].to_numpy()

print(f"running free permutation ({B} reps, {N_JOBS} jobs)...")
t0 = time.time()
D_free, p_free, null_free = topo_twosample_perm_par(X, Y, n_perm=B)
print(f"  D_obs={D_free:.4f}  p_free={p_free:.4f}  null_med={np.median(null_free):.1f}  ({time.time()-t0:.0f}s)")

print(f"running county-restricted permutation ({B} reps, {N_JOBS} jobs)...")
t0 = time.time()
D_restr, p_restr, null_restr = topo_twosample_restr(X, Y, cX, cY, n_perm=B)
print(f"  D_obs={D_restr:.4f}  p_restr={p_restr:.4f}  null_med={np.median(null_restr):.1f}  "
      f"null_max={null_restr.max():.1f}  ({time.time()-t0:.0f}s)")

assert abs(D_free - 669.8175) < 0.01 and abs(D_restr - 669.8175) < 0.01
if B == 9999:
    assert abs(np.median(null_free) - 39.06) < 2
    assert abs(np.median(null_restr) - 216.5) < 5
    assert abs(null_restr.max() - 408.4) < 15

pd.DataFrame({"iteration": range(1, B + 1), "D_null_free": null_free}).to_csv(
    TABLE_DIR / "fig_A4_null_free.csv", index=False)
pd.DataFrame({"iteration": range(1, B + 1), "D_null_restricted": null_restr}).to_csv(
    TABLE_DIR / "fig_A4_null_restricted.csv", index=False)

# ---------------------------------------------------------------------
# 4) Render
# ---------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5.5))
ax.hist(null_free, bins=50, alpha=0.55, color="#4393c3", label=f"free permutation null (median={np.median(null_free):.1f})")
ax.hist(null_restr, bins=50, alpha=0.55, color="#d6604d", label=f"county-restricted null (median={np.median(null_restr):.1f}, max={null_restr.max():.1f})")
ax.axvline(D_free, color="black", lw=2, ls="-", label=f"observed D = {D_free:.1f}")
ax.set_xlabel("D (Euler-characteristic two-sample statistic)")
ax.set_ylabel("frequency")
ax.set_title("Q1_low vs Q4_high: free vs. county-restricted permutation null\n(TopoTest, B=9,999 each)")
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(FIGURE_DIR / "fig_A4_county_restricted_null.png", dpi=300)
fig.savefig(FIGURE_DIR / "fig_A4_county_restricted_null.pdf")
plt.close(fig)
print("\nsaved fig_A4_county_restricted_null.png/.pdf and CSV backing tables")
