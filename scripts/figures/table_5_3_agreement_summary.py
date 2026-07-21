"""
Table 5.3 -- Two-level agreement judgement summary
ERP Sec 5.3 ("The two measurement families fail to align") and the
figure/table production checklist (Sec "Table 5.3").

No new data is read here: every value is a verified pipeline output
already reproduced and asserted in this project's other scripts --
rho/p from fig_5_2_sector_agreement.py, ARI/NMI from
fig_5_4_ball_mapper.R / fig_5_5_bootstrap_agreement.R. This script only
assembles them into the checklist's summary table; it performs no
statistical computation of its own.
"""
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def render_table(df, title, out_png, out_pdf, fontsize=10, header_wrap=16):
    import textwrap
    wrapped_cols = ["\n".join(textwrap.wrap(str(c), header_wrap)) for c in df.columns]
    n_rows, n_cols = df.shape
    fig_w = max(9, 1.55 * n_cols)
    fig_h = 0.55 * (n_rows + 1) + 0.9
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    tbl = ax.table(cellText=df.values, colLabels=wrapped_cols, cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)
    tbl.auto_set_column_width(col=list(range(n_cols)))
    tbl.scale(1, 2.0)
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

table = pd.DataFrame([
    {"Level": "Sector", "Metric": "Spearman rho", "Observed": 0.150,
     "Criterion": 0.50, "Statistical reference": "one-tailed p = 0.276",
     "Decision": "Criterion not met"},
    {"Level": "Region", "Metric": "ARI", "Observed": 0.027,
     "Criterion": 0.30, "Statistical reference": "county-restricted p = 0.0005",
     "Decision": "Criterion not met"},
    {"Level": "Region", "Metric": "NMI", "Observed": 0.037,
     "Criterion": "Supplementary", "Statistical reference": "county-restricted p = 0.0005",
     "Decision": "Very low overlap"},
])
table.to_csv(TABLE_DIR / "table_5_3_agreement_summary.csv", index=False)
print(table.to_string(index=False))
print(f"\nwrote table_5_3_agreement_summary.csv ({len(table)} rows)")

fmt = table.copy()
fmt["Observed"] = fmt["Observed"].map(lambda v: f"{v:.3f}")
render_table(fmt, "Table 5.3. Two-level agreement judgement summary",
             FIGURE_DIR / "table_5_3_agreement_summary.png",
             FIGURE_DIR / "table_5_3_agreement_summary.pdf")
print("wrote table_5_3_agreement_summary.png/.pdf")
