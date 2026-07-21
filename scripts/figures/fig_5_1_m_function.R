# Figure 5.1: sector-stratum M-function curves and an agriculture envelope.

suppressMessages({
  library(arrow)
  library(dbmss)
  library(spatstat.geom)
})

root <- normalizePath(Sys.getenv("ERP_REPO_ROOT", unset = "."), mustWork = TRUE)
processed_dir <- file.path(root, "data", "processed")
table_dir <- file.path(root, "results", "tables")
figure_dir <- file.path(root, "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
envelope_simulations <- as.integer(Sys.getenv("ERP_ENVELOPE_SIMULATIONS", unset = "999"))
strata  <- c("Q1_low", "Q2", "Q3", "Q4_high")
strata_col <- c(Q1_low = "#2c7bb6", Q2 = "#abd9e9", Q3 = "#fdae61", Q4_high = "#d7191c")

# ---------------------------------------------------------------------
# 1) Load block-level weighted point pattern (same projection as RQ1)
# ---------------------------------------------------------------------
pts <- read_parquet(file.path(processed_dir, "ma_block_points_proj_2023.parquet"))
stopifnot(nrow(pts) == 119181, sum(pts$jobs) == 3203251)

W <- owin(xrange = range(pts$x_m), yrange = range(pts$y_m))
P_by_q <- lapply(strata, function(q) {
  keep <- pts$density_q == q
  wmppp(data.frame(X = pts$x_m[keep], Y = pts$y_m[keep],
                    PointType = pts$sector[keep], PointWeight = pts$jobs[keep]),
        window = W)
})
names(P_by_q) <- strata

# ---------------------------------------------------------------------
# 2) Panel A -- five principal sectors, M(r) for r = 0..5000m by 50m
# ---------------------------------------------------------------------
principal_sectors <- c(CNS05 = "Manufacturing", CNS12 = "Professional, Scientific & Technical Services",
                        CNS06 = "Wholesale Trade", CNS07 = "Retail Trade", CNS19 = "Other Services")
r_grid_main <- seq(0, 5000, by = 50)

curve_list <- list()
for (s in names(principal_sectors)) {
  for (q in strata) {
    Ms <- Mhat(P_by_q[[q]], ReferenceType = s, NeighborType = s, r = r_grid_main)
    curve_list[[length(curve_list) + 1]] <- data.frame(
      sector = s, industry = principal_sectors[[s]], density_q = q,
      r_m = r_grid_main, M_obs = Ms$M
    )
  }
}
curves_main <- do.call(rbind, curve_list)
write.csv(curves_main, file.path(table_dir, "fig_5_1_mfunction_curves.csv"), row.names = FALSE)

# spot-check against the pipeline's r*=1000m snapshot (rq1_Msq_by_stratum.csv)
chk <- curves_main[curves_main$r_m == 1000 & curves_main$sector == "CNS12", ]
cat("CNS12 M(1000):\n"); print(chk[, c("density_q", "M_obs")])
stopifnot(abs(chk$M_obs[chk$density_q == "Q1_low"]  - 1.389) < 0.005,
          abs(chk$M_obs[chk$density_q == "Q4_high"] - 0.874) < 0.005)

# ---------------------------------------------------------------------
# 3) Panel B -- Agriculture (CNS01) in Q4_high, with Monte Carlo envelope
#    Exactly reproduces Mfunction_fixed.R's random-labelling check.
# ---------------------------------------------------------------------
r_grid_agri <- seq(0, 1500, by = 50)
set.seed(20260612)
env_CNS01_Q4 <- MEnvelope(
  P_by_q[["Q4_high"]],
  ReferenceType       = "CNS01",
  NeighborType        = "CNS01",
  r                   = r_grid_agri,
  NumberOfSimulations = envelope_simulations,
  Alpha               = 0.05,
  SimulationType      = "RandomLabeling"
)
i_star <- which.min(abs(r_grid_agri - 1000))
cat(sprintf("\nAgriculture (CNS01) Q4_high @ r=1000m: obs=%.3f  lo=%.3f  hi=%.3f\n",
            env_CNS01_Q4$obs[i_star], env_CNS01_Q4$lo[i_star], env_CNS01_Q4$hi[i_star]))
stopifnot(abs(env_CNS01_Q4$obs[i_star] - 23.049) < 0.01)

agri_tab <- data.frame(r_m = r_grid_agri, M_obs = env_CNS01_Q4$obs,
                        envelope_lower = env_CNS01_Q4$lo, envelope_upper = env_CNS01_Q4$hi)
write.csv(agri_tab, file.path(table_dir, "fig_5_1_agriculture_envelope.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 4) Render: two-panel figure (300dpi PNG + vector PDF)
# ---------------------------------------------------------------------
draw_fig <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.6, 1))
  par(mar = c(4.2, 4.4, 3, 1), mgp = c(2.5, 0.7, 0))

  # Panel A: 5 principal sectors, small multiples, shared y-axis
  par(mfrow = c(1, 1))
  layout(matrix(1:6, nrow = 2, byrow = TRUE))
  ylim_main <- range(curves_main$M_obs, 1)
  for (s in names(principal_sectors)) {
    sub <- curves_main[curves_main$sector == s, ]
    plot(NA, xlim = c(0, 5000), ylim = ylim_main, xlab = "r (m)", ylab = "M(r)",
         main = sprintf("%s (%s)", s, principal_sectors[[s]]), cex.main = 0.85)
    abline(h = 1, lty = 2, col = "grey50")
    abline(v = 1000, lty = 3, col = "grey50")
    for (q in strata) {
      d <- sub[sub$density_q == q, ]
      lines(d$r_m, d$M_obs, col = strata_col[q], lwd = 1.8)
    }
  }
  plot.new()
  legend("center", legend = strata, col = strata_col, lwd = 2, bty = "n",
         title = "Density stratum", cex = 0.9)

  # Panel B (separate device call below) -- Agriculture with envelope
}

draw_agri_panel <- function() {
  par(mar = c(4.2, 4.4, 3, 1.2), mgp = c(2.5, 0.7, 0))
  ylim_a <- range(c(env_CNS01_Q4$obs, env_CNS01_Q4$lo, env_CNS01_Q4$hi))
  plot(NA, xlim = range(r_grid_agri), ylim = ylim_a, xlab = "r (m)", ylab = "M(r)",
       main = "CNS01 Agriculture in Q4_high\n(observed vs. random-labelling envelope)", cex.main = 0.85)
  polygon(c(r_grid_agri, rev(r_grid_agri)), c(env_CNS01_Q4$hi, rev(env_CNS01_Q4$lo)),
          col = adjustcolor("grey70", alpha.f = 0.5), border = NA)
  lines(r_grid_agri, env_CNS01_Q4$obs, col = "firebrick", lwd = 2)
  abline(h = 1, lty = 2, col = "grey50")
  abline(v = 1000, lty = 3, col = "grey50")
  legend("topleft", bty = "n", cex = 0.75,
         legend = c("observed M(r)", sprintf("%d-sim random-labelling envelope (95%%)", envelope_simulations)),
         col = c("firebrick", "grey70"), lwd = c(2, 8))
}

png(file.path(figure_dir, "fig_5_1_m_function_main.png"), width = 3000, height = 2000, res = 300)
draw_fig()
dev.off()
pdf(file.path(figure_dir, "fig_5_1_m_function_main.pdf"), width = 10, height = 6.67)
draw_fig()
dev.off()

png(file.path(figure_dir, "fig_5_1_m_function_agriculture.png"), width = 1800, height = 1500, res = 300)
draw_agri_panel()
dev.off()
pdf(file.path(figure_dir, "fig_5_1_m_function_agriculture.pdf"), width = 6, height = 5)
draw_agri_panel()
dev.off()

cat("\nsaved fig_5_1_m_function_main.png/.pdf, fig_5_1_m_function_agriculture.png/.pdf, and CSV backing tables\n")
