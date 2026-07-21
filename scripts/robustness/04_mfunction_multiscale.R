# =====================================================================
# Figure A.1 -- Multi-scale M-function sensitivity
# ERP Sec 4.1 / 5.1 and the figure/table production checklist (Sec "Figure A.1").
#
# Uses the processed block points and the archived sector result table.
#
# IMPORTANT SCOPE NOTE: this is a DESCRIPTIVE robustness check only. The
# paper's actual inferential result is Ds at r*=1000m with a B=9,999
# permutation test (RQ1_Mfunction.R, ~5 hours on 4 cores). Recomputing
# that full permutation test at 500m/2000m/5000m as well would cost
# several more multi-hour runs, which is out of scope for a descriptive
# appendix figure. Instead, this script computes the same deterministic
# Ds(r) = Var_q[M_s,q(r)] formula (no permutation, no p-value) at each
# radius and asks only whether the SECTOR RANKING is stable across scale
# -- it does not claim significance at any radius but r*=1000m.
#
# Verification target: this script's own Ds(r=1000) ranking must
# reproduce rq1_results.csv's D_rank exactly (18 sectors, CNS02 excluded).
# =====================================================================

suppressMessages({ library(arrow); library(dbmss) })
root <- normalizePath(Sys.getenv("ERP_REPO_ROOT", unset = "."), mustWork = TRUE)
processed_dir <- file.path(root, "data", "processed")
table_dir <- file.path(root, "results", "tables")
figure_dir <- file.path(root, "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
strata_levels <- c("Q1_low", "Q2", "Q3", "Q4_high")

pts <- read_parquet(file.path(processed_dir, "ma_block_points_proj_2023.parquet"))
stopifnot(nrow(pts) == 119181, sum(pts$jobs) == 3203251)

rq1 <- read.csv(file.path(table_dir, "rq1_results.csv"))
Family1 <- rq1$sector  # 18 sectors, CNS02 already excluded upstream

W <- spatstat.geom::owin(xrange = range(pts$x_m), yrange = range(pts$y_m))
P_by_q <- lapply(strata_levels, function(q) {
  keep <- pts$density_q == q
  wmppp(data.frame(X = pts$x_m[keep], Y = pts$y_m[keep],
                    PointType = pts$sector[keep], PointWeight = pts$jobs[keep]),
        window = W)
})
names(P_by_q) <- strata_levels

varpop <- function(v) mean((v - mean(v))^2)
radii <- c(500, 1000, 2000, 5000)

Ds_by_r <- sapply(radii, function(r_target) {
  r_vec <- sort(unique(c(0, r_target)))
  i_target <- which(r_vec == r_target)
  vapply(Family1, function(s) {
    Ms <- vapply(strata_levels, function(q) {
      Mhat(P_by_q[[q]], ReferenceType = s, NeighborType = s, r = r_vec)$M[i_target]
    }, numeric(1))
    varpop(Ms)
  }, numeric(1))
})
colnames(Ds_by_r) <- paste0("r", radii, "m")
rownames(Ds_by_r) <- Family1

# consistency check: r=1000m column must reproduce rq1_results.csv ranking
# round to the pipeline's own stored precision (4dp) before ranking, since
# rq1_results.csv itself rounds D_obs to 4dp and ties at that precision
# (e.g. CNS18/CNS11/CNS15 all 0.0079) would otherwise break differently
# under this script's full-precision recomputation
rank_1000 <- rank(-round(Ds_by_r[, "r1000m"], 4), ties.method = "min")
rq1_rank  <- rank(-rq1$D_obs, ties.method = "min"); names(rq1_rank) <- rq1$sector
cat("r=1000m descriptive Ds vs rq1_results.csv D_obs (should match to rounding):\n")
print(round(cbind(Ds_here = Ds_by_r[, "r1000m"], D_obs_pipeline = rq1$D_obs[match(Family1, rq1$sector)]), 4))
stopifnot(all(rank_1000 == rq1_rank[Family1]))
cat("rank(r=1000m) matches rq1_results.csv D_rank exactly: TRUE\n\n")

rank_by_r <- apply(-Ds_by_r, 2, rank, ties.method = "min")
write.csv(data.frame(sector = Family1, Ds_by_r, rank_by_r), file.path(table_dir, "fig_A1_multiscale_sensitivity.csv"), row.names = FALSE)

rho_vs_1000 <- sapply(colnames(rank_by_r), function(cn) cor(rank_by_r[, "r1000m"], rank_by_r[, cn], method = "spearman"))
cat("Spearman rho of sector ranking at each radius vs. the r*=1000m ranking:\n")
print(round(rho_vs_1000, 3))

# ---------------------------------------------------------------------
# Render: rank heatmap (sectors x radii), 1 = strongest concentration
# ---------------------------------------------------------------------
ord <- order(rank_by_r[, "r1000m"])
draw_fig <- function() {
  par(mar = c(4, 7, 3, 5))
  n_s <- length(Family1); n_r <- length(radii)
  image(1:n_r, 1:n_s, t(rank_by_r[ord, ]), col = rev(heat.colors(18)),
        axes = FALSE, xlab = "", ylab = "", main = "Sector rank by radius (1 = strongest concentration, Ds descriptive only)")
  axis(1, at = 1:n_r, labels = colnames(rank_by_r)); axis(2, at = 1:n_s, labels = Family1[ord], las = 1, cex.axis = 0.8)
  for (i in 1:n_s) for (j in 1:n_r) text(j, i, rank_by_r[ord, ][i, j], cex = 0.7)
  mtext(paste0("Spearman rho vs r*=1000m ranking: ",
               paste(sprintf("%s=%.2f", colnames(rank_by_r), rho_vs_1000), collapse = "  ")),
        side = 1, line = 2.8, cex = 0.75)
}
png(file.path(figure_dir, "fig_A1_multiscale_sensitivity.png"), width = 2200, height = 2200, res = 300)
draw_fig(); dev.off()
pdf(file.path(figure_dir, "fig_A1_multiscale_sensitivity.pdf"), width = 7.3, height = 7.3)
draw_fig(); dev.off()

cat("\nsaved fig_A1_multiscale_sensitivity.png/.pdf and CSV\n")
