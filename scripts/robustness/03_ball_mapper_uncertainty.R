# Ball Mapper landmark bootstrap and county-restricted regional null.

suppressMessages({
  library(arrow)
  library(dplyr)
  library(igraph)
  library(BallMapper)
})

root <- normalizePath(Sys.getenv("ERP_REPO_ROOT", unset = "."), mustWork = TRUE)
processed_dir <- file.path(root, "data", "processed")
table_dir <- file.path(root, "results", "tables")
figure_dir <- file.path(root, "results", "figures")
checkpoint_dir <- file.path(root, "results", "checkpoints")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
strata_levels <- c("Q1_low", "Q2", "Q3", "Q4_high")

# ---------------------------------------------------------------------
# 1) Load point cloud + density_q labels (same as Figures 5.1/5.4)
# ---------------------------------------------------------------------
clr <- read_parquet(file.path(processed_dir, "ma_tract_clr_2023.parquet"))
clr_cols <- grep("^clr_CNS", names(clr), value = TRUE)
X        <- as.matrix(clr[, clr_cols])
points   <- as.data.frame(X)
tract_id <- as.character(clr$tract)

blk <- read_parquet(file.path(processed_dir, "ma_block_points_proj_2023.parquet"))
lut <- unique(blk[, c("tract", "density_q")])
strat <- lut$density_q[match(tract_id, as.character(lut$tract))]
strat <- factor(strat, levels = strata_levels)
y <- as.integer(strat)
stopifnot(nrow(points) == 1598, sum(is.na(strat)) == 0)

# ---------------------------------------------------------------------
# 2) ARI / NMI helpers (self-computed, matching RQ3_TDABM.R exactly)
# ---------------------------------------------------------------------
adj_rand <- function(a, b) {
  tab <- table(a, b); ni <- rowSums(tab); nj <- colSums(tab); n <- sum(tab)
  c2 <- function(x) sum(choose(x, 2))
  (c2(as.vector(tab)) - c2(ni) * c2(nj) / choose(n, 2)) /
    ((c2(ni) + c2(nj)) / 2 - c2(ni) * c2(nj) / choose(n, 2))
}
nmi <- function(a, b) {
  tab <- table(a, b); n <- sum(tab); pij <- tab / n; pi <- rowSums(pij); pj <- colSums(pij)
  H <- function(p) -sum(p[p > 0] * log(p[p > 0]))
  MI <- sum(ifelse(pij > 0, pij * log(pij / outer(pi, pj)), 0)); MI / sqrt(H(pi) * H(pj))
}

# ---------------------------------------------------------------------
# 3) Single-run reference: original-package BallMapper (reuse Figure 5.4
#    cache; same eps*=7.0, seed=1 run) -> Louvain -> majority-vote tract
#    community -> ARI/NMI vs density_q
# ---------------------------------------------------------------------
cache_path <- file.path(checkpoint_dir, "ballmapper_eps7_seed1.rds")
if (file.exists(cache_path)) {
  l <- readRDS(cache_path)$l
  cat("loaded cached original-package BallMapper() result from Figure 5.4\n")
} else {
  values <- data.frame(stratum = as.numeric(strat))
  set.seed(1)
  l <- BallMapper(points, values, 7.0)
}
stopifnot(nrow(l$vertices) == 158, nrow(l$edges) == 3205)

g_pkg <- graph_from_data_frame(l$edges, directed = FALSE,
                                vertices = data.frame(id = seq_len(nrow(l$vertices))))
set.seed(1)
ball_comm <- membership(cluster_louvain(g_pkg))

n <- nrow(points)
votes <- matrix(0, n, length(unique(ball_comm)))
for (b in seq_along(l$points_covered_by_landmarks)) {
  pts <- l$points_covered_by_landmarks[[b]]
  votes[pts, ball_comm[b]] <- votes[pts, ball_comm[b]] + 1
}
tract_comm <- max.col(votes, ties.method = "first")

ari_main <- adj_rand(tract_comm, y)
nmi_main <- nmi(tract_comm, y)
cat(sprintf("single-run (package, eps*=7.0, seed=1):  ARI = %.4f  NMI = %.4f\n", ari_main, nmi_main))
stopifnot(abs(ari_main - 0.0274) < 0.001, abs(nmi_main - 0.0366) < 0.001)

# ---------------------------------------------------------------------
# 4) Landmark bootstrap (B=2,000): VECTORISED reimplementation, an
#    independent landmark order per seed -- matching RQ3_TDABM.R exactly.
#    Deliberately NOT the original package (too slow for B=2,000 reps).
# ---------------------------------------------------------------------
D <- as.matrix(dist(X))

ball_mapper_vec <- function(D, eps, seed = 1) {
  n <- nrow(D); set.seed(seed)
  ord <- sample(n); covered <- logical(n); L <- integer(0)
  for (i in ord) if (!covered[i]) { L <- c(L, i); covered <- covered | (D[i, ] <= eps) }
  B  <- D[, L, drop = FALSE] <= eps
  Co <- crossprod(B + 0); k <- length(L)
  ut <- which(upper.tri(Co) & Co > 0, arr.ind = TRUE)
  edges <- data.frame(from = ut[, 1], to = ut[, 2])
  pcbl  <- lapply(seq_len(k), function(j) which(B[, j]))
  list(vertices = data.frame(id = seq_len(k), size = colSums(B)),
       edges = edges, landmarks = L, points_covered_by_landmarks = pcbl)
}

eval_once <- function(seed, eps = 7.0) {
  bm <- ball_mapper_vec(D, eps, seed = seed)
  g  <- graph_from_data_frame(bm$edges, directed = FALSE,
                              vertices = data.frame(id = seq_len(nrow(bm$vertices))))
  set.seed(seed)
  bc <- membership(cluster_louvain(g))
  votes <- matrix(0, n, max(bc))
  for (b in seq_along(bm$points_covered_by_landmarks)) {
    pts <- bm$points_covered_by_landmarks[[b]]
    votes[pts, bc[b]] <- votes[pts, bc[b]] + 1
  }
  tc <- max.col(votes, ties.method = "first")
  c(seed = seed, ari = adj_rand(tc, y), nmi = nmi(tc, y),
    n_vertices = nrow(bm$vertices), n_communities = length(unique(bc)))
}

B <- as.integer(Sys.getenv("ERP_BOOTSTRAPS", unset = "2000"))
t0 <- Sys.time()
boot <- t(sapply(1:B, eval_once))
cat(sprintf("landmark bootstrap: %d reps in %.1f min\n", B, as.numeric(Sys.time() - t0, units = "mins")))

boot_df <- as.data.frame(boot)
write.csv(boot_df, file.path(table_dir, "fig_5_5_bootstrap_landmark.csv"), row.names = FALSE)

ari_med <- median(boot_df$ari); ari_ci <- quantile(boot_df$ari, c(.025, .975))
nmi_med <- median(boot_df$nmi); nmi_ci <- quantile(boot_df$nmi, c(.025, .975))
cat(sprintf("ARI: median=%.4f  95%% CI=[%.4f, %.4f]\n", ari_med, ari_ci[1], ari_ci[2]))
cat(sprintf("NMI: median=%.4f  95%% CI=[%.4f, %.4f]\n", nmi_med, nmi_ci[1], nmi_ci[2]))
if (B == 2000) {
  stopifnot(abs(ari_med - 0.0249) < 0.003, abs(nmi_med - 0.0411) < 0.003)
  stopifnot(abs(ari_ci[1] - 0.0099) < 0.005, abs(ari_ci[2] - 0.0468) < 0.005)
  stopifnot(abs(nmi_ci[1] - 0.0229) < 0.005, abs(nmi_ci[2] - 0.0609) < 0.005)
}

# ---------------------------------------------------------------------
# 5) County-restricted permutation null (B=1,999): fixed package Louvain
#    partition (tract_comm), density_q shuffled within county only --
#    matching RQ3_county_restricted.R exactly.
# ---------------------------------------------------------------------
strata_tab <- read_parquet(file.path(processed_dir, "ma_tract_strata_2023.parquet"))
county <- as.factor(strata_tab$county_fips[match(tract_id, as.character(strata_tab$tract))])
stopifnot(sum(is.na(county)) == 0)

set.seed(1)
Bn <- as.integer(Sys.getenv("ERP_REGIONAL_PERMUTATIONS", unset = "1999"))
idx_by_county <- split(seq_along(y), county)
perm_within <- function() {
  yp <- y
  for (ix in idx_by_county) if (length(ix) > 1) yp[ix] <- sample(y[ix])
  yp
}
perm <- replicate(Bn, {
  yp <- perm_within()
  c(adj_rand(tract_comm, yp), nmi(tract_comm, yp))
})
null_df <- data.frame(iteration = seq_len(Bn), ARI_null = perm[1, ], NMI_null = perm[2, ])
write.csv(null_df, file.path(table_dir, "fig_5_5_county_null.csv"), row.names = FALSE)

p_ari <- (1 + sum(null_df$ARI_null >= ari_main)) / (Bn + 1)
p_nmi <- (1 + sum(null_df$NMI_null >= nmi_main)) / (Bn + 1)
cat(sprintf("county-restricted permutation:  p_ARI = %.4f  p_NMI = %.4f\n", p_ari, p_nmi))
if (Bn == 1999) stopifnot(abs(p_ari - 0.0005) < 1e-9 || p_ari <= 0.0005 + 1e-9)

# ---------------------------------------------------------------------
# 6) Render: two-panel figure (300dpi PNG + vector PDF)
# ---------------------------------------------------------------------
draw_panel <- function(boot_x, null_x, med, ci, single_run, thresh, lab, col_fill) {
  h <- hist(boot_x, breaks = 40, plot = FALSE)
  hn <- hist(null_x, breaks = 40, plot = FALSE)
  # Zoom the axis to the data actually being compared (bootstrap + null);
  # the ARI>=0.30 criterion sits far outside this range and is marked with
  # an off-scale arrow below rather than stretched into the axis, which
  # would otherwise squash the whole bootstrap distribution into a sliver.
  xr <- range(c(boot_x, null_x, single_run), na.rm = TRUE)
  xr <- xr + c(-1, 1) * 0.05 * diff(xr)
  plot(h, col = adjustcolor(col_fill, alpha.f = 0.75), border = "white",
       freq = FALSE, main = lab, xlab = paste0(lab, "  (bootstrap B=2,000 vs. county-null B=1,999)"),
       ylab = "density", xlim = xr)
  plot(hn, col = adjustcolor("grey40", alpha.f = 0.35), border = "white", freq = FALSE, add = TRUE)
  abline(v = ci, col = "grey30", lty = 2, lwd = 1.6)
  abline(v = med, col = "black", lty = 1, lwd = 2)
  abline(v = single_run, col = "royalblue", lty = 3, lwd = 2)
  if (!is.null(thresh)) {
    usr <- par("usr")
    arrows(x0 = xr[2] * 0.82, y0 = usr[4] * 0.9, x1 = xr[2] * 0.97, y1 = usr[4] * 0.9,
           length = 0.08, col = "firebrick", lwd = 2)
    text(xr[2] * 0.80, usr[4] * 0.9, sprintf("agreement criterion\nARI %s 0.30 (off-scale)", "≥"),
         col = "firebrick", cex = 0.62, adj = c(1, 0.5))
  }
  legend("topright", bty = "n", cex = 0.68, seg.len = 1.6, inset = c(0, 0.12),
         legend = c(sprintf("bootstrap median = %.4f", med),
                    sprintf("bootstrap 95%% CI [%.4f, %.4f]", ci[1], ci[2]),
                    sprintf("single-run estimate = %.4f", single_run),
                    "county-restricted null (B=1,999)"),
         col = c("black", "grey30", "royalblue", "grey40"),
         lty = c(1, 2, 3, NA), lwd = c(2, 1.6, 2, NA), pch = c(NA, NA, NA, 15))
}

draw_fig <- function() {
  par(mfrow = c(1, 2), mar = c(4.2, 4.4, 3, 1.2), mgp = c(2.5, 0.7, 0))
  draw_panel(boot_df$ari, null_df$ARI_null, ari_med, ari_ci, ari_main, 0.30, "ARI", "#cfe0f0")
  draw_panel(boot_df$nmi, null_df$NMI_null, nmi_med, nmi_ci, nmi_main, NULL, "NMI", "#d6ecd6")
  mtext("Landmark-seed bootstrap vs. county-restricted permutation null: two distinct uncertainty analyses (not the same distribution)",
        side = 3, outer = TRUE, line = -1.3, cex = 0.8)
}

png(file.path(figure_dir, "fig_5_5_bootstrap_agreement.png"), width = 3400, height = 1600, res = 300)
draw_fig()
dev.off()

pdf(file.path(figure_dir, "fig_5_5_bootstrap_agreement.pdf"), width = 11.3, height = 5.3)
draw_fig()
dev.off()

cat("\nsaved fig_5_5_bootstrap_agreement.png/.pdf and CSV backing tables\n")
