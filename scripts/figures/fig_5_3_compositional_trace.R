# Figure 5.3: PERMDISP distribution and TopoTest heatmap.

suppressMessages({
  library(arrow)
  library(dplyr)
  library(vegan)
})

root <- normalizePath(Sys.getenv("ERP_REPO_ROOT", unset = "."), mustWork = TRUE)
processed_dir <- file.path(root, "data", "processed")
table_dir <- file.path(root, "results", "tables")
figure_dir <- file.path(root, "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
B <- as.integer(Sys.getenv("ERP_PERMUTATIONS", unset = "9999"))
SEED <- 14196142
strata_levels <- c("Q1_low", "Q2", "Q3", "Q4_high")
strata_col <- c(Q1_low = "#2c7bb6", Q2 = "#abd9e9", Q3 = "#fdae61", Q4_high = "#d7191c")

# ---------------------------------------------------------------------
# 1) Load CLR + density_q (mirrors PERMANOVA_clean.R Steps 1-4)
# ---------------------------------------------------------------------
clr <- as.data.frame(read_parquet(file.path(processed_dir, "ma_tract_clr_2023.parquet")))
if (!"tract" %in% names(clr)) names(clr)[1] <- "tract"
clr_cols <- grep("^clr_CNS(0[1-9]|1[0-9])$", names(clr), value = TRUE)
stopifnot(length(clr_cols) == 19)

blk <- as.data.frame(read_parquet(file.path(processed_dir, "ma_block_points_proj_2023.parquet")))
dq  <- blk %>% distinct(tract, density_q)

dat <- clr %>% left_join(dq, by = "tract")
stopifnot(nrow(dat) == 1598, !anyNA(dat$density_q))

Y <- as.matrix(dat[, clr_cols]); rownames(Y) <- dat$tract
g <- factor(dat$density_q, levels = strata_levels)

# ---------------------------------------------------------------------
# 2) Aitchison distance (= CLR-Euclidean) -> PERMDISP
# ---------------------------------------------------------------------
D  <- vegdist(Y, method = "euclidean")
bd <- betadisper(D, g)

mean_dist <- tapply(bd$distances, g, mean)
cat("mean distance to centroid:\n"); print(round(mean_dist, 4))
stopifnot(all(abs(mean_dist - c(5.6143, 5.8698, 5.9693, 6.7968)) < 0.001))

set.seed(SEED)
bd_test <- permutest(bd, permutations = B)
cat("\nPERMDISP omnibus:\n"); print(bd_test)
Fobs <- bd_test$tab$F[1]; pobs <- bd_test$tab[["Pr(>F)"]][1]
stopifnot(abs(Fobs - 55.764) < 0.01)
if (B == 9999) stopifnot(pobs <= 0.0001 + 1e-9)

tract_tab <- data.frame(tract = dat$tract, density_q = as.character(g),
                         distance_to_centroid = bd$distances)
write.csv(tract_tab, file.path(table_dir, "fig_5_3_distance_to_centroid.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 3) TopoTest pairwise shape statistics (fixed input, see header note)
#    Source: RQ2_TopoTest_results.md, Sec 2 (free-permutation, B=9,999, BH)
# ---------------------------------------------------------------------
topotest_pairs <- data.frame(
  pair   = c("Q1_low vs Q2", "Q1_low vs Q3", "Q1_low vs Q4_high",
             "Q2 vs Q3", "Q2 vs Q4_high", "Q3 vs Q4_high"),
  D      = c(85.749, 167.271, 669.817, 53.509, 398.495, 293.032),
  p_BH   = c(0.1150, 0.0012, 0.0002, 0.2827, 0.0002, 0.0002)
)
write.csv(topotest_pairs, file.path(table_dir, "fig_5_3_topotest_pairs.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 4) Render: two-panel figure (300dpi PNG + vector PDF)
# ---------------------------------------------------------------------
draw_5_3a <- function() {
  par(mar = c(4.2, 4.6, 3, 1), mgp = c(2.6, 0.7, 0))
  set.seed(SEED)  # jitter reproducibility only; no inferential randomness
  bx <- boxplot(distance_to_centroid ~ density_q, data = tract_tab,
                col = adjustcolor(strata_col[strata_levels], alpha.f = 0.35),
                border = strata_col[strata_levels], outline = FALSE, lwd = 1.6,
                xlab = "Employment-density stratum",
                ylab = "Aitchison distance to stratum centroid",
                main = "Compositional dispersion by density stratum")
  for (i in seq_along(strata_levels)) {
    d <- tract_tab$distance_to_centroid[tract_tab$density_q == strata_levels[i]]
    points(jitter(rep(i, length(d)), amount = 0.28), d,
           pch = 16, cex = 0.35, col = adjustcolor(strata_col[strata_levels[i]], alpha.f = 0.35))
  }
  lines(seq_along(strata_levels), mean_dist[strata_levels], type = "b",
        pch = 18, cex = 1.4, lwd = 2, col = "black")
  mtext(sprintf("PERMDISP omnibus: F(3,1594) = %.3f, p = %.4f  (B = %s permutations)",
                Fobs, pobs, format(B, big.mark = ",")), side = 3, line = 0.3, cex = 0.8)
}

draw_5_3b <- function() {
  par(mar = c(1, 5.5, 3, 5.5))
  n <- length(strata_levels)
  Dmat <- matrix(NA_real_, n, n, dimnames = list(strata_levels, strata_levels))
  Pmat <- matrix(NA_real_, n, n, dimnames = list(strata_levels, strata_levels))
  for (k in seq_len(nrow(topotest_pairs))) {
    ab <- strsplit(topotest_pairs$pair[k], " vs ")[[1]]
    Dmat[ab[1], ab[2]] <- Dmat[ab[2], ab[1]] <- topotest_pairs$D[k]
    Pmat[ab[1], ab[2]] <- Pmat[ab[2], ab[1]] <- topotest_pairs$p_BH[k]
  }
  col_ramp <- colorRampPalette(c("white", "#fee08b", "#d73027"))(100)
  rng <- range(Dmat, na.rm = TRUE)
  plot(NA, xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5), axes = FALSE, xlab = "", ylab = "",
       main = "TopoTest pairwise shape difference (D, Euler characteristic)")
  for (i in 1:n) for (j in 1:n) if (i != j) {
    v <- Dmat[i, j]
    colr <- if (is.na(v)) "grey95" else col_ramp[max(1, ceiling(100 * (v - rng[1]) / (rng[2] - rng[1])))]
    involves_q4 <- (strata_levels[i] == "Q4_high" || strata_levels[j] == "Q4_high")
    rect(j - 0.5, n - i + 0.5, j + 0.5, n - i + 1.5, col = colr,
         border = if (involves_q4) "black" else "grey60",
         lwd = if (involves_q4) 2.4 else 1)
    if (!is.na(v)) {
      sig <- Pmat[i, j] < 0.05
      lab <- sprintf("D=%.1f\np_BH=%.4f%s", v, Pmat[i, j], if (!sig) "\n(n.s.)" else "")
      text(j, n - i + 1, lab, cex = 0.68, font = if (sig) 2 else 1)
    }
  }
  axis(3, at = 1:n, labels = strata_levels, tick = FALSE, cex.axis = 0.85, line = -0.5)
  axis(2, at = n:1, labels = strata_levels, tick = FALSE, las = 1, cex.axis = 0.85, line = -0.5)
}

png(file.path(figure_dir, "fig_5_3_compositional_trace.png"), width = 3400, height = 1600, res = 300)
layout(matrix(1:2, nrow = 1), widths = c(1, 1.15))
draw_5_3a(); draw_5_3b()
dev.off()

pdf(file.path(figure_dir, "fig_5_3_compositional_trace.pdf"), width = 11.3, height = 5.3)
layout(matrix(1:2, nrow = 1), widths = c(1, 1.15))
draw_5_3a(); draw_5_3b()
dev.off()

cat("\nsaved fig_5_3_compositional_trace.png/.pdf and CSV backing tables\n")
