# =====================================================================
# Figure A.5 -- Ball Mapper epsilon sensitivity
# ERP Sec 4.4 and the figure/table production checklist (Sec "Figure A.5").
#
# Uses the processed tract-level CLR vectors.
#
# Uses the VECTORISED Ball Mapper reimplementation (same as the Figure 5.5
# landmark bootstrap and RQ3_TDABM.R's own eps-grid diagnostic step),
# seed = 1, eps = 6.0 to 8.5 by 0.25 -- this is a deterministic single
# construction per eps value, not a permutation test, so it is cheap.
#
# Verification target (ERP Sec 4.4): eps* = 7.0 gives V ~ 156 (package
# run gives 158; the two differ only by RNG/landmark order per the
# paper's own note), iso_frac ~ 0.058, max_frac ~ 0.268.
# =====================================================================

suppressMessages({ library(arrow); library(igraph) })
root <- normalizePath(Sys.getenv("ERP_REPO_ROOT", unset = "."), mustWork = TRUE)
processed_dir <- file.path(root, "data", "processed")
table_dir <- file.path(root, "results", "tables")
figure_dir <- file.path(root, "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

clr <- read_parquet(file.path(processed_dir, "ma_tract_clr_2023.parquet"))
clr_cols <- grep("^clr_CNS", names(clr), value = TRUE)
X <- as.matrix(clr[, clr_cols])
D <- as.matrix(dist(X))
n <- nrow(D)

ball_mapper <- function(D, eps, seed = 1) {
  set.seed(seed)
  ord <- sample(nrow(D)); covered <- logical(nrow(D)); L <- integer(0)
  for (i in ord) if (!covered[i]) { L <- c(L, i); covered <- covered | (D[i, ] <= eps) }
  B  <- D[, L, drop = FALSE] <= eps
  Co <- crossprod(B + 0); k <- length(L)
  ut <- which(upper.tri(Co) & Co > 0, arr.ind = TRUE)
  edges <- data.frame(from = ut[, 1], to = ut[, 2])
  deg <- tabulate(c(edges$from, edges$to), nbins = k)
  list(vertices = data.frame(id = seq_len(k), size = colSums(B)), edges = edges, n_iso = sum(deg == 0))
}

eps_grid <- seq(6.0, 8.5, by = 0.25)
diag_tbl <- do.call(rbind, lapply(eps_grid, function(e) {
  bm <- ball_mapper(D, e, seed = 1)
  sizes <- bm$vertices$size
  data.frame(eps = e, vertices = nrow(bm$vertices), edges = nrow(bm$edges),
             isolated = bm$n_iso, iso_frac = round(bm$n_iso / nrow(bm$vertices), 3),
             med_size = median(sizes), max_size = max(sizes),
             max_frac = round(max(sizes) / n, 3))
}))
print(diag_tbl, row.names = FALSE)
write.csv(diag_tbl, file.path(table_dir, "fig_A5_epsilon_sensitivity.csv"), row.names = FALSE)

ref <- diag_tbl[diag_tbl$eps == 7.0, ]
cat(sprintf("\neps*=7.0 check: V=%d iso_frac=%.3f max_frac=%.3f (paper: V~156, iso_frac~0.058, max_frac~0.268)\n",
            ref$vertices, ref$iso_frac, ref$max_frac))

# draw_panels() draws the 3-panel row; oma reserves a dedicated bottom strip
# so the caption (drawn via mtext(..., outer=TRUE)) never overlaps the
# per-panel "epsilon" axis titles -- the previous version used a negative
# line offset with no outer margin reserved, which drew the caption on top
# of the middle panel's own x-axis label.
draw_panels <- function() {
  par(mfrow = c(1, 3), mar = c(4.2, 4.4, 2.5, 1), oma = c(3, 0, 0, 0), mgp = c(2.5, 0.7, 0))
  plot(diag_tbl$eps, diag_tbl$vertices, type = "b", pch = 16, xlab = "epsilon", ylab = "vertices (V)",
       main = "Graph size"); abline(v = 7.0, lty = 2, col = "firebrick")
  plot(diag_tbl$eps, diag_tbl$iso_frac, type = "b", pch = 16, xlab = "epsilon", ylab = "isolated fraction",
       main = "Isolated-vertex fraction"); abline(v = 7.0, lty = 2, col = "firebrick")
  plot(diag_tbl$eps, diag_tbl$max_frac, type = "b", pch = 16, xlab = "epsilon", ylab = "max covered fraction",
       main = "Largest-ball fraction"); abline(v = 7.0, lty = 2, col = "firebrick")
  mtext("Ball Mapper epsilon-grid diagnostics (vectorised reimplementation, seed=1); dashed line = eps* = 7.0",
        side = 1, outer = TRUE, line = 1.3, cex = 0.75)
}

png(file.path(figure_dir, "fig_A5_epsilon_sensitivity.png"), width = 2600, height = 1650, res = 300)
draw_panels()
dev.off()

pdf(file.path(figure_dir, "fig_A5_epsilon_sensitivity.pdf"), width = 8.7, height = 5.5)
draw_panels()
dev.off()

cat("\nsaved fig_A5_epsilon_sensitivity.png/.pdf and CSV\n")
