#!/usr/bin/env Rscript

# Sector-level M-function profiles and free-label permutation test.

suppressPackageStartupMessages({
  library(arrow)
  library(dbmss)
  library(dplyr)
  library(parallel)
})

root <- normalizePath(Sys.getenv("ERP_REPO_ROOT", unset = "."), mustWork = TRUE)
data_file <- file.path(root, "data", "processed", "ma_block_points_proj_2023.parquet")
table_dir <- file.path(root, "results", "tables")
null_dir <- file.path(root, "results", "null_distributions")
checkpoint_dir <- file.path(root, "results", "checkpoints", "mfunction_free")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(null_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

R_STAR <- 1000
MIN_POINTS <- 30
B <- as.integer(Sys.getenv("ERP_PERMUTATIONS", unset = "9999"))
SEED <- 20260612
NCORES <- as.integer(Sys.getenv("ERP_CORES", unset = "4"))
BATCH <- as.integer(Sys.getenv("ERP_BATCH", unset = "1000"))
stopifnot(B >= 1, NCORES >= 1, BATCH >= 1)

sectors <- sprintf("CNS%02d", 1:19)
strata_levels <- c("Q1_low", "Q2", "Q3", "Q4_high")

pts <- read_parquet(data_file)
stopifnot(nrow(pts) == 119181, length(unique(pts$sector)) == 19,
          sum(pts$jobs) == 3203251)

window <- spatstat.geom::owin(xrange = range(pts$x_m), yrange = range(pts$y_m))
pattern <- wmppp(data.frame(
  X = pts$x_m,
  Y = pts$y_m,
  PointType = pts$sector,
  PointWeight = pts$jobs
), window = window)

patterns_by_stratum <- lapply(strata_levels, function(q) {
  keep <- pts$density_q == q & pts$sector %in% sectors
  wmppp(data.frame(
    X = pts$x_m[keep], Y = pts$y_m[keep],
    PointType = pts$sector[keep], PointWeight = pts$jobs[keep]
  ), window = Window(pattern))
})
names(patterns_by_stratum) <- strata_levels

counts <- sapply(strata_levels, function(q) {
  table(factor(marks(patterns_by_stratum[[q]])$PointType, levels = sectors))
})
analysis_family <- sectors[apply(counts, 1, min) >= MIN_POINTS]
stopifnot(length(analysis_family) == 18, !"CNS02" %in% analysis_family)

population_variance <- function(x) mean((x - mean(x))^2)
dense_grid <- seq(0, 1500, by = 50)
star_index <- which.min(abs(dense_grid - R_STAR))
m_values <- matrix(NA_real_, nrow = length(analysis_family),
                   ncol = length(strata_levels),
                   dimnames = list(analysis_family, strata_levels))
for (q in strata_levels) {
  for (sector in analysis_family) {
    curve <- Mhat(patterns_by_stratum[[q]], ReferenceType = sector,
                  NeighborType = sector, r = dense_grid)
    m_values[sector, q] <- curve$M[star_index]
  }
}
d_observed <- apply(m_values, 1, population_variance)

write.csv(
  data.frame(sector = rownames(m_values), m_values, check.names = FALSE),
  file.path(table_dir, "rq1_Msq_by_stratum.csv"), row.names = FALSE
)
write.csv(
  data.frame(sector = rownames(counts), counts, check.names = FALSE),
  file.path(table_dir, "rq1_block_counts_by_stratum.csv"), row.names = FALSE
)

sparse_grid <- c(0, 250, 500, 750, 1000)
sparse_star <- length(sparse_grid)
tract_q <- pts |> distinct(tract, density_q)
stopifnot(nrow(tract_q) == 1598)
bx <- pts$x_m
by <- pts$y_m
bsec <- pts$sector
bw <- pts$jobs
btract <- pts$tract
global_window <- Window(pattern)

one_permutation <- function(permuted_labels) {
  label_map <- setNames(permuted_labels, tract_q$tract)
  block_labels <- label_map[btract]
  result <- matrix(NA_real_, nrow = length(analysis_family),
                   ncol = length(strata_levels),
                   dimnames = list(analysis_family, strata_levels))
  for (q in strata_levels) {
    keep <- which(block_labels == q)
    p_q <- wmppp(data.frame(
      X = bx[keep], Y = by[keep], PointType = bsec[keep], PointWeight = bw[keep]
    ), window = global_window)
    marks_q <- marks(p_q)$PointType
    for (sector in analysis_family) {
      if (sum(marks_q == sector) >= 1) {
        curve <- Mhat(p_q, ReferenceType = sector, NeighborType = sector,
                      r = sparse_grid)
        result[sector, q] <- curve$M[sparse_star]
      }
    }
  }
  apply(result, 1, population_variance)
}

cluster <- makeCluster(NCORES)
on.exit(try(stopCluster(cluster), silent = TRUE), add = TRUE)
clusterEvalQ(cluster, suppressPackageStartupMessages(library(dbmss)))
clusterExport(cluster, c(
  "analysis_family", "strata_levels", "population_variance", "tract_q",
  "bx", "by", "bsec", "bw", "btract", "global_window", "sparse_grid",
  "sparse_star", "one_permutation"
))
clusterSetRNGStream(cluster, iseed = SEED)

checkpoint <- file.path(checkpoint_dir, sprintf("free_B%d.rds", B))
d_null <- if (file.exists(checkpoint)) readRDS(checkpoint) else {
  matrix(NA_real_, nrow = length(analysis_family), ncol = 0,
         dimnames = list(analysis_family, NULL))
}
completed <- ncol(d_null)
cat(sprintf("M-function free permutations: %d/%d already available\n", completed, B))

while (completed < B) {
  current <- min(BATCH, B - completed)
  batch_values <- parSapply(cluster, seq_len(current), function(i) {
    one_permutation(sample(tract_q$density_q))
  })
  d_null <- cbind(d_null, batch_values)
  completed <- ncol(d_null)
  saveRDS(d_null, checkpoint)
  cat(sprintf("completed %d/%d\n", completed, B))
}
stopCluster(cluster)
stopifnot(identical(dim(d_null), c(length(analysis_family), B)), !anyNA(d_null))

archive <- list(
  D_null = d_null,
  D_obs = d_observed,
  sectors = analysis_family,
  M_by_stratum = m_values,
  seed = SEED,
  permutations = B,
  radius_m = R_STAR,
  workers = NCORES
)
saveRDS(archive, file.path(null_dir, sprintf("rq1_perm_null_B%d.rds", B)))

p_raw <- (1 + rowSums(d_null >= d_observed)) / (B + 1)
p_bh <- p.adjust(p_raw, method = "BH")
result <- data.frame(
  sector = analysis_family,
  D_obs = round(d_observed, 4),
  p_raw = round(p_raw, 4),
  p_BH = round(p_bh, 4),
  sig_BH = p_bh < 0.05
)
result <- result[order(result$p_BH, -result$D_obs), ]
write.csv(result, file.path(table_dir, "rq1_results.csv"), row.names = FALSE)
print(result, row.names = FALSE)
