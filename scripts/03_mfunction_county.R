#!/usr/bin/env Rscript

# County-restricted M-function permutation test.

suppressPackageStartupMessages({
  library(arrow)
  library(dbmss)
  library(dplyr)
  library(parallel)
  library(tidyr)
})

root <- normalizePath(Sys.getenv("ERP_REPO_ROOT", unset = "."), mustWork = TRUE)
data_file <- file.path(root, "data", "processed", "ma_block_points_proj_2023.parquet")
table_dir <- file.path(root, "results", "tables")
null_dir <- file.path(root, "results", "null_distributions")
checkpoint_dir <- file.path(root, "results", "checkpoints", "mfunction_county")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(null_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

B <- as.integer(Sys.getenv("ERP_PERMUTATIONS", unset = "9999"))
NCORES <- as.integer(Sys.getenv("ERP_CORES", unset = "4"))
CHUNK_SIZE <- as.integer(Sys.getenv("ERP_CHUNK", unset = "250"))
SEED <- 20260612
stopifnot(B >= 1, NCORES >= 1, CHUNK_SIZE >= 1)

pts <- read_parquet(data_file)
stopifnot(nrow(pts) == 119181, sum(pts$jobs) == 3203251)
sectors <- sprintf("CNS%02d", 1:19)
strata_levels <- c("Q1_low", "Q2", "Q3", "Q4_high")

window <- spatstat.geom::owin(xrange = range(pts$x_m), yrange = range(pts$y_m))
pattern <- wmppp(data.frame(
  X = pts$x_m, Y = pts$y_m, PointType = pts$sector, PointWeight = pts$jobs
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
analysis_family <- sectors[apply(counts, 1, min) >= 30]
stopifnot(length(analysis_family) == 18, !"CNS02" %in% analysis_family)

population_variance <- function(x) mean((x - mean(x))^2)
dense_grid <- seq(0, 1500, by = 50)
star_index <- which.min(abs(dense_grid - 1000))
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

tract_q <- pts |> distinct(tract, density_q)
tract_q$county <- substr(as.character(tract_q$tract), 1, 5)
stopifnot(nrow(tract_q) == 1598)
county_levels <- tract_q |> group_by(county) |> summarise(n = n_distinct(density_q))
stopifnot(nrow(county_levels) == 14, all(county_levels$n >= 2))

sparse_grid <- c(0, 250, 500, 750, 1000)
sparse_star <- length(sparse_grid)
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

density_labels <- as.character(tract_q$density_q)
permute_within_county <- function() {
  ave(density_labels, tract_q$county, FUN = sample)
}

n_chunks <- ceiling(B / CHUNK_SIZE)
cluster <- makeCluster(NCORES)
on.exit(try(stopCluster(cluster), silent = TRUE), add = TRUE)
clusterEvalQ(cluster, suppressPackageStartupMessages(library(dbmss)))
clusterExport(cluster, c(
  "analysis_family", "strata_levels", "population_variance", "tract_q",
  "density_labels", "bx", "by", "bsec", "bw", "btract", "global_window",
  "sparse_grid", "sparse_star", "one_permutation", "permute_within_county"
))
clusterSetRNGStream(cluster, iseed = SEED)

for (chunk in seq_len(n_chunks)) {
  checkpoint <- file.path(checkpoint_dir, sprintf("chunk_%03d.rds", chunk))
  if (file.exists(checkpoint)) next
  current <- min(CHUNK_SIZE, B - (chunk - 1) * CHUNK_SIZE)
  values <- parSapply(cluster, seq_len(current), function(i) {
    one_permutation(permute_within_county())
  })
  stopifnot(!anyNA(values))
  saveRDS(values, checkpoint)
  cat(sprintf("completed chunk %d/%d\n", chunk, n_chunks))
}
stopCluster(cluster)

checkpoint_files <- file.path(checkpoint_dir, sprintf("chunk_%03d.rds", seq_len(n_chunks)))
stopifnot(all(file.exists(checkpoint_files)))
d_null <- do.call(cbind, lapply(checkpoint_files, readRDS))
stopifnot(identical(dim(d_null), c(length(analysis_family), B)), !anyNA(d_null))
saveRDS(d_null, file.path(null_dir, sprintf("rq1_perm_null_B%d_countyrestricted.rds", B)))

p_restricted <- (1 + rowSums(d_null >= d_observed[analysis_family])) / (B + 1)
p_restricted_bh <- p.adjust(p_restricted, method = "BH")
result <- data.frame(
  sector = analysis_family,
  D_obs = round(d_observed[analysis_family], 4),
  null_median_cr = round(apply(d_null, 1, median), 4),
  p_restr = round(p_restricted, 4),
  p_restr_BH = round(p_restricted_bh, 4)
)

free_file <- file.path(table_dir, "rq1_results.csv")
if (file.exists(free_file)) {
  free <- read.csv(free_file)
  result <- merge(result, free[, c("sector", "p_BH")], by = "sector", all.x = TRUE)
  names(result)[names(result) == "p_BH"] <- "p_free_BH"
  result$sig_free <- result$p_free_BH < 0.05
  result$sig_restr <- result$p_restr_BH < 0.05
  result$flipped <- result$sig_free != result$sig_restr
}
result <- result[order(result$p_restr_BH, result$sector), ]
write.csv(result, file.path(table_dir, "rq1_county_restricted_results.csv"),
          row.names = FALSE)
print(result, row.names = FALSE)
