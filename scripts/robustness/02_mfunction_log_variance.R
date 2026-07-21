#!/usr/bin/env Rscript

# Compare sector rankings using Var[M] and the level-standardised Var[log M].

suppressPackageStartupMessages({
  library(arrow)
  library(dbmss)
  library(dplyr)
})

root <- normalizePath(Sys.getenv("ERP_REPO_ROOT", unset = "."), mustWork = TRUE)
processed <- file.path(root, "data", "processed")
table_dir <- file.path(root, "results", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

points <- read_parquet(file.path(processed, "ma_block_points_proj_2023.parquet"))
clr <- read_parquet(file.path(processed, "ma_tract_clr_2023.parquet"))
sectors <- sprintf("CNS%02d", 1:19)
strata_levels <- c("Q1_low", "Q2", "Q3", "Q4_high")

window <- spatstat.geom::owin(xrange = range(points$x_m), yrange = range(points$y_m))
patterns <- lapply(strata_levels, function(q) {
  keep <- points$density_q == q
  wmppp(data.frame(
    X = points$x_m[keep], Y = points$y_m[keep],
    PointType = points$sector[keep], PointWeight = points$jobs[keep]
  ), window = window)
})
names(patterns) <- strata_levels
counts <- sapply(patterns, function(pattern) {
  table(factor(marks(pattern)$PointType, levels = sectors))
})
analysis_family <- sectors[apply(counts, 1, min) >= 30]

m_values <- sapply(strata_levels, function(q) {
  vapply(analysis_family, function(sector) {
    Mhat(patterns[[q]], ReferenceType = sector, NeighborType = sector,
         r = c(0, 1000))$M[2]
  }, numeric(1))
})
stopifnot(all(m_values > 0))
population_variance <- function(x) mean((x - mean(x))^2)
d_raw <- apply(m_values, 1, population_variance)
d_log <- apply(log(m_values), 1, population_variance)

if ("__index_level_0__" %in% names(clr)) {
  clr <- rename(clr, tract = `__index_level_0__`)
}
density <- points |> distinct(tract, density_q)
dat <- clr |> inner_join(density, by = "tract")
clr_cols <- grep("^clr_CNS", names(dat), value = TRUE)
grand_mean <- colMeans(dat[clr_cols])
groups <- split(dat, dat$density_q)
between_ss <- sapply(clr_cols, function(column) {
  sum(sapply(groups, function(group) {
    nrow(group) * (mean(group[[column]]) - grand_mean[column])^2
  }))
})
contribution <- between_ss / sum(between_ss)
names(contribution) <- sub("^clr_", "", names(contribution))

result <- data.frame(
  sector = analysis_family,
  D_raw = d_raw,
  D_log = d_log,
  compositional_contribution = contribution[analysis_family]
)
raw_test <- cor.test(result$D_raw, result$compositional_contribution,
                     method = "spearman", alternative = "greater", exact = FALSE)
log_test <- cor.test(result$D_log, result$compositional_contribution,
                     method = "spearman", alternative = "greater", exact = FALSE)
summary <- data.frame(
  variant = c("variance_M", "variance_log_M"),
  rho = c(unname(raw_test$estimate), unname(log_test$estimate)),
  p_value_one_sided = c(raw_test$p.value, log_test$p.value)
)
write.csv(result, file.path(table_dir, "rq3_mfunction_level_sensitivity.csv"),
          row.names = FALSE)
write.csv(summary, file.path(table_dir, "rq3_mfunction_level_sensitivity_summary.csv"),
          row.names = FALSE)
print(summary, row.names = FALSE)
