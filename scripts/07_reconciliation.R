#!/usr/bin/env Rscript

# Compare sector-level M-function variation with CLR between-stratum contribution.

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

root <- normalizePath(Sys.getenv("ERP_REPO_ROOT", unset = "."), mustWork = TRUE)
processed <- file.path(root, "data", "processed")
table_dir <- file.path(root, "results", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

clr <- read_parquet(file.path(processed, "ma_tract_clr_2023.parquet"))
blocks <- read_parquet(file.path(processed, "ma_block_points_proj_2023.parquet"))
if ("__index_level_0__" %in% names(clr)) {
  clr <- rename(clr, tract = `__index_level_0__`)
}
density <- blocks |> distinct(tract, density_q)
dat <- clr |> inner_join(density, by = "tract")
clr_cols <- grep("^clr_CNS", names(dat), value = TRUE)
stopifnot(nrow(dat) == 1598, length(clr_cols) == 19, !anyNA(dat$density_q))

grand_mean <- colMeans(dat[clr_cols])
groups <- split(dat, dat$density_q)
between_ss <- sapply(clr_cols, function(column) {
  sum(sapply(groups, function(group) {
    nrow(group) * (mean(group[[column]]) - grand_mean[column])^2
  }))
})
contribution <- between_ss / sum(between_ss)

rq1 <- read.csv(file.path(table_dir, "rq1_results.csv"))
comparison <- data.frame(
  sector = sub("^clr_", "", clr_cols),
  compositional_contribution = as.numeric(contribution)
) |> inner_join(rq1[, c("sector", "D_obs", "p_BH", "sig_BH")], by = "sector")
stopifnot(nrow(comparison) == 18)

test <- cor.test(comparison$D_obs, comparison$compositional_contribution,
                 method = "spearman", alternative = "greater", exact = FALSE)
comparison$D_rank <- rank(-comparison$D_obs, ties.method = "average")
comparison$composition_rank <- rank(-comparison$compositional_contribution,
                                    ties.method = "average")
comparison <- comparison[order(comparison$D_rank), ]

summary <- data.frame(
  statistic = "sector_spearman",
  rho = unname(test$estimate),
  p_value_one_sided = test$p.value,
  comparison_value = 0.50,
  meets_comparison = unname(test$estimate) >= 0.50 && test$p.value < 0.05,
  sectors = nrow(comparison),
  between_stratum_ss = sum(between_ss)
)
write.csv(comparison, file.path(table_dir, "rq3_sector_reconciliation.csv"),
          row.names = FALSE)
write.csv(summary, file.path(table_dir, "rq3_sector_spearman.csv"), row.names = FALSE)
print(comparison, row.names = FALSE)
print(summary, row.names = FALSE)
