#!/usr/bin/env Rscript

# Refit PERMDISP after removing isolated Ball Mapper tracts.

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(vegan)
})

root <- normalizePath(Sys.getenv("ERP_REPO_ROOT", unset = "."), mustWork = TRUE)
processed <- file.path(root, "data", "processed")
table_dir <- file.path(root, "results", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

B <- as.integer(Sys.getenv("ERP_PERMUTATIONS", unset = "9999"))
SEED <- 14196142

clr <- as.data.frame(read_parquet(file.path(processed, "ma_tract_clr_2023.parquet")))
blocks <- as.data.frame(read_parquet(file.path(processed, "ma_block_points_proj_2023.parquet")))
if (!"tract" %in% names(clr)) names(clr)[1] <- "tract"
clr_cols <- grep("^clr_CNS", names(clr), value = TRUE)
density <- blocks |> distinct(tract, density_q)
dat <- clr |> left_join(density, by = "tract")
stopifnot(nrow(dat) == 1598, length(clr_cols) == 19, !anyNA(dat$density_q))

tract_id <- sprintf("%011s", trimws(as.character(dat$tract)))
tract_id <- gsub(" ", "0", tract_id)
group <- factor(dat$density_q, levels = c("Q1_low", "Q2", "Q3", "Q4_high"))
Y <- as.matrix(dat[, clr_cols])
D <- vegdist(Y, method = "euclidean")
distance_matrix <- as.matrix(D)

isolated_nine <- c(
  "25025981800", "25025981300", "25025081001", "25025981900",
  "25025981502", "25027732001", "25013812104", "25005651900",
  "25025080401"
)
special_use_four <- c("25025981800", "25025981300", "25025981900", "25025981502")
stopifnot(all(isolated_nine %in% tract_id))

run_scenario <- function(drop_ids, scenario) {
  keep <- !tract_id %in% drop_ids
  g_sub <- droplevels(group[keep])
  d_sub <- as.dist(distance_matrix[keep, keep])
  stopifnot(max(abs(d_sub - vegdist(Y[keep, ], method = "euclidean"))) < 1e-10)
  fit <- betadisper(d_sub, g_sub)
  set.seed(SEED)
  test <- permutest(fit, permutations = B, pairwise = TRUE)
  spread <- tapply(fit$distances, g_sub, mean)
  pair_p <- test$pairwise$permuted
  pairwise <- data.frame(
    scenario = scenario,
    pair = names(pair_p),
    p_raw = as.numeric(pair_p)
  )
  pairwise$p_BH <- p.adjust(pairwise$p_raw, method = "BH")
  pairwise$sig_BH <- pairwise$p_BH < 0.05
  list(
    spread = data.frame(scenario = scenario, stratum = names(spread),
                        mean_distance = as.numeric(spread)),
    pairwise = pairwise,
    omnibus = data.frame(scenario = scenario, n = sum(keep),
                         F = test$tab$F[1],
                         p_value = test$tab[["Pr(>F)"]][1])
  )
}

scenarios <- list(
  baseline = run_scenario(character(0), "baseline"),
  remove_isolated_nine = run_scenario(isolated_nine, "remove_isolated_nine"),
  remove_special_use_four = run_scenario(special_use_four, "remove_special_use_four")
)

spread <- do.call(rbind, lapply(scenarios, `[[`, "spread"))
pairwise <- do.call(rbind, lapply(scenarios, `[[`, "pairwise"))
omnibus <- do.call(rbind, lapply(scenarios, `[[`, "omnibus"))
baseline <- spread$mean_distance[spread$scenario == "baseline"]
names(baseline) <- spread$stratum[spread$scenario == "baseline"]
spread$delta_from_baseline <- spread$mean_distance - baseline[spread$stratum]

write.csv(spread, file.path(table_dir, "rq2_permdisp_leaveout_spreads.csv"),
          row.names = FALSE)
write.csv(pairwise, file.path(table_dir, "rq2_permdisp_leaveout_pairwise.csv"),
          row.names = FALSE)
write.csv(omnibus, file.path(table_dir, "rq2_permdisp_leaveout_omnibus.csv"),
          row.names = FALSE)
print(spread, row.names = FALSE)
print(omnibus, row.names = FALSE)
