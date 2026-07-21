#!/usr/bin/env Rscript

# Compositional location and dispersion analyses in Aitchison geometry.

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(permute)
  library(vegan)
})

root <- normalizePath(Sys.getenv("ERP_REPO_ROOT", unset = "."), mustWork = TRUE)
processed <- file.path(root, "data", "processed")
table_dir <- file.path(root, "results", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

B <- as.integer(Sys.getenv("ERP_PERMUTATIONS", unset = "9999"))
SEED <- 14196142
stopifnot(B >= 1)

clr <- as.data.frame(read_parquet(file.path(processed, "ma_tract_clr_2023.parquet")))
strata <- as.data.frame(read_parquet(file.path(processed, "ma_tract_strata_2023.parquet")))
blocks <- as.data.frame(read_parquet(file.path(processed, "ma_block_points_proj_2023.parquet")))
if (!"tract" %in% names(clr)) names(clr)[1] <- "tract"

clr_cols <- grep("^clr_CNS(0[1-9]|1[0-9])$", names(clr), value = TRUE)
stopifnot(nrow(clr) == 1598, length(clr_cols) == 19,
          max(abs(rowSums(clr[, clr_cols]))) < 1e-8)

density <- blocks |> distinct(tract, density_q)
dat <- clr |> left_join(density, by = "tract") |>
  left_join(strata[, c("tract", "county_fips", "county_name")], by = "tract")
stopifnot(nrow(dat) == 1598, !anyNA(dat$density_q), !anyNA(dat$county_fips))

Y <- as.matrix(dat[, clr_cols])
rownames(Y) <- dat$tract
g <- factor(dat$density_q, levels = c("Q1_low", "Q2", "Q3", "Q4_high"))
county <- factor(dat$county_fips)
D <- vegdist(Y, method = "euclidean")
stopifnot(max(abs(D - dist(Y))) < 1e-10)

set.seed(SEED)
permanova_free <- adonis2(D ~ g, permutations = B, by = NULL)

dispersion <- betadisper(D, g)
set.seed(SEED)
permdisp_free <- permutest(dispersion, permutations = B, pairwise = TRUE)

county_control <- how(blocks = county, nperm = B)
set.seed(SEED)
permanova_county <- adonis2(D ~ g, permutations = county_control, by = NULL)
set.seed(SEED)
permdisp_county <- permutest(dispersion, permutations = county_control)

omnibus <- data.frame(
  test = c("PERMANOVA", "PERMANOVA", "PERMDISP", "PERMDISP"),
  permutation = c("free", "county_restricted", "free", "county_restricted"),
  statistic = c(permanova_free$F[1], permanova_county$F[1],
                permdisp_free$tab$F[1], permdisp_county$tab$F[1]),
  R2 = c(permanova_free$R2[1], permanova_county$R2[1], NA_real_, NA_real_),
  p_value = c(permanova_free[["Pr(>F)"]][1],
              permanova_county[["Pr(>F)"]][1],
              permdisp_free$tab[["Pr(>F)"]][1],
              permdisp_county$tab[["Pr(>F)"]][1]),
  permutations = B,
  seed = SEED
)
write.csv(omnibus, file.path(table_dir, "rq2_omnibus_results.csv"), row.names = FALSE)

distance_matrix <- as.matrix(D)
pairs <- combn(levels(g), 2)
location_rows <- lapply(seq_len(ncol(pairs)), function(i) {
  pair <- pairs[, i]
  keep <- g %in% pair
  subset_distance <- as.dist(distance_matrix[keep, keep])
  subset_group <- droplevels(g[keep])
  set.seed(SEED)
  fit <- adonis2(subset_distance ~ subset_group, permutations = B, by = NULL)
  data.frame(
    pair = paste(pair, collapse = " vs "),
    n = sum(keep),
    F = fit$F[1],
    t = sqrt(fit$F[1]),
    R2 = fit$R2[1],
    p_raw = fit[["Pr(>F)"]][1]
  )
})
pairwise_location <- do.call(rbind, location_rows)
pairwise_location$p_BH <- p.adjust(pairwise_location$p_raw, method = "BH")
pairwise_location$sig_BH <- pairwise_location$p_BH < 0.05
write.csv(pairwise_location, file.path(table_dir, "rq2_permanova_pairwise.csv"),
          row.names = FALSE)

pair_p <- permdisp_free$pairwise$permuted
pairwise_dispersion <- data.frame(pair = names(pair_p), p_raw = as.numeric(pair_p))
pairwise_dispersion$p_BH <- p.adjust(pairwise_dispersion$p_raw, method = "BH")
pairwise_dispersion$sig_BH <- pairwise_dispersion$p_BH < 0.05
write.csv(pairwise_dispersion, file.path(table_dir, "rq2_permdisp_pairwise.csv"),
          row.names = FALSE)

mean_distance <- tapply(dispersion$distances, g, mean)
write.csv(data.frame(stratum = names(mean_distance), mean_distance = as.numeric(mean_distance)),
          file.path(table_dir, "rq2_dispersion_by_stratum.csv"), row.names = FALSE)

tukey <- as.data.frame(TukeyHSD(dispersion)$group)
tukey$pair <- rownames(tukey)
rownames(tukey) <- NULL
write.csv(tukey, file.path(table_dir, "rq2_permdisp_tukey.csv"), row.names = FALSE)

county_table <- as.data.frame.matrix(table(county, g))
county_table$county_fips <- rownames(county_table)
rownames(county_table) <- NULL
county_table <- county_table[, c("county_fips", levels(g))]
write.csv(county_table, file.path(table_dir, "rq2_county_by_quartile.csv"), row.names = FALSE)

print(omnibus, row.names = FALSE)
print(pairwise_location, row.names = FALSE)
print(pairwise_dispersion, row.names = FALSE)
