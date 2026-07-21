# Ball Mapper network analysis of 1,598 tract CLR vectors.
# Uses the original BallMapper package with epsilon 7.0 and seed 1.

suppressMessages({
  library(arrow)
  library(BallMapper)
  library(igraph)
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
strata_col <- c(Q1_low = "#2c7bb6", Q2 = "#abd9e9", Q3 = "#fdae61", Q4_high = "#d7191c")

# ---------------------------------------------------------------------
# 1) Load point cloud (CLR vectors) + density_q labels
# ---------------------------------------------------------------------
clr <- read_parquet(file.path(processed_dir, "ma_tract_clr_2023.parquet"))
clr_cols <- grep("^clr_CNS", names(clr), value = TRUE)
points   <- as.data.frame(clr[, clr_cols])
tract_id <- as.character(clr$tract)

blk <- read_parquet(file.path(processed_dir, "ma_block_points_proj_2023.parquet"))
lut <- unique(blk[, c("tract", "density_q")])
strat <- lut$density_q[match(tract_id, as.character(lut$tract))]
strat <- factor(strat, levels = strata_levels)
stopifnot(nrow(points) == 1598, ncol(points) == 19, sum(is.na(strat)) == 0)

# ---------------------------------------------------------------------
# 2) Original-package Ball Mapper run (eps* = 7.0, seed = 1)
#    Cached after the first run: the pure-R BallMapper package call takes
#    ~8-9 minutes on this machine, and caching lets cosmetic replotting
#    (colours, layout, labels) iterate without re-running the identical
#    fixed-seed computation. Delete the cache file to force a clean rerun.
# ---------------------------------------------------------------------
eps <- 7.0
cache_path <- file.path(checkpoint_dir, "ballmapper_eps7_seed1.rds")

if (file.exists(cache_path)) {
  cached <- readRDS(cache_path)
  l <- cached$l
  cat("loaded cached BallMapper() result from", cache_path, "\n")
} else {
  values <- data.frame(stratum = as.numeric(strat))
  set.seed(1)
  l <- BallMapper(points, values, eps)
  saveRDS(list(l = l), cache_path)
}

n_vertices <- nrow(l$vertices)
n_edges    <- nrow(l$edges)
cat(sprintf("eps* = %.2f  vertices = %d  edges = %d\n", eps, n_vertices, n_edges))
stopifnot(n_vertices == 158, n_edges == 3205)

g <- graph_from_data_frame(l$edges, directed = FALSE,
                            vertices = data.frame(id = seq_len(n_vertices)))
deg <- degree(g)
is_isolated <- deg == 0
n_iso <- sum(is_isolated)
cat("isolated balls (degree 0):", n_iso, "\n")
stopifnot(n_iso == 9)

# ---------------------------------------------------------------------
# 3) Louvain community structure (same call as RQ3_TDABM_package.R)
# ---------------------------------------------------------------------
set.seed(1)
comm <- cluster_louvain(g)
ball_comm <- membership(comm)
n_comm <- length(unique(ball_comm))
mod <- modularity(comm)
largest_comm <- max(table(ball_comm))
cat(sprintf("Louvain communities = %d  largest = %d  modularity = %.3f\n",
            n_comm, largest_comm, mod))
stopifnot(n_comm == 13, largest_comm == 47, abs(mod - 0.158) < 0.001)

singleton_comms <- as.integer(names(table(ball_comm))[table(ball_comm) == 1])
singleton_balls <- which(ball_comm %in% singleton_comms)
cat("singleton-community balls == isolated balls? ",
    setequal(singleton_balls, which(is_isolated)), "\n")
stopifnot(setequal(singleton_balls, which(is_isolated)))

# ---------------------------------------------------------------------
# 4) Per-ball majority density_q + tract-level community/isolated tables
# ---------------------------------------------------------------------
sizes <- l$vertices[, 2]
majority_q <- vapply(l$points_covered_by_landmarks, function(pts) {
  as.character(strat[pts][which.max(table(strat[pts]))])
}, character(1))
q_shares <- t(vapply(l$points_covered_by_landmarks, function(pts) {
  prop.table(table(factor(strat[pts], levels = strata_levels)))
}, numeric(4)))
colnames(q_shares) <- paste0(strata_levels, "_share")

node_tab <- data.frame(
  node_id = seq_len(n_vertices), n_tract = sizes, community = ball_comm,
  is_isolated = is_isolated, majority_density_q = majority_q
)
node_tab <- cbind(node_tab, q_shares)
write.csv(node_tab, file.path(table_dir, "fig_5_4_nodes.csv"), row.names = FALSE)
write.csv(l$edges, file.path(table_dir, "fig_5_4_edges.csv"), row.names = FALSE)

iso_ids <- which(is_isolated)
iso_tract_idx <- unlist(l$points_covered_by_landmarks[iso_ids])
write.csv(
  data.frame(ball_id = iso_ids, tract = tract_id[iso_tract_idx],
             density_q = as.character(strat[iso_tract_idx])),
  file.path(table_dir, "fig_5_4_isolated_balls.csv"), row.names = FALSE
)

adjusted_rand <- function(a, b) {
  tab <- table(a, b); row_n <- rowSums(tab); col_n <- colSums(tab); n <- sum(tab)
  choose2 <- function(x) sum(choose(x, 2))
  expected <- choose2(row_n) * choose2(col_n) / choose(n, 2)
  (choose2(as.vector(tab)) - expected) /
    ((choose2(row_n) + choose2(col_n)) / 2 - expected)
}
normalized_mutual_information <- function(a, b) {
  tab <- table(a, b); probabilities <- tab / sum(tab)
  row_p <- rowSums(probabilities); col_p <- colSums(probabilities)
  entropy <- function(p) -sum(p[p > 0] * log(p[p > 0]))
  mutual_information <- sum(ifelse(
    probabilities > 0,
    probabilities * log(probabilities / outer(row_p, col_p)), 0
  ))
  mutual_information / sqrt(entropy(row_p) * entropy(col_p))
}

votes <- matrix(0, nrow(points), max(ball_comm))
for (ball in seq_along(l$points_covered_by_landmarks)) {
  covered <- l$points_covered_by_landmarks[[ball]]
  votes[covered, ball_comm[ball]] <- votes[covered, ball_comm[ball]] + 1
}
tract_community <- max.col(votes, ties.method = "first")
ari <- adjusted_rand(tract_community, as.integer(strat))
nmi_value <- normalized_mutual_information(tract_community, as.integer(strat))
write.csv(data.frame(tract = tract_id, density_q = as.character(strat),
                     community = tract_community),
          file.path(table_dir, "ball_mapper_tract_communities.csv"), row.names = FALSE)
write.csv(data.frame(epsilon = eps, seed = 1, vertices = n_vertices,
                     edges = n_edges, isolated = n_iso,
                     communities = n_comm, modularity = mod,
                     ARI = ari, NMI = nmi_value),
          file.path(table_dir, "regional_agreement.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 5) Render: node fill = Louvain community, node border colour = majority
#    density_q, thick black border + square shape = isolated balls.
# ---------------------------------------------------------------------
draw_network <- function(seed_for_plotting = 1) {
  comm_palette <- colorRampPalette(c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3",
                                      "#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd",
                                      "#ccebc5","#ffed6f","#1f78b4"))(n_comm)
  V(g)$size  <- 4.2 * sqrt(sizes / max(sizes)) + 3.2
  V(g)$color <- comm_palette[ball_comm]
  V(g)$frame.color <- ifelse(is_isolated, "black", strata_col[majority_q])
  V(g)$frame.width <- ifelse(is_isolated, 3.5, 2.2)
  V(g)$shape <- ifelse(is_isolated, "square", "circle")
  V(g)$label <- NA
  E(g)$color <- adjustcolor("grey50", alpha.f = 0.10)
  E(g)$width <- 0.4

  set.seed(seed_for_plotting)
  layout_xy <- layout_with_fr(g, niter = 2000)
  plot(g, layout = layout_xy,
       main = sprintf("Ball Mapper graph of the 1,598 tract CLR vectors (eps* = %.1f)\nV=%d  edges=%d  isolated=%d  Louvain communities=%d  modularity=%.3f",
                       eps, n_vertices, n_edges, n_iso, n_comm, mod))
  legend("bottomleft", bty = "n", cex = 0.7,
         legend = c("isolated ball (degree 0)", paste0("majority ", strata_levels)),
         pch = c(22, rep(21, 4)),
         pt.bg = c("grey80", strata_col[strata_levels]),
         col = c("black", strata_col[strata_levels]), pt.cex = 1.3)
}

png(file.path(figure_dir, "fig_5_4_ball_mapper.png"), width = 2600, height = 2600, res = 300)
draw_network(seed_for_plotting = 1)
dev.off()

pdf(file.path(figure_dir, "fig_5_4_ball_mapper.pdf"), width = 8.7, height = 8.7)
draw_network(seed_for_plotting = 1)
dev.off()

cat("\nsaved fig_5_4_ball_mapper.png/.pdf and CSV backing tables\n")
