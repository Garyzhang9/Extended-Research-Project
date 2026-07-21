#!/usr/bin/env Rscript

required <- c(
  arrow = "20.0.0.2",
  BallMapper = "0.2.0",
  dbmss = "2.11.0",
  dplyr = "1.1.4",
  ggplot2 = "3.5.2",
  igraph = "2.1.4",
  permute = "0.9.8",
  spatstat.geom = "3.5.0",
  tidyr = "1.3.1",
  vegan = "2.7.1"
)

missing <- names(required)[!vapply(names(required), requireNamespace,
                                   logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "))
}

installed <- vapply(names(required), function(pkg) {
  as.character(packageVersion(pkg))
}, character(1))

report <- data.frame(package = names(required), expected = unname(required),
                     installed = unname(installed),
                     exact = unname(installed) == unname(required))
print(report, row.names = FALSE)
cat("R version:", as.character(getRversion()), "\n")

if (!all(report$exact)) {
  warning("Package versions differ from the reference environment.")
}
