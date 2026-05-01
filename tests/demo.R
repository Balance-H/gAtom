#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(igraph))

# Prefer installed package workflow:
#   R CMD INSTALL .
#   Rscript tests/demo.R
if (requireNamespace("gAtom", quietly = TRUE)) {
  suppressPackageStartupMessages(library(gAtom))
} else {
  # Fallback for source-only debugging: load R wrappers + compiled shared lib if present
  repo_root <- normalizePath(file.path(getwd()), mustWork = TRUE)
  if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
    repo_root <- normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
  }
  src_dir <- file.path(repo_root, "src")

  dynlib <- Sys.glob(file.path(src_dir, "gAtom.*"))
  dynlib <- dynlib[grepl("\\.(so|dll|dylib)$", dynlib, ignore.case = TRUE)]

  if (length(dynlib) == 0L) {
    stop(
      "gAtom is not installed and no compiled shared library found in src/. ",
      "Please run: R CMD INSTALL ."
    )
  }

  dyn.load(dynlib[[1]], local = FALSE)
  source(file.path(repo_root, "R", "gAtom.R"))
}

g <- make_ring(6)
print(g)

out_min <- get_minimal_collapsible(g, c(2L, 3L))
cat("out_min:", out_min, "\n")

out_dec <- decompose_atoms(g)
str(out_dec)
