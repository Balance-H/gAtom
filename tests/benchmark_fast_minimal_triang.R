library(igraph)
library(gRbase)

# Stub for static analysis
if (FALSE) fast_minimal_triang <- function(...) NULL

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required to load source package.")
}

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0L) {
  script_path <- normalizePath(
    sub("^--file=", "", file_arg[[1L]]),
    winslash = "/",
    mustWork = TRUE
  )
  project_root <- dirname(dirname(script_path))
} else {
  project_root <- getwd()
}

devtools::load_all(project_root)

edge_key_set <- function(graph) {
  edge_mat <- igraph::as_edgelist(graph, names = FALSE)

  if (length(edge_mat) == 0L) {
    return(character())
  }

  edge_mat <- matrix(as.integer(edge_mat), ncol = 2L)
  edge_mat <- cbind(
    pmin(edge_mat[, 1L], edge_mat[, 2L]),
    pmax(edge_mat[, 1L], edge_mat[, 2L])
  )

  sort(unique(apply(edge_mat, 1L, function(x) paste(x, collapse = "-"))))
}

compare_triangulation <- function(g) {
  if (is.null(igraph::V(g)$name)) {
    igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
  }

  ref_elapsed <- system.time(ref_tri <- gRbase::minimal_triang(g))[["elapsed"]]
  fast_elapsed <- system.time(fast_tri <- fast_minimal_triang(g))[["elapsed"]]

  ref_graph <- as(ref_tri, "igraph")

  ref_edges <- edge_key_set(ref_graph)
  fast_edges <- edge_key_set(fast_tri)
  original_edges <- edge_key_set(g)

  same_edges <- identical(ref_edges, fast_edges)
  fast_chordal <- isTRUE(igraph::is_chordal(fast_tri)$chordal)
  fast_supergraph <- all(original_edges %in% fast_edges)

  data.frame(
    vertices = igraph::vcount(g),
    edges = igraph::ecount(g),
    ref_elapsed = unname(ref_elapsed),
    fast_elapsed = unname(fast_elapsed),
    speedup = ifelse(fast_elapsed > 0, ref_elapsed / fast_elapsed, Inf),
    same_edges = same_edges,
    fast_chordal = fast_chordal,
    fast_supergraph = fast_supergraph,
    stringsAsFactors = FALSE
  )
}

set.seed(20260426)

cat("=== Correctness check (small graph) ===\n")
g_small <- igraph::sample_gnp(n = 500L, p = 0.08, directed = FALSE)
small_result <- compare_triangulation(g_small)
print(small_result)

small_ok <- isTRUE(small_result$same_edges) ||
  (isTRUE(small_result$fast_chordal) && isTRUE(small_result$fast_supergraph))

if (!small_ok) {
  stop("Correctness check failed: edge sets differ and fast result is not a valid chordal supergraph.")
}

cat("\n=== Performance check (large random graphs, 500 vertices) ===\n")
n_runs <- 3L
large_results <- vector("list", n_runs)

for (i in seq_len(n_runs)) {
  g_large <- igraph::sample_gnp(n = 500L, p = 0.08, directed = FALSE)
  run_result <- compare_triangulation(g_large)
  run_result$run <- i
  large_results[[i]] <- run_result
}

large_table <- do.call(rbind, large_results)
large_table <- large_table[
  ,
  c(
    "run",
    "vertices",
    "edges",
    "ref_elapsed",
    "fast_elapsed",
    "speedup",
    "same_edges",
    "fast_chordal",
    "fast_supergraph"
  )
]

print(large_table)

if (!all(large_table$fast_chordal & large_table$fast_supergraph)) {
  stop("At least one large-graph run produced an invalid triangulation.")
}

cat(
  sprintf(
    "\nMean speedup over %d runs: %.3fx\n",
    n_runs,
    mean(large_table$speedup)
  )
)