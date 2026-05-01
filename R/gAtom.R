#' Validate igraph input
#' @keywords internal
.validate_graph <- function(graph) {
  if (!igraph::is_igraph(graph)) {
    cli::cli_abort("{.arg graph} must be an igraph object.")
  }
  invisible(TRUE)
}

#' Prepare graph inputs for native calls
#' @keywords internal
.prepare_graph_inputs <- function(graph) {
  edge_mat <- igraph::as_edgelist(graph, names = FALSE)
  if (length(edge_mat) == 0L) {
    edge_mat <- matrix(integer(), ncol = 2)
  }

  storage.mode(edge_mat) <- "integer"

  # Convert from R's 1-based to C's 0-based indexing
  if (length(edge_mat) > 0L) {
    edge_mat <- edge_mat - 1L
  }

  list(
    edge_mat = edge_mat,
    directed = isTRUE(igraph::is_directed(graph)),
    n_vertices = as.integer(igraph::vcount(graph))
  )
}

#' Validate vertex IDs (1-based)
#' @keywords internal
.validate_vertex_ids <- function(vertices, n_vertices, arg_name) {
  vertices <- as.integer(vertices)

  if (length(vertices) == 0L) {
    cli::cli_abort("{.arg {arg_name}} must not be empty.")
  }
  if (anyNA(vertices)) {
    cli::cli_abort("{.arg {arg_name}} contains NA values.")
  }
  if (any(vertices < 1L | vertices > n_vertices)) {
    cli::cli_abort("{.arg {arg_name}} must be within [1, {.val {n_vertices}}].")
  }

  unique(vertices)
}

#' Compute minimal collapsible set from seed vertices
#'
#' @param graph An igraph graph object.
#' @param r_nodes Integer vector of seed vertices (1-based indices).
#'
#' @return Integer vector of vertices in the minimal collapsible set (1-based indices).
#' @export
get_minimal_collapsible <- function(graph, r_nodes) {
  .validate_graph(graph)
  graph_inputs <- .prepare_graph_inputs(graph)
  r_nodes <- .validate_vertex_ids(r_nodes, graph_inputs$n_vertices, "r_nodes")
  # Convert from R's 1-based to C's 0-based indexing
  r_nodes <- r_nodes - 1L

  .Call(
    "_gAtom_get_minimal_collapsible_impl",
    graph_inputs$edge_mat,
    graph_inputs$directed,
    graph_inputs$n_vertices,
    r_nodes
  )
}

#' Decompose graph into atoms and separators
#'
#' @param graph An igraph graph object.
#'
#' @return A named list with three fields:
#'   \describe{
#'     \item{atoms}{List of integer vectors, each representing an atom (1-based indices).}
#'     \item{separators}{List of integer vectors, each representing a separator (1-based indices).}
#'     \item{tree}{Integer matrix with two columns representing edges of the atom tree.}
#'   }
#' @export
decompose_atoms <- function(graph) {
  .validate_graph(graph)
  graph_inputs <- .prepare_graph_inputs(graph)

  result <- .Call(
    "_gAtom_decompose_atoms_impl",
    graph_inputs$edge_mat,
    graph_inputs$directed,
    graph_inputs$n_vertices
  )

  # Convert tree edge list to proper matrix format if present
  if (!is.null(result$tree) && length(result$tree) > 0) {
    result$tree <- matrix(result$tree, ncol = 2, byrow = TRUE)
    colnames(result$tree) <- c("from", "to")
  } else {
    result$tree <- matrix(integer(), ncol = 2)
  }

  result
}
