# Stubs for static analysis
if (FALSE) {
  .validate_graph <- function(...) NULL
  decompose_atoms <- function(...) NULL
}

#' Accelerated minimal triangulation via atom decomposition
#'
#' @param g An undirected igraph graph object.
#'
#' @return An igraph graph object obtained by triangulating each atom and
#'   merging their global edges.
#' @export
fast_minimal_triang <- function(g) {
  .validate_graph(g)
  if (isTRUE(igraph::is_directed(g))) {
    cli::cli_abort("{.arg g} must be an undirected graph.")
  }

  atoms_obj <- decompose_atoms(g)
  atoms_list <- atoms_obj$atoms
  if (is.null(atoms_list)) {
    atoms_list <- atoms_obj
  }

  if (!is.list(atoms_list)) {
    cli::cli_abort("`decompose_atoms(g)` must return a list of atom vertex sets.")
  }

  n_vertices <- igraph::vcount(g)
  edge_blocks <- vector("list", length(atoms_list))

  for (i in seq_along(atoms_list)) {
    atom <- as.integer(atoms_list[[i]])

    if (length(atom) == 0L) {
      next
    }
    if (anyNA(atom) || any(atom < 1L | atom > n_vertices)) {
      cli::cli_abort("Atom {.val {i}} contains invalid vertex indices.")
    }

    # induced_subgraph() reindexes vertices to 1..|atom|; sorting guarantees
    # that atom[local_id] is a stable local->global mapping.
    atom <- sort(atom)
    sub_g <- igraph::induced_subgraph(g, vids = atom)

    # gRbase::minimal_triang() may fail on unnamed igraph objects.
    igraph::V(sub_g)$name <- as.character(seq_len(igraph::vcount(sub_g)))

    sub_triangulated <- gRbase::minimal_triang(sub_g)
    sub_edges <- igraph::as_edgelist(sub_triangulated, names = FALSE)

    if (length(sub_edges) == 0L) {
      next
    }

    storage.mode(sub_edges) <- "integer"

    # Required local->global mapping:
    # atom[V(sub_triangulated)]
    global_ids <- atom[as.integer(igraph::V(sub_triangulated))]
    u <- global_ids[sub_edges[, 1L]]
    v <- global_ids[sub_edges[, 2L]]

    # Canonicalize undirected edges to (min, max)
    global_edges <- cbind(pmin(u, v), pmax(u, v))
    storage.mode(global_edges) <- "integer"

    edge_blocks[[i]] <- global_edges
  }

  merged_graph <- igraph::make_empty_graph(n = n_vertices, directed = FALSE)
  non_empty_blocks <- edge_blocks[!vapply(edge_blocks, is.null, logical(1))]

  if (length(non_empty_blocks) == 0L) {
    return(merged_graph)
  }

  merged_edges <- do.call(rbind, non_empty_blocks)
  if (!is.matrix(merged_edges)) {
    merged_edges <- matrix(merged_edges, ncol = 2L)
  }
  storage.mode(merged_edges) <- "integer"

  # Remove duplicated separator edges across atoms
  merged_edges <- unique(merged_edges)

  merged_graph <- igraph::add_edges(
    merged_graph,
    edges = as.vector(t(merged_edges))
  )

  igraph::simplify(merged_graph, remove.multiple = TRUE, remove.loops = TRUE)
}