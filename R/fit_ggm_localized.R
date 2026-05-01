#' Localized Gaussian graphical model fitting via atom decomposition
#'
#' Fits local concentration models on graph atoms and assembles a global
#' precision estimate using separator correction:
#' \deqn{
#' \hat{K} = \sum_{V \in \mathcal{V}} [\hat{K}_{VV}]^\Delta
#' - n \sum_{S \in \tilde{\mathcal{S}}} [(W_{SS})^{-1}]^\Delta,
#' \quad W_{SS} = n S_{SS}.
#' }
#'
#' @param data Numeric matrix/data.frame. Rows are observations, columns are
#'   variables.
#' @param graph An undirected \code{igraph} object with
#'   \code{igraph::vcount(graph) == ncol(data)}.
#' @param method Local estimator used on each atom. One of \code{"ggm"} or
#'   \code{"glasso"}.
#' @param ... Additional arguments passed to \code{ggm::fitConGraph()} or
#'   \code{glasso::glasso()}.
#'
#' @return A named list with:
#' \describe{
#'   \item{K}{Global precision estimate.}
#'   \item{Sigma}{Inverse of \code{K}.}
#'   \item{S}{Sample covariance used for fitting.}
#'   \item{n}{Effective sample size after complete-case filtering.}
#'   \item{atoms}{Atom vertex sets (1-based).}
#'   \item{separators}{Separator vertex sets (1-based).}
#'   \item{method}{Selected fitting method.}
#' }
# Stubs for static analysis
if (FALSE) {
  .validate_graph <- function(...) NULL
  decompose_atoms <- function(...) NULL
}

fit_ggm_localized <- function(data, graph, method = c("ggm", "glasso"), ...) {
  .validate_graph(graph)

  if (isTRUE(igraph::is_directed(graph))) {
    cli::cli_abort("{.arg graph} must be an undirected igraph object.")
  }

  method <- match.arg(method)

  if (!is.matrix(data) && !is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a matrix or data.frame.")
  }

  if (is.data.frame(data)) {
    is_num_col <- vapply(data, is.numeric, logical(1))
    if (!all(is_num_col)) {
      cli::cli_abort("{.arg data} must contain only numeric columns.")
    }
  }

  data <- as.matrix(data)
  if (!is.numeric(data)) {
    cli::cli_abort("{.arg data} must be numeric.")
  }
  storage.mode(data) <- "double"

  finite_or_na <- is.finite(data) | is.na(data)
  if (!all(finite_or_na)) {
    cli::cli_abort("{.arg data} contains non-finite values.")
  }

  p <- as.integer(igraph::vcount(graph))
  if (ncol(data) != p) {
    cli::cli_abort(
      "{.arg ncol(data)} must equal number of graph vertices ({.val {p}})."
    )
  }

  # Filter marginally complete data for each atom and separator
  filter_marginal_data <- function(data, indices) {
    rowSums(is.na(data[, indices, drop = FALSE])) == 0
  }

  n <- nrow(data)

  if (n < 2L) {
    cli::cli_abort("Need at least 2 complete observations to compute covariance.")
  }

  decomp <- decompose_atoms(graph)
  atoms <- decomp$atoms
  separators <- decomp$separators

  if (!is.list(atoms)) {
    cli::cli_abort("`decompose_atoms(graph)$atoms` must be a list.")
  }
  if (is.null(separators)) {
    separators <- list()
  }
  if (!is.list(separators)) {
    cli::cli_abort("`decompose_atoms(graph)$separators` must be a list.")
  }

  K <- matrix(0, nrow = p, ncol = p)

  for (i in seq_along(atoms)) {
    atom <- as.integer(atoms[[i]])
    if (length(atom) == 0L) {
      next
    }
    if (anyNA(atom) || any(atom < 1L | atom > p)) {
      cli::cli_abort("Atom {.val {i}} contains invalid vertex indices.")
    }
    atom <- sort(unique(atom))

    atom_complete <- filter_marginal_data(data, atom)
    if (sum(atom_complete) < 2L) {
      next
    }
    data_atom <- data[atom_complete, atom, drop = FALSE]
    S_sub <- stats::cov(data_atom)

    g_sub <- igraph::induced_subgraph(graph, vids = atom)
    adj_sub <- as.matrix(igraph::as_adjacency_matrix(g_sub, sparse = FALSE))
    storage.mode(adj_sub) <- "integer"

    fit_result <- .fit_local_precision(
      S_sub = S_sub,
      adj_sub = adj_sub,
      n = n,
      method = method,
      ...
    )

    K_sub <- fit_result$K
    K[atom, atom] <- K[atom, atom] + K_sub
  }

  # Separator correction:
  # n * (W_SS)^(-1), W_SS = n * S_SS  =>  solve(S_SS)
  for (i in seq_along(separators)) {
    sep <- as.integer(separators[[i]])
    if (length(sep) == 0L) {
      next
    }
    if (anyNA(sep) || any(sep < 1L | sep > p)) {
      cli::cli_abort("Separator {.val {i}} contains invalid vertex indices.")
    }
    sep <- sort(unique(sep))

    sep_complete <- filter_marginal_data(data, sep)
    if (sum(sep_complete) < 2L) {
      next
    }
    data_sep <- data[sep_complete, sep, drop = FALSE]
    S_sep <- stats::cov(data_sep)

    W_sep <- n * S_sep
    K_sep <- n * .safe_solve(W_sep, label = paste0("separator ", i))
    K[sep, sep] <- K[sep, sep] - K_sep
  }

  Sigma <- .safe_solve(K, label = "global precision")

  list(
    K = K,
    Sigma = Sigma,
    S = NULL,
    n = n,
    atoms = atoms,
    separators = separators,
    method = method
  )
}

#' Safe matrix inverse with small ridge fallback
#' @keywords internal
.safe_solve <- function(mat, label, ridge = 1e-8) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"

  inv <- tryCatch(solve(mat), error = function(e) NULL)
  if (!is.null(inv)) {
    return(inv)
  }

  p <- nrow(mat)
  inv <- tryCatch(
    solve(mat + diag(ridge, p)),
    error = function(e) NULL
  )
  if (!is.null(inv)) {
    return(inv)
  }

  cli::cli_abort("Failed to invert {.field {label}} matrix.")
}

#' Extract precision matrix from fitter output
#' @keywords internal
.extract_precision <- function(fit_obj) {
  if (is.matrix(fit_obj)) {
    return(fit_obj)
  }

  if (!is.list(fit_obj)) {
    cli::cli_abort("Unsupported local fit output type.")
  }

  for (nm in c("K", "Khat", "iSigma", "invSigma", "wi", "Omega", "precision")) {
    val <- fit_obj[[nm]]
    if (is.matrix(val)) {
      return(val)
    }
  }

  for (nm in c("Shat", "Sigma", "cov", "S")) {
    val <- fit_obj[[nm]]
    if (is.matrix(val)) {
      return(.safe_solve(val, label = paste0("fit output field ", nm)))
    }
  }

  cli::cli_abort("Cannot extract precision matrix from local fit result.")
}

#' Fit local precision matrix on one atom
#' @keywords internal
#' @return A list with \code{K} (precision matrix) and \code{S} (covariance matrix).
.fit_local_precision <- function(S_sub, adj_sub, n, method, ...) {
  p_sub <- nrow(S_sub)

  if (p_sub == 1L) {
    K <- matrix(1 / S_sub[1L, 1L], nrow = 1L, ncol = 1L)
    return(list(K = K, S = S_sub))
  }

  if (p_sub == 2L) {
    # For 2-node atoms: disconnected -> diagonal inverse;
    # connected -> full 2x2 inverse (avoids glasso subscript-out-of-bounds).
    if (adj_sub[1L, 2L] == 0L) {
      K <- diag(1 / diag(S_sub), nrow = 2L)
    } else {
      K <- .safe_solve(S_sub, label = "2-node atom")
    }
    return(list(K = K, S = S_sub))
  }

  if (identical(method, "ggm")) {
    if (!requireNamespace("ggm", quietly = TRUE)) {
      cli::cli_abort(
        "Package {.pkg ggm} is required for {.arg method = \"ggm\"}."
      )
    }

    amat <- adj_sub

    # Ensure adjacency matrix has vertex labels for ggm::fitConGraph
    if (is.null(colnames(S_sub))) {
      rownames(amat) <- colnames(amat) <- paste0("V", seq_len(p_sub))
    } else {
      rownames(amat) <- colnames(amat) <- colnames(S_sub)
    }
    diag(amat) <- 1L

    fit_obj <- tryCatch(
      do.call(
        ggm::fitConGraph,
        c(
          list(
            amat = amat,
            S = S_sub,
            n = n,
            cli = NULL
          ),
          list(...)
        )
      ),
      error = function(e) {
        # 若报错为奇异矩阵，自动加微小ridge再试一次
        msg <- conditionMessage(e)
        if (grepl("奇异|singular|倒条件数", msg, ignore.case = TRUE)) {
          S_ridge <- S_sub + diag(1e-6, nrow(S_sub))
          fit_obj2 <- tryCatch(
            do.call(
              ggm::fitConGraph,
              c(
                list(
                  amat = amat,
                  S = S_ridge,
                  n = n,
                  cli = NULL
                ),
                list(...)
              )
            ),
            error = function(e2) NULL
          )
          if (!is.null(fit_obj2)) {
            return(fit_obj2)
          }
        }
        cli::cli_abort(
          "ggm::fitConGraph failed on atom of size {.val {p_sub}}: {msg}"
        )
      }
    )

    K_sub <- .extract_precision(fit_obj)
    return(list(K = as.matrix(K_sub), S = S_sub))
  }

  if (!requireNamespace("glasso", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg glasso} is required for {.arg method = \"glasso\"}."
    )
  }

  non_edges_idx <- which(upper.tri(adj_sub) & (adj_sub == 0L), arr.ind = TRUE)
  if (length(non_edges_idx) == 0L) {
    zero_arg <- NULL
  } else {
    storage.mode(non_edges_idx) <- "integer"
    zero_arg <- non_edges_idx
  }

  # rho = 0 with graph constraints is constrained MLE (expected and correct).
  # Suppress the fixed "convergence problems" warning that glasso emits for rho=0.
  fit_obj <- tryCatch(
    suppressWarnings(
      do.call(
        glasso::glasso,
        c(
          list(
            s = S_sub,
            rho = 0,
            zero = zero_arg
          ),
          list(...)
        )
      )
    ),
    error = function(e) {
      cli::cli_abort(
        "glasso::glasso failed on atom of size {.val {p_sub}}: {conditionMessage(e)}"
      )
    }
  )

  if (is.null(fit_obj$wi) || !is.matrix(fit_obj$wi)) {
    cli::cli_abort("glasso output does not contain a valid precision matrix in `$wi`.")
  }

  list(K = as.matrix(fit_obj$wi), S = S_sub)
}