library(igraph)

utils::globalVariables(c("Method", "FitMethod", "MissingRate", "K_Fnorm", "Sigma_Fnorm", "Time"))

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required for multivariate normal sampling.")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required to load source package.")
}

args      <- commandArgs(trailingOnly = FALSE)
file_arg  <- grep("^--file=", args, value = TRUE)

project_root <- NULL
if (length(file_arg) > 0L) {
  file_val <- sub("^--file=", "", file_arg[[1L]])
  if (nzchar(file_val)) {
    script_path <- tryCatch(
      normalizePath(file_val, winslash = "/", mustWork = TRUE),
      error = function(e) NULL
    )
    if (!is.null(script_path) && nzchar(script_path)) {
      project_root <- dirname(dirname(script_path))
    }
  }
}
if (is.null(project_root)) {
  project_root <- getwd()
}

# Suppress "masks gAtom::xxx" masking messages from load_all
suppressMessages(devtools::load_all(project_root, quiet = TRUE))

# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

read_animal_graph <- function(file_path) {
  edge_df <- utils::read.table(file_path, header = FALSE, sep = "\t")
  if (ncol(edge_df) < 2L) {
    stop("Animal-Network.txt must contain at least two columns.")
  }

  edges <- as.matrix(edge_df[, 1:2, drop = FALSE])
  storage.mode(edges) <- "integer"

  # File uses 0-based indices -> convert to R/igraph 1-based
  if (min(edges) == 0L) {
    edges <- edges + 1L
  }

  n_nodes <- max(edges)
  g <- igraph::make_empty_graph(n = n_nodes, directed = FALSE)
  g <- igraph::add_edges(g, as.vector(t(edges)))
  igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
}

ensure_connected <- function(g, seed = 20260427L) {
  comp <- igraph::components(g)
  if (comp$no <= 1L) {
    return(g)
  }
  set.seed(seed)
  rep_nodes <- vapply(
    seq_len(comp$no),
    function(i) sample(which(comp$membership == i), 1L),
    integer(1)
  )
  added <- cbind(rep_nodes[-length(rep_nodes)], rep_nodes[-1L])
  g <- igraph::add_edges(g, as.vector(t(added)))
  igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
}

generate_sparse_precision <- function(g, seed = 20260427L) {
  set.seed(seed)
  adj <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
  p   <- nrow(adj)

  K       <- matrix(0, nrow = p, ncol = p)
  eidx    <- which(upper.tri(adj) & (adj == 1L), arr.ind = TRUE)

  if (nrow(eidx) > 0L) {
    vals <- sample(c(-1, 1), nrow(eidx), replace = TRUE) *
      runif(nrow(eidx), min = 0.05, max = 0.25)
    for (k in seq_len(nrow(eidx))) {
      K[eidx[k, 1], eidx[k, 2]] <- vals[k]
      K[eidx[k, 2], eidx[k, 1]] <- vals[k]
    }
  }

  # Strict diagonal dominance => SPD
  diag(K) <- rowSums(abs(K)) + runif(p, min = 0.5, max = 1.5)
  K
}

extract_precision_from_fit <- function(fit_obj) {
  if (is.matrix(fit_obj)) return(fit_obj)
  if (!is.list(fit_obj)) stop("Unsupported fit output type.")
  for (nm in c("K", "Khat", "iSigma", "invSigma", "wi", "Omega", "precision")) {
    val <- fit_obj[[nm]]
    if (is.matrix(val)) return(val)
  }
  for (nm in c("Shat", "Sigma", "cov", "S")) {
    val <- fit_obj[[nm]]
    if (is.matrix(val)) return(solve(val))
  }
  stop("Cannot extract precision matrix from fit output.")
}

# Direct global fit (whole-graph S, single call to fitter)
fit_direct_global <- function(S, g, n, method) {
  adj <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
  storage.mode(adj) <- "integer"

  if (identical(method, "ggm")) {

    # Ensure adjacency matrix has vertex labels for ggm::fitConGraph
    amat <- adj
    diag(amat) <- 1L
    if (is.null(colnames(S))) {
      rownames(amat) <- colnames(amat) <- paste0("V", seq_len(nrow(amat)))
    } else {
      rownames(amat) <- colnames(amat) <- colnames(S)
    }
    fit_obj <- ggm::fitConGraph(
      amat = amat, S = S, n = n, cli = NULL
    )
    K <- extract_precision_from_fit(fit_obj)

  } else {                             # glasso
    zero_idx <- which(upper.tri(adj) & (adj == 0L), arr.ind = TRUE)
    zero_arg <- if (length(zero_idx) == 0L) NULL else {
      storage.mode(zero_idx) <- "integer"; zero_idx
    }
    fit_obj <- suppressWarnings(
      glasso::glasso(s = S, rho = 0, zero = zero_arg)
    )
    if (is.null(fit_obj$wi) || !is.matrix(fit_obj$wi)) {
      stop("glasso output missing precision matrix.")
    }
    K <- fit_obj$wi
  }

  K <- as.matrix(K)
  (K + t(K)) / 2        # symmetrize numerical noise
}

# ---------------------------------------------------------------------------
# Availability check for each method
# ---------------------------------------------------------------------------

available_methods <- character()
if (requireNamespace("glasso", quietly = TRUE)) {
  available_methods <- c(available_methods, "glasso")
}

if (requireNamespace("ggm", quietly = TRUE)) {
  available_methods <- c(available_methods, "ggm")
}
if (length(available_methods) == 0L) {
  stop("Need at least one of packages 'glasso' or 'ggm' installed.")
}

cat(sprintf("Available fitting methods: %s\n\n", paste(available_methods, collapse = ", ")))

# ---------------------------------------------------------------------------
# Load graph
# ---------------------------------------------------------------------------

network_path <- file.path(project_root, "example", "Animal-Network.txt")
if (!file.exists(network_path)) {
  stop("Cannot find example/Animal-Network.txt under project root.")
}

g <- read_animal_graph(network_path)
g <- ensure_connected(g, seed = 20260427L)

cat("=== Localized vs Direct fitting benchmark ===\n")
cat(sprintf("Graph:  |V|=%d, |E|=%d\n", igraph::vcount(g), igraph::ecount(g)))
cat(sprintf("Methods: %s\n", paste(available_methods, collapse = " / ")))

n_rep <- 3L
n_obs <- 1000L
tol_rel <- 5e-3


# ====== 数据生成与变量定义 ======
set.seed(20260427L)
K_true <- generate_sparse_precision(g, seed = 20260427L)
Sigma_true <- solve(K_true)
data <- MASS::mvrnorm(n = n_obs, mu = rep(0, nrow(K_true)), Sigma = Sigma_true)
data <- as.matrix(data)
storage.mode(data) <- "double"
colnames(data) <- paste0("V", seq_len(ncol(data)))
graph <- g

all_rows <- vector("list", n_rep * length(available_methods) * 2L)
row_id   <- 1L

















simulate_missing_data <- function(data_input, missing_rate) {
  if (!is.matrix(data_input) && !is.data.frame(data_input)) {
    stop("data_input must be a matrix or data.frame.")
  }
  if (is.function(data_input) || inherits(data_input, "closure")) {
    stop("data_input cannot be a function or closure.")
  }
  if (missing_rate < 0 || missing_rate > 1) {
    stop("missing_rate must be between 0 and 1.")
  }
  if (is.null(dim(data_input))) {
    stop("data_input must have dimensions (e.g., rows and columns).")
  }
  if (any(dim(data_input) == 0)) {
    stop("data_input must have non-zero rows and columns.")
  }
  if (is.environment(data_input)) {
    stop("data_input cannot be an environment.")
  }
  if (!is.numeric(data_input)) {
    stop("data_input must contain numeric values.")
  }
  if (is.character(data_input)) {
    stop("data_input cannot be a character matrix or data.frame.")
  }
  set.seed(20260427L)
  n <- nrow(data_input)
  p <- ncol(data_input)
  missing_indices <- which(runif(n * p) < missing_rate, arr.ind = TRUE)
  data_input[missing_indices] <- NA
  data_input
}

# 按用户要求：缺失率为x时，先选x比例的样本行，再对这些行的每个变量以0.01概率独立置NA
simulate_missing_data <- function(data_input, missing_rate, col_missing_prob = 0.01) {
  set.seed(20260427L)
  n <- nrow(data_input)
  p <- ncol(data_input)
  n_miss <- ceiling(n * missing_rate)
  if (n_miss == 0) return(data_input)
  miss_rows <- sample(seq_len(n), n_miss)
  for (i in miss_rows) {
    miss_cols <- which(runif(p) < col_missing_prob)
    if (length(miss_cols) > 0) {
      data_input[i, miss_cols] <- NA
    }
  }
  data_input
}

test_missing_data <- function() {
  missing_rates <- c(0, 0.2, 0.4, 0.6, 0.8)
  n_rep <- 3
  n_obs <- 1000
  fit_methods <- c("ggm", "glasso")
  results <- data.frame(
    Method = character(),
    FitMethod = character(),
    MissingRate = numeric(),
    K_Fnorm = numeric(),
    Sigma_Fnorm = numeric(),
    Time = numeric(),
    stringsAsFactors = FALSE
  )
  for (rate in missing_rates) {
    for (rep in 1:n_rep) {
      set.seed(20260427L + rep)
      # 重新生成数据，保证每次独立
      K_true <- generate_sparse_precision(g, seed = 20260427L)
      Sigma_true <- solve(K_true)
      data <- MASS::mvrnorm(n = n_obs, mu = rep(0, nrow(K_true)), Sigma = Sigma_true)
      data <- as.matrix(data)
      storage.mode(data) <- "double"
      colnames(data) <- paste0("V", seq_len(ncol(data)))
      if (rate == 0) {
        data_with_missing <- data
      } else {
        data_with_missing <- simulate_missing_data(data, rate)
      }
      for (fit_method in fit_methods) {
        # Direct method
        t1 <- proc.time()
        complete_data <- data_with_missing[complete.cases(data_with_missing), ]
        direct_K <- NA; Sigma_hat_direct <- NA; direct_K_err <- NA; direct_Sigma_err <- NA; direct_time <- NA
        if (nrow(complete_data) >= 2) {
          S_cc <- cov(complete_data)
          direct_K <- tryCatch(
            fit_direct_global(S_cc, graph, nrow(complete_data), method = fit_method),
            error = function(e) NA
          )
          Sigma_hat_direct <- tryCatch(solve(direct_K), error=function(e) NA)
          direct_K_err <- if (all(!is.na(direct_K))) sqrt(sum((direct_K-K_true)^2)) else NA
          direct_Sigma_err <- if (all(!is.na(Sigma_hat_direct))) sqrt(sum((Sigma_hat_direct-Sigma_true)^2)) else NA
          direct_time <- (proc.time()-t1)[[3]]
          results <- rbind(results, data.frame(Method="Direct", FitMethod=fit_method, MissingRate=rate, K_Fnorm=direct_K_err, Sigma_Fnorm=direct_Sigma_err, Time=direct_time))
        } else {
          results <- rbind(results, data.frame(Method="Direct", FitMethod=fit_method, MissingRate=rate, K_Fnorm=NA, Sigma_Fnorm=NA, Time=NA))
        }
        # Localized method
        t2 <- proc.time()
        localized_result <- tryCatch(
          fit_ggm_localized(data_with_missing, graph, method = fit_method),
          error = function(e) NA
        )
        local_K_err <- NA; local_Sigma_err <- NA; local_time <- NA
        if (is.list(localized_result) && !is.null(localized_result$K) && !is.null(localized_result$Sigma)) {
          local_K_err <- sqrt(sum((localized_result$K-K_true)^2))
          local_Sigma_err <- sqrt(sum((localized_result$Sigma-Sigma_true)^2))
          local_time <- (proc.time()-t2)[[3]]
          results <- rbind(results, data.frame(Method="Localized", FitMethod=fit_method, MissingRate=rate, K_Fnorm=local_K_err, Sigma_Fnorm=local_Sigma_err, Time=local_time))
        } else {
          results <- rbind(results, data.frame(Method="Localized", FitMethod=fit_method, MissingRate=rate, K_Fnorm=NA, Sigma_Fnorm=NA, Time=NA))
        }
      }
    }
  }
  # 汇总平均
  library(dplyr)
  summary_tbl <- results %>%
    group_by(Method, FitMethod, MissingRate) %>%
    summarise(K_Fnorm=mean(K_Fnorm, na.rm=TRUE),
              Sigma_Fnorm=mean(Sigma_Fnorm, na.rm=TRUE),
              Time=mean(Time, na.rm=TRUE), .groups="drop")
  print(as.data.frame(summary_tbl), row.names=FALSE)
}
test_missing_data()