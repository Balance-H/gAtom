# decom_h.dll depends on igraph.dll, which is installed via vcpkg-master and stored in this directory
Sys.setenv(PATH = paste("D:/vcpkg-master/installed/x64-windows/bin", Sys.getenv("PATH"), sep=";"))

# decom_h.dll contains the recursive_decom function used for igraph graph decomposition
dyn.load("decom_h.dll")

# Verify whether the recursive_decom function has been successfully loaded
is.loaded("recursive_decom")  # TRUE

# Load the igraph package
library(igraph)

# Download grips-main from https://github.com/hojsgaard/gRips, extract it, and install/load the grips package locally
devtools::clean_dll("./gRips")
devtools::load_all("./gRips")
library(gRbase)
library(tidyverse)
library(dplyr)
library(MASS)
library(htmlTable)
library(doBy)

# Decomposition function interface (calls C/C++ implementation via .Call)
create_igraph_external_ptr <- function(edges, n_vertices, directed = FALSE) {
  .Call("create_igraph_external_ptr",
        as.integer(edges),        # edge list as integer vector
        as.integer(n_vertices),   # number of vertices
        as.logical(directed))     # whether the graph is directed
}

# Recursive graph decomposition function (CMSA or other supported methods)
recursive_decom_r <- function(graph_ptr, method = "cmsa") {
  .Call("recursive_decom_r",
        graph_ptr,                # external pointer to igraph object
        as.character(method))     # decomposition method name
}


# Ensure that the graph is connected by adding edges if necessary
make_connect_graph <- function(edges, n_nodes) {
  
  # Create an undirected graph from the edge list
  g <- graph_from_edgelist(matrix(edges, ncol = 2, byrow = TRUE), directed = FALSE)
  
  # Identify connected components
  comp <- components(g)
  k <- comp$no                # number of connected components
  membership <- comp$membership  # component membership of each vertex
  
  # If already connected, return the original edge matrix
  if (k == 1) {
    return(matrix(edges, ncol = 2, byrow = TRUE))
  }
  
  # Select one representative vertex from each component
  rep_nodes <- sapply(1:k, function(i) which(membership == i)[1])
  
  # Construct additional edges to connect components sequentially
  added_edges <- matrix(NA, nrow = k - 1, ncol = 2)
  for (i in 1:(k - 1)) {
    added_edges[i, ] <- c(rep_nodes[i], rep_nodes[i + 1])
  }
  
  # Return the complete edge list (original edges + added edges)
  edges_full <- rbind(matrix(edges, ncol = 2, byrow = TRUE), added_edges)
  return(edges_full)
}

generate_sigma_and_theta <- function(g, n, seed = 123) {
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Get adjacency matrix (non-sparse form)
  adj <- as_adjacency_matrix(g, sparse = FALSE)
  
  # Initialize precision matrix theta
  theta <- matrix(0, n, n)
  
  # Randomly assign positive values to diagonal elements
  diag(theta) <- runif(n, 1, 2)
  
  # Randomly assign non-zero values to adjacent edges (non-zero precision matrix)
  edges <- which(adj != 0, arr.ind = TRUE)
  for(k in 1:nrow(edges)){
    i <- edges[k,1]
    j <- edges[k,2]
    if(i < j){  # Assign to upper triangle to ensure symmetry
      val <- runif(1, 0.1, 0.5)
      theta[i,j] <- val
      theta[j,i] <- val
    }
  }
  
  # Adjust diagonal elements to ensure positive definiteness
  diag(theta) <- diag(theta) + rowSums(abs(theta))
  
  # Invert to get covariance matrix sigma
  sigma <- solve(theta)
  
  # Return both theta and sigma
  return(list(theta = theta, sigma = sigma))
}

generate_masked_rows <- function(data, n_vars, mask_fraction, seed = 123, n_na_per_row = 40) {
  set.seed(seed)
  
  data <- as.matrix(data)
  
  data_masked <- data
  
  # Only randomly set NA if mask_fraction > 0
  n_mask_rows <- round(mask_fraction * 1093)
  if (n_mask_rows > 0) {
    rows_to_mask <- sample(seq_len(1093), n_mask_rows)
    
    for (r in rows_to_mask) {
      na_cols <- sample(seq_len(n_vars), n_na_per_row)
      data_masked[r, na_cols] <- NA
    }
  }
  return(data_masked)
}

# Specify the folder path containing network files
folder_path <- "./examples"

# List of network files to be processed
file_list <- c(
  "Animal-Network.txt"
)


# Global parameters
global_args <- list(
  eps   = 1e-3,
  nobs  = 1093,
  maxit = 50000
)

# Experimental parameter setup
design <- expand.grid(
  method = c("cov", "ncd"),
  approach = c("Decom", "Direct"),
  mask  = seq(0.1, 0.7, by = 0.15),
  rep    = 1:20,
  marg   = "edge",
  dat    = c("n01"),#
  eng    = "cpp",
  bench  = "B3",
  network = file_list,
  stringsAsFactors = FALSE
)

# Generate a unique identifier for each row
rownames(design) <- with(
  design,
  interaction(rep, approach, method, marg, dat, mask, bench, network, sep = "_")
)

# Split by row into a list (each element is one design combination)
des_lst <- split(design, seq_len(nrow(design)))

des_lst <- lapply(des_lst, function(df){
  rownames(df) <- seq_len(nrow(df))
  df
})

# Inspect the split list
des_lst

# Setting of the training function
train_network <- function(experiment_row, folder_path, global_args) {
  fname <- experiment_row$network
  file_path <- file.path(folder_path, fname)
  edge_df <- read.table(file_path, header = FALSE, sep = "\t")
  colnames(edge_df) <- c("from", "to")
  edges <- as.integer(t(edge_df))
  n_nodes <- max(edges)
  edges_full <- make_connect_graph(edges, n_nodes)
  edges_full_vec <- as.integer(t(edges_full))
  
  edges_list_full <- split(edges_full, seq(nrow(edges_full))) 
  edges_list_full <- lapply(edges_list_full, function(x) unname(as.numeric(x)))
  
  g_ptr <- create_igraph_external_ptr(edges_full_vec, n_nodes, directed = FALSE)
  g_igraph <- make_graph(edges = edges_full_vec, n = n_nodes, directed = FALSE)
  
  res_matrix <- generate_sigma_and_theta(g_igraph, n_nodes, seed = 2026 + experiment_row$rep)
    
  theta <- res_matrix$theta  # precision matrix
  sigma <- res_matrix$sigma  # covariance matrix
  
    
  # Set parameters
  n_obs <- 1093        # number of observations (rows)
  n_vars <- n_nodes        # number of variables (columns)
  mu <- rep(0, n_vars) # mean vector N(0, sigma)
  
  # Generate multivariate normal distribution sample
  sample_matrix <- mvrnorm(n = n_obs, mu = mu, Sigma = sigma)
  
  # Convert to data frame
  data <- as.data.frame(sample_matrix)
  data_masked <- generate_masked_rows(data, n_vars = n_nodes, mask_fraction = experiment_row$mask, seed = 2026+experiment_row$rep)
  
  set.seed(2026 + experiment_row$rep)
  
  if (!is.null(experiment_row$approach) && tolower(experiment_row$approach) == "direct") {
    data_complete <- data_masked[
      complete.cases(data_masked),
      , drop = FALSE
    ]
    
    cov_data <- cov(data_complete)
    
    gl <- edges_list_full  # list of length-2 vectors
    arg_list <- list(
      S = cov_data,
      formula = gl,
      method = experiment_row$method,
      nobs = global_args$nobs,
      eps  = global_args$eps,
      maxit = global_args$maxit,
      aux = list(engine = experiment_row$eng)
    )
    
    fit_fn <- set_default(fit_ggm, c(arg_list, global_args))
    
    t0 <- Sys.time()
    fit_out <- do.call(fit_fn, list())
    total_time <- as.numeric(Sys.time() - t0, units = "secs")
    sigma_hat <- fit_out$Sigma
    
    sigma_diff <- norm(sigma_hat - sigma, type = "F")
    
    cat(sprintf("approach=%s, method=%s, mask=%.2f, time=%.2f secs, cov-F-norm=%.4f\n",
                experiment_row$approach,
                experiment_row$method,
                experiment_row$mask,
                as.numeric(total_time, units = "secs"),
                sigma_diff))
    
    # adj: 0-1 adjacency matrix, 1 indicates allowed edges
    
    return(data.frame(
      mask_fraction = experiment_row$mask,
      approach = experiment_row$approach,
      method   = experiment_row$method,
      rep      = experiment_row$rep,
      marg     = experiment_row$marg,
      dat      = experiment_row$dat,
      eng      = experiment_row$eng,
      bench    = experiment_row$bench,
      network  = experiment_row$network,
      Sigma        = I(list(sigma_hat)),
      sigma_diff = sigma_diff,
      time     = total_time
    ))
  }
  
  t0 <- Sys.time()
  
  # Decomposition of undirected graphs into atoms and separators
  res <- recursive_decom_r(g_ptr, method = "cmsa")
  atoms <- res[[1]]
  separators <- res[[2]]
  
  # Construct training atoms (size ≥3)
  atom_edges_list <- lapply(atoms, function(atom_nodes) {
    nodes_global <- sort(atom_nodes)
    if(length(nodes_global) < 4) return(NULL)  
    sub_edges <- edge_df[edge_df$from %in% atom_nodes & edge_df$to %in% atom_nodes, ]
    node_map <- setNames(seq_along(atom_nodes), as.character(atom_nodes))
    edges_local <- t(apply(sub_edges, 1, function(e) as.numeric(node_map[as.character(e)])))
    list(edges_local = edges_local, nodes_global = nodes_global)
  })
  atom_edges_list <- Filter(Negate(is.null), atom_edges_list)
  
  make_arg_list <- function(edges_list, dat_matrix, type_label) {
    lapply(edges_list, function(sub) {
      sub_data_full <- dat_matrix[
        complete.cases(dat_matrix[, sub$nodes_global]),
        sub$nodes_global,
        drop = FALSE
      ]
      
      S <- cov(sub_data_full)
      formula <- lapply(seq_len(nrow(sub$edges_local)), function(i) sub$edges_local[i, ])
      list(
        S = S,
        formula = formula,
        method = experiment_row$method,
        type = type_label,
        nodes_global = sub$nodes_global
      )
    })
  }
  
  atom_arg_lst <- make_arg_list(atom_edges_list, data_masked, "atom")
  
  fn_lst <- lapply(atom_arg_lst, function(lsti){
    set_default(fit_ggm, c(lsti, global_args))
  })
  
  res_lst <- lapply(fn_lst, function(fit_fn) {do.call(fit_fn, list())})
  
  # Constructing the global precision matrix Q
  Q <- matrix(0, n_nodes, n_nodes)
  
  # atom size ≥3:
  for (i in seq_along(res_lst)) {
    nodes <- atom_arg_lst[[i]]$nodes_global
    Kinv_local <- res_lst[[i]]$K
    Q[nodes, nodes] <- Q[nodes, nodes] + Kinv_local
  }
  
  # Handle small atoms
  small_atoms <- atoms[sapply(atoms, length) < 4]
  for(nodes in small_atoms) {
    if(length(nodes) == 1) {
      Q[nodes, nodes] <- Q[nodes, nodes] + 1 / (var(data_masked[, nodes], na.rm = TRUE) + 1e-8)
    } else {
      sub_data <- data_masked[, nodes, drop = FALSE]  # extract multiple columns
      sub_data_full <- sub_data[complete.cases(sub_data), ]  # keep rows without NA
      
      S <- cov(sub_data_full)
      S_sub <- as.matrix(S)
      storage.mode(S_sub) <- "double"
      K_sub <- solve(S_sub + 1e-8 * diag(length(nodes)))
      Q[nodes, nodes] <- Q[nodes, nodes] + K_sub
    }
  }
  
  # Handle separators
  for(nodes in separators) {
    sub_data <- data_masked[, nodes, drop = FALSE]
    sub_data_full <- sub_data[complete.cases(sub_data), , drop = FALSE]
    
    S <- cov(sub_data_full)
    S <- as.matrix(S)
    storage.mode(S) <- "double"
    
    K_sub <- solve(S + 1e-8 * diag(ncol(S)))
    
    Q[nodes, nodes] <- Q[nodes, nodes] - K_sub
  }
  
  total_time <- Sys.time() - t0
  sigma_hat = solve(Q)
  
  sigma_diff <- norm(sigma_hat - sigma, type = "F")
  
  cat(sprintf("approach=%s, method=%s, mask=%.2f, time=%.2f secs, cov-F-norm=%.5f\n",
              experiment_row$approach,
              experiment_row$method,
              experiment_row$mask,
              as.numeric(total_time, units = "secs"),
              sigma_diff))
  
  data.frame(
    mask_fraction = experiment_row$mask,
    approach = experiment_row$approach,
    method = experiment_row$method,
    rep    = experiment_row$rep,
    marg   = experiment_row$marg,
    dat    = experiment_row$dat,
    eng    = experiment_row$eng,
    bench  = experiment_row$bench,
    network = experiment_row$network,
    Sigma  = I(list(sigma_hat)), 
    sigma_diff = sigma_diff,
    time   = as.numeric(total_time, units="secs")
  )
}

res_lst <- lapply(seq_along(des_lst), function(i) {
  train_network(des_lst[[i]][1, ], folder_path, global_args)
})
#res_lst

all_df <- do.call(rbind, res_lst)

result <- all_df %>%
  group_by(network, mask_fraction,  method, dat, approach) %>%
  summarise(
    F_norm     = round(mean(sigma_diff, na.rm = TRUE), 6), 
    n_reps     = n(),
    .groups    = 'drop'
  ) %>%
  arrange(mask_fraction, network, dat, method, approach)

options(pillar.sigfig = 7)
print(result)


