library(igraph)

dyn.load(file.path(deps_dir, "gAtom.dll"), local = FALSE)
source(file.path(repo_root, "R", "gAtom.R"))

g <- make_ring(6)
#g <- make_full_graph(5)
print(g)

out_min <- get_minimal_collapsible(g, c(2L, 3L))
cat("out_min:", out_min, "\n")

out_dec <- decompose_atoms(g)
cat("out_dec:", str(out_dec), "\n")