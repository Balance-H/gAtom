# gAtom

gAtom is an R package that exposes two graph decomposition routines backed by C + igraph.

## Exported R interfaces

- `get_minimal_collapsible(graph, r_nodes)`
- `decompose_atoms(graph)`

Both functions expect `graph` to be an `igraph` object.

## Build requirements

- R (with compilation toolchain)
- CRAN package igraph (provides bundled igraph C library)
- OpenMP support in compiler

## Install

```r
install.packages("igraph")
install.packages(".", repos = NULL, type = "source")
```

Or with devtools:

```r
devtools::load_all(".")
```

## Quick smoke test

```r
library(igraph)
library(gAtom)

g <- make_ring(6)
get_minimal_collapsible(g, c(1L, 3L))
out <- decompose_atoms(g)
str(out)
```
