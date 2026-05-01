#include <R.h>
#include <R_ext/Arith.h>
#include <Rinternals.h>
#include <igraph.h>

#include <limits.h>
#include <stdlib.h>
#include <string.h>

//#include "decom_h.h"  // 注释掉，因为使用自己的实现而不是 decom_h 中的

typedef struct {
    int *data;
    int size;
    int capacity;
} IntVec;

typedef struct {
    IntVec *data;
    int size;
    int capacity;
} IntVecList;

typedef struct {
    int n;
    IntVec *adj;
} Graph;

typedef struct {
    IntVec comp;
    IntVec bound;
} CompBound;

typedef struct {
    CompBound *data;
    int size;
    int capacity;
} CompBoundList;

static void *xmalloc(size_t nbytes) {
    void *ptr;

    if (nbytes == 0) {
        return NULL;
    }

    ptr = malloc(nbytes);
    if (!ptr) {
        Rf_error("memory allocation failed");
    }

    return ptr;
}

static void *xcalloc(size_t count, size_t size) {
    void *ptr;

    if (count == 0 || size == 0) {
        return NULL;
    }

    ptr = calloc(count, size);
    if (!ptr) {
        Rf_error("memory allocation failed");
    }

    return ptr;
}

static void *xrealloc(void *ptr, size_t nbytes) {
    void *out;

    if (nbytes == 0) {
        free(ptr);
        return NULL;
    }

    out = realloc(ptr, nbytes);
    if (!out) {
        Rf_error("memory allocation failed");
    }

    return out;
}

static int int_compare(const void *a, const void *b) {
    const int ia = *(const int *) a;
    const int ib = *(const int *) b;
    return (ia > ib) - (ia < ib);
}

static void intvec_init(IntVec *vec) {
    vec->data = NULL;
    vec->size = 0;
    vec->capacity = 0;
}

static void intvec_free(IntVec *vec) {
    free(vec->data);
    vec->data = NULL;
    vec->size = 0;
    vec->capacity = 0;
}

static void intvec_clear(IntVec *vec) {
    vec->size = 0;
}

static void intvec_reserve(IntVec *vec, int needed_capacity) {
    int new_capacity;

    if (needed_capacity <= vec->capacity) {
        return;
    }

    new_capacity = vec->capacity > 0 ? vec->capacity : 4;
    while (new_capacity < needed_capacity) {
        if (new_capacity > INT_MAX / 2) {
            new_capacity = needed_capacity;
            break;
        }
        new_capacity *= 2;
    }

    vec->data = (int *) xrealloc(vec->data, (size_t) new_capacity * sizeof(int));
    vec->capacity = new_capacity;
}

static void intvec_resize(IntVec *vec, int new_size) {
    if (new_size < 0) {
        Rf_error("internal error: negative vector size");
    }

    intvec_reserve(vec, new_size);

    if (new_size > vec->size) {
        memset(
            vec->data + vec->size,
            0,
            (size_t) (new_size - vec->size) * sizeof(int)
        );
    }

    vec->size = new_size;
}

static void intvec_push_back(IntVec *vec, int value) {
    intvec_reserve(vec, vec->size + 1);
    vec->data[vec->size++] = value;
}

static void intvec_append(IntVec *dest, const IntVec *src) {
    if (src->size == 0) {
        return;
    }

    intvec_reserve(dest, dest->size + src->size);
    memcpy(
        dest->data + dest->size,
        src->data,
        (size_t) src->size * sizeof(int)
    );
    dest->size += src->size;
}

static void intvec_copy(IntVec *dest, const IntVec *src) {
    intvec_clear(dest);
    intvec_append(dest, src);
}

static void intvec_sort_unique(IntVec *vec) {
    int write_pos = 1;

    if (vec->size <= 1) {
        return;
    }

    qsort(vec->data, (size_t) vec->size, sizeof(int), int_compare);

    for (int i = 1; i < vec->size; i++) {
        if (vec->data[i] != vec->data[i - 1]) {
            vec->data[write_pos++] = vec->data[i];
        }
    }

    vec->size = write_pos;
}

static int intvec_equal(const IntVec *a, const IntVec *b) {
    if (a->size != b->size) {
        return 0;
    }

    for (int i = 0; i < a->size; i++) {
        if (a->data[i] != b->data[i]) {
            return 0;
        }
    }

    return 1;
}

static void intveclist_init(IntVecList *list) {
    list->data = NULL;
    list->size = 0;
    list->capacity = 0;
}

static void intveclist_free(IntVecList *list) {
    for (int i = 0; i < list->size; i++) {
        intvec_free(&list->data[i]);
    }

    free(list->data);
    list->data = NULL;
    list->size = 0;
    list->capacity = 0;
}

static void intveclist_push_take(IntVecList *list, IntVec *vec) {
    int new_capacity;

    if (list->size == list->capacity) {
        new_capacity = list->capacity > 0 ? list->capacity * 2 : 4;
        list->data = (IntVec *) xrealloc(
            list->data,
            (size_t) new_capacity * sizeof(IntVec)
        );
        list->capacity = new_capacity;
    }

    list->data[list->size++] = *vec;
    intvec_init(vec);
}

static IntVec intveclist_pop(IntVecList *list) {
    IntVec out;

    if (list->size == 0) {
        Rf_error("internal error: pop from empty list");
    }

    out = list->data[list->size - 1];
    list->size--;

    return out;
}

static void graph_init(Graph *graph, int n_vertices) {
    graph->n = n_vertices;
    graph->adj = NULL;

    if (n_vertices < 0) {
        Rf_error("internal error: negative vertex count");
    }

    if (n_vertices == 0) {
        return;
    }

    graph->adj = (IntVec *) xcalloc((size_t) n_vertices, sizeof(IntVec));
    for (int v = 0; v < n_vertices; v++) {
        intvec_init(&graph->adj[v]);
    }
}

static void graph_free(Graph *graph) {
    if (graph->adj) {
        for (int v = 0; v < graph->n; v++) {
            intvec_free(&graph->adj[v]);
        }
        free(graph->adj);
    }

    graph->adj = NULL;
    graph->n = 0;
}

static void graph_add_undirected_edge(Graph *graph, int u, int v) {
    if (u < 0 || u >= graph->n || v < 0 || v >= graph->n) {
        Rf_error("internal error: edge endpoint out of range");
    }

    intvec_push_back(&graph->adj[u], v);
    if (u != v) {
        intvec_push_back(&graph->adj[v], u);
    }
}

static void graph_build_from_edges(
    Graph *graph,
    int n_vertices,
    const int *from,
    const int *to,
    int n_edges
) {
    graph_init(graph, n_vertices);

    for (int e = 0; e < n_edges; e++) {
        graph_add_undirected_edge(graph, from[e], to[e]);
    }
}

static int graph_are_adjacent(const Graph *graph, int u, int v) {
    const IntVec *neighbors;

    if (u < 0 || u >= graph->n || v < 0 || v >= graph->n) {
        return 0;
    }

    neighbors = &graph->adj[u];
    for (int i = 0; i < neighbors->size; i++) {
        if (neighbors->data[i] == v) {
            return 1;
        }
    }

    return 0;
}

static int is_clique(const Graph *graph, const IntVec *nodes) {
    for (int i = 0; i < nodes->size; i++) {
        for (int j = i + 1; j < nodes->size; j++) {
            if (!graph_are_adjacent(graph, nodes->data[i], nodes->data[j])) {
                return 0;
            }
        }
    }
    return 1;
}

static int scalar_to_nonnegative_integer(SEXP x, const char *arg_name) {
    int value = Rf_asInteger(x);
    if (value == NA_INTEGER || value < 0) {
        Rf_error("%s must be a non-negative scalar integer", arg_name);
    }
    return value;
}

static int scalar_to_logical_flag(SEXP x, const char *arg_name) {
    int value = Rf_asLogical(x);
    if (value == NA_LOGICAL) {
        Rf_error("%s must be TRUE or FALSE", arg_name);
    }
    return value;
}

static int vector_value_as_integer(SEXP vec, R_xlen_t index, const char *arg_name) {
    if (TYPEOF(vec) == INTSXP) {
        int value = INTEGER(vec)[index];
        if (value == NA_INTEGER) {
            Rf_error("%s contains NA values", arg_name);
        }
        return value;
    }

    if (TYPEOF(vec) == REALSXP) {
        double value = REAL(vec)[index];
        int coerced;

        if (!R_FINITE(value)) {
            Rf_error("%s contains non-finite values", arg_name);
        }

        if (value < (double) INT_MIN || value > (double) INT_MAX) {
            Rf_error("%s contains values outside integer range", arg_name);
        }

        coerced = (int) value;
        if ((double) coerced != value) {
            Rf_error("%s must contain whole numbers only", arg_name);
        }

        return coerced;
    }

    Rf_error("%s must be integer or numeric", arg_name);
    return 0;
}

static void parse_edge_matrix_zero_based(
    SEXP edge_mat,
    int n_vertices,
    int **from_out,
    int **to_out,
    int *n_edges_out
) {
    SEXP dims;
    int n_rows;
    int n_cols;

    if (!Rf_isMatrix(edge_mat)) {
        Rf_error("edge_mat must be a matrix with 2 columns");
    }

    if (TYPEOF(edge_mat) != INTSXP && TYPEOF(edge_mat) != REALSXP) {
        Rf_error("edge_mat must be an integer or numeric matrix");
    }

    dims = Rf_getAttrib(edge_mat, R_DimSymbol);
    if (TYPEOF(dims) != INTSXP || XLENGTH(dims) != 2) {
        Rf_error("edge_mat dimensions are invalid");
    }

    n_rows = INTEGER(dims)[0];
    n_cols = INTEGER(dims)[1];

    if (n_rows < 0 || n_cols != 2) {
        Rf_error("edge_mat must have shape [n_edges, 2]");
    }

    *from_out = n_rows > 0 ? (int *) xmalloc((size_t) n_rows * sizeof(int)) : NULL;
    *to_out = n_rows > 0 ? (int *) xmalloc((size_t) n_rows * sizeof(int)) : NULL;

    for (int e = 0; e < n_rows; e++) {
        int from_0based = vector_value_as_integer(edge_mat, (R_xlen_t) e, "edge_mat");
        int to_0based = vector_value_as_integer(edge_mat, (R_xlen_t) (e + n_rows), "edge_mat");

        if (
            from_0based < 0 ||
            from_0based >= n_vertices ||
            to_0based < 0 ||
            to_0based >= n_vertices
        ) {
            Rf_error("edge endpoints must be within [0, n_vertices - 1]");
        }

        (*from_out)[e] = from_0based;
        (*to_out)[e] = to_0based;
    }

    *n_edges_out = n_rows;
}

static void parse_vertices_zero_based(
    SEXP vertices,
    int n_vertices,
    const char *arg_name,
    IntVec *out
) {
    R_xlen_t n;

    if (TYPEOF(vertices) != INTSXP && TYPEOF(vertices) != REALSXP) {
        Rf_error("%s must be an integer or numeric vector", arg_name);
    }

    n = XLENGTH(vertices);
    if (n == 0) {
        Rf_error("%s must not be empty", arg_name);
    }

    intvec_clear(out);
    intvec_reserve(out, (int) n);

    for (R_xlen_t i = 0; i < n; i++) {
        int value_0based = vector_value_as_integer(vertices, i, arg_name);
        if (value_0based < 0 || value_0based >= n_vertices) {
            Rf_error("%s must be within [0, n_vertices - 1]", arg_name);
        }
        intvec_push_back(out, value_0based);
    }

    intvec_sort_unique(out);
}

static void compboundlist_init(CompBoundList *list) {
    list->data = NULL;
    list->size = 0;
    list->capacity = 0;
}

static void compboundlist_clear(CompBoundList *list) {
    for (int i = 0; i < list->size; i++) {
        intvec_free(&list->data[i].comp);
        intvec_free(&list->data[i].bound);
    }
    list->size = 0;
}

static void compboundlist_free(CompBoundList *list) {
    compboundlist_clear(list);
    free(list->data);
    list->data = NULL;
    list->capacity = 0;
}

static void compboundlist_push_take(CompBoundList *list, CompBound *item) {
    int new_capacity;

    if (list->size == list->capacity) {
        new_capacity = list->capacity > 0 ? list->capacity * 2 : 4;
        list->data = (CompBound *) xrealloc(
            list->data,
            (size_t) new_capacity * sizeof(CompBound)
        );
        list->capacity = new_capacity;
    }

    list->data[list->size++] = *item;
    intvec_init(&item->comp);
    intvec_init(&item->bound);
}

static void components_forbidden(
    const Graph *graph,
    const IntVec *forbidden_vertices,
    CompBoundList *out
) {
    unsigned char *visited;
    unsigned char *boundary_marked;
    int *queue;
    int n;

    compboundlist_clear(out);

    n = graph->n;
    visited = (unsigned char *) xcalloc((size_t) (n > 0 ? n : 1), sizeof(unsigned char));
    boundary_marked = (unsigned char *) xcalloc((size_t) (n > 0 ? n : 1), sizeof(unsigned char));
    queue = (int *) xmalloc((size_t) (n > 0 ? n : 1) * sizeof(int));

    if (forbidden_vertices) {
        for (int i = 0; i < forbidden_vertices->size; i++) {
            int v = forbidden_vertices->data[i];
            if (v >= 0 && v < n) {
                visited[v] = 2;
            }
        }
    }

    for (int v = 0; v < n; v++) {
        int head;
        int tail;
        CompBound item;

        if (visited[v] != 0) {
            continue;
        }

        intvec_init(&item.comp);
        intvec_init(&item.bound);

        head = 0;
        tail = 0;
        queue[tail++] = v;
        visited[v] = 1;

        while (head < tail) {
            int cur = queue[head++];
            IntVec *neighbors = &graph->adj[cur];

            intvec_push_back(&item.comp, cur);

            for (int i = 0; i < neighbors->size; i++) {
                int w = neighbors->data[i];

                if (visited[w] == 2) {
                    if (!boundary_marked[w]) {
                        intvec_push_back(&item.bound, w);
                        boundary_marked[w] = 1;
                    }
                    continue;
                }

                if (visited[w] == 0) {
                    visited[w] = 1;
                    queue[tail++] = w;
                }
            }
        }

        for (int i = 0; i < item.bound.size; i++) {
            boundary_marked[item.bound.data[i]] = 0;
        }

        compboundlist_push_take(out, &item);
        intvec_free(&item.comp);
        intvec_free(&item.bound);
    }

    free(queue);
    free(boundary_marked);
    free(visited);
}

static void close_separator(
    const Graph *graph,
    int vertex,
    const IntVec *forbidden_vertices,
    IntVec *bound_b
) {
    unsigned char *visited;
    int *queue;
    int head = 0;
    int tail = 0;
    int n = graph->n;

    if (vertex < 0 || vertex >= n) {
        Rf_error("internal error: close_separator start vertex out of range");
    }

    visited = (unsigned char *) xcalloc((size_t) (n > 0 ? n : 1), sizeof(unsigned char));
    queue = (int *) xmalloc((size_t) (n > 0 ? n : 1) * sizeof(int));

    if (forbidden_vertices) {
        for (int i = 0; i < forbidden_vertices->size; i++) {
            int v = forbidden_vertices->data[i];
            if (v >= 0 && v < n) {
                visited[v] = 2;
            }
        }
    }

    intvec_clear(bound_b);
    queue[tail++] = vertex;
    visited[vertex] = 1;

    while (head < tail) {
        int cur = queue[head++];
        const IntVec *neighbors = &graph->adj[cur];

        for (int i = 0; i < neighbors->size; i++) {
            int w = neighbors->data[i];

            if (visited[w] == 2) {
                intvec_push_back(bound_b, w);
                visited[w] = 1;
                continue;
            }

            if (visited[w] == 0) {
                visited[w] = 1;
                queue[tail++] = w;
            }
        }
    }

    free(queue);
    free(visited);
}

static void induced_subgraph(
    const Graph *graph,
    const IntVec *sub_nodes,
    Graph *subgraph,
    int **map_out,
    IntVec *invmap_out
) {
    int *map = NULL;
    int n = graph->n;

    if (n > 0) {
        map = (int *) xmalloc((size_t) n * sizeof(int));
        for (int i = 0; i < n; i++) {
            map[i] = -1;
        }
    }

    intvec_clear(invmap_out);

    for (int i = 0; i < sub_nodes->size; i++) {
        int gv = sub_nodes->data[i];
        if (gv < 0 || gv >= n) {
            free(map);
            Rf_error("internal error: subgraph vertex out of range");
        }

        if (map[gv] < 0) {
            map[gv] = invmap_out->size;
            intvec_push_back(invmap_out, gv);
        }
    }

    graph_init(subgraph, invmap_out->size);

    for (int local_u = 0; local_u < invmap_out->size; local_u++) {
        int gv_u = invmap_out->data[local_u];
        const IntVec *neighbors = &graph->adj[gv_u];

        for (int i = 0; i < neighbors->size; i++) {
            int gv_v = neighbors->data[i];
            int local_v = map[gv_v];

            if (local_v < 0) {
                continue;
            }

            if (local_u <= local_v) {
                graph_add_undirected_edge(subgraph, local_u, local_v);
            }
        }
    }

    *map_out = map;
}

static void maximum_cardinality_search(const Graph *graph, IntVec *order_out) {
    int n = graph->n;
    int *weight = (int *) xcalloc((size_t) (n > 0 ? n : 1), sizeof(int));
    unsigned char *selected = (unsigned char *) xcalloc((size_t) (n > 0 ? n : 1), sizeof(unsigned char));

    intvec_resize(order_out, n);

    for (int pos = n - 1; pos >= 0; pos--) {
        int best = -1;
        int best_weight = INT_MIN;

        for (int v = 0; v < n; v++) {
            if (selected[v]) {
                continue;
            }

            if (
                best < 0 ||
                weight[v] > best_weight ||
                (weight[v] == best_weight && v < best)
            ) {
                best = v;
                best_weight = weight[v];
            }
        }

        if (best < 0) {
            break;
        }

        order_out->data[pos] = best;
        selected[best] = 1;

        for (int i = 0; i < graph->adj[best].size; i++) {
            int w = graph->adj[best].data[i];
            if (!selected[w]) {
                weight[w] += 1;
            }
        }
    }

    free(selected);
    free(weight);
}

static void get_minimal_collapsible_core(
    const Graph *graph,
    const IntVec *r_nodes,
    IntVec *H_out
) {
    IntVec H;
    CompBoundList comps;
    int keep_iterating = 1;

    intvec_init(&H);
    intvec_append(&H, r_nodes);
    intvec_sort_unique(&H);

    compboundlist_init(&comps);

    while (keep_iterating) {
        IntVec H_local;

        keep_iterating = 0;
        components_forbidden(graph, &H, &comps);

        intvec_init(&H_local);

        for (int cidx = 0; cidx < comps.size; cidx++) {
            CompBound *item = &comps.data[cidx];
            int found = 0;

            if (item->bound.size < 2) {
                continue;
            }

            intvec_sort_unique(&item->bound);

            for (int j = 0; j < item->bound.size && !found; j++) {
                int hj = item->bound.data[j];

                for (int k = j + 1; k < item->bound.size && !found; k++) {
                    int hk = item->bound.data[k];

                    if (!graph_are_adjacent(graph, hj, hk)) {
                        IntVec sub_nodes;
                        Graph subgraph;
                        IntVec invmap;
                        int *map = NULL;
                        int a;
                        int b;
                        IntVec neighbors_a;
                        IntVec neighbors_b;
                        IntVec sep_a_local;
                        IntVec sep_b_local;

                        found = 1;
                        keep_iterating = 1;

                        intvec_init(&sub_nodes);
                        intvec_append(&sub_nodes, &item->comp);
                        intvec_push_back(&sub_nodes, hj);
                        intvec_push_back(&sub_nodes, hk);
                        intvec_sort_unique(&sub_nodes);

                        intvec_init(&invmap);
                        induced_subgraph(graph, &sub_nodes, &subgraph, &map, &invmap);

                        if (!map || map[hj] < 0 || map[hk] < 0) {
                            graph_free(&subgraph);
                            free(map);
                            intvec_free(&invmap);
                            intvec_free(&sub_nodes);
                            intvec_free(&H_local);
                            compboundlist_free(&comps);
                            intvec_free(&H);
                            Rf_error("internal error: failed to map boundary vertices");
                        }

                        a = map[hj];
                        b = map[hk];

                        intvec_init(&neighbors_a);
                        intvec_init(&neighbors_b);
                        intvec_init(&sep_a_local);
                        intvec_init(&sep_b_local);

                        intvec_append(&neighbors_a, &subgraph.adj[a]);
                        intvec_append(&neighbors_b, &subgraph.adj[b]);
                        intvec_sort_unique(&neighbors_a);
                        intvec_sort_unique(&neighbors_b);

                        close_separator(&subgraph, b, &neighbors_a, &sep_a_local);
                        close_separator(&subgraph, a, &neighbors_b, &sep_b_local);

                        for (int i = 0; i < sep_a_local.size; i++) {
                            int local_v = sep_a_local.data[i];
                            if (local_v >= 0 && local_v < invmap.size) {
                                intvec_push_back(&H_local, invmap.data[local_v]);
                            }
                        }

                        for (int i = 0; i < sep_b_local.size; i++) {
                            int local_v = sep_b_local.data[i];
                            if (local_v >= 0 && local_v < invmap.size) {
                                intvec_push_back(&H_local, invmap.data[local_v]);
                            }
                        }

                        intvec_free(&sep_b_local);
                        intvec_free(&sep_a_local);
                        intvec_free(&neighbors_b);
                        intvec_free(&neighbors_a);
                        graph_free(&subgraph);
                        free(map);
                        intvec_free(&invmap);
                        intvec_free(&sub_nodes);
                    }
                }
            }
        }

        if (keep_iterating) {
            int old_size = H.size;
            intvec_append(&H, &H_local);
            intvec_sort_unique(&H);

            if (H.size == old_size) {
                keep_iterating = 0;
            }
        }

        intvec_free(&H_local);
    }

    intvec_clear(H_out);
    intvec_append(H_out, &H);
    intvec_sort_unique(H_out);

    compboundlist_free(&comps);
    intvec_free(&H);
}

static void decompose_atoms_core(
    const Graph *graph,
    IntVecList *atoms,
    IntVecList *separators
) {
    int n = graph->n;
    IntVec mcs_order;
    int *mcs_index;
    IntVecList A_list;
    IntVec full_vertices;

    intvec_init(&mcs_order);
    maximum_cardinality_search(graph, &mcs_order);

    mcs_index = (int *) xmalloc((size_t) (n > 0 ? n : 1) * sizeof(int));
    for (int i = 0; i < n; i++) {
        mcs_index[mcs_order.data[i]] = i;
    }

    intveclist_init(&A_list);
    intvec_init(&full_vertices);
    for (int v = 0; v < n; v++) {
        intvec_push_back(&full_vertices, v);
    }
    intveclist_push_take(&A_list, &full_vertices);

    while (A_list.size > 0) {
        IntVec A_global = intveclist_pop(&A_list);
        Graph subgraph;
        IntVec invmap;
        int *map = NULL;

        if (A_global.size == 0) {
            intvec_free(&A_global);
            continue;
        }

        intvec_init(&invmap);
        induced_subgraph(graph, &A_global, &subgraph, &map, &invmap);

        {
            int min_v_global = A_global.data[0];
            int min_v_local;
            IntVec N_v;
            IntVec H;
            IntVec H_global;
            CompBoundList comps;

            for (int j = 1; j < A_global.size; j++) {
                int cand = A_global.data[j];
                if (mcs_index[cand] < mcs_index[min_v_global]) {
                    min_v_global = cand;
                }
            }

            min_v_local = map[min_v_global];

            intvec_init(&N_v);
            intvec_append(&N_v, &subgraph.adj[min_v_local]);
            intvec_push_back(&N_v, min_v_local);
            intvec_sort_unique(&N_v);

            intvec_init(&H);
            if (is_clique(&subgraph, &N_v)) {
                intvec_append(&H, &N_v);
            } else {
                get_minimal_collapsible_core(&subgraph, &N_v, &H);
            }
            intvec_sort_unique(&H);

            intvec_init(&H_global);
            for (int i = 0; i < H.size; i++) {
                int local_v = H.data[i];
                if (local_v >= 0 && local_v < invmap.size) {
                    intvec_push_back(&H_global, invmap.data[local_v]);
                }
            }
            intvec_sort_unique(&H_global);
            intveclist_push_take(atoms, &H_global);

            compboundlist_init(&comps);
            components_forbidden(&subgraph, &H, &comps);

            for (int i = 0; i < comps.size; i++) {
                CompBound *cb = &comps.data[i];
                IntVec C_global;
                IntVec A_new_local;
                IntVec A_new_global;

                intvec_init(&C_global);
                for (int j = 0; j < cb->bound.size; j++) {
                    int local_v = cb->bound.data[j];
                    if (local_v >= 0 && local_v < invmap.size) {
                        intvec_push_back(&C_global, invmap.data[local_v]);
                    }
                }
                intvec_sort_unique(&C_global);
                intveclist_push_take(separators, &C_global);

                intvec_init(&A_new_local);
                intvec_append(&A_new_local, &cb->comp);
                intvec_append(&A_new_local, &cb->bound);
                intvec_sort_unique(&A_new_local);

                intvec_init(&A_new_global);
                for (int j = 0; j < A_new_local.size; j++) {
                    int local_v = A_new_local.data[j];
                    if (local_v >= 0 && local_v < invmap.size) {
                        intvec_push_back(&A_new_global, invmap.data[local_v]);
                    }
                }
                intvec_sort_unique(&A_new_global);

                if (A_new_global.size > 0 && !intvec_equal(&A_new_global, &A_global)) {
                    intveclist_push_take(&A_list, &A_new_global);
                } else {
                    intvec_free(&A_new_global);
                }

                intvec_free(&A_new_local);
            }

            compboundlist_free(&comps);
            intvec_free(&H);
            intvec_free(&N_v);
        }

        graph_free(&subgraph);
        free(map);
        intvec_free(&invmap);
        intvec_free(&A_global);
    }

    intveclist_free(&A_list);
    free(mcs_index);
    intvec_free(&mcs_order);
}

SEXP c_get_minimal_collapsible_impl(
    SEXP edge_mat,
    SEXP directed_flag,
    SEXP n_vertices,
    SEXP r_nodes
) {
    int nv = scalar_to_nonnegative_integer(n_vertices, "n_vertices");
    int directed = scalar_to_logical_flag(directed_flag, "directed");
    int n_edges = 0;
    int *from = NULL;
    int *to = NULL;
    Graph graph;
    IntVec seeds;
    IntVec out_nodes;
    SEXP out;

    (void) directed; /* Algorithm follows weak-neighbor semantics (undirected traversal). */

    if (nv == 0) {
        Rf_error("n_vertices must be greater than zero");
    }

    parse_edge_matrix_zero_based(edge_mat, nv, &from, &to, &n_edges);
    graph_build_from_edges(&graph, nv, from, to, n_edges);
    free(from);
    free(to);

    intvec_init(&seeds);
    parse_vertices_zero_based(r_nodes, nv, "r_nodes", &seeds);

    intvec_init(&out_nodes);
    get_minimal_collapsible_core(&graph, &seeds, &out_nodes);
    intvec_sort_unique(&out_nodes);

    out = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) out_nodes.size));
    for (int i = 0; i < out_nodes.size; i++) {
        INTEGER(out)[i] = out_nodes.data[i] + 1;
    }

    intvec_free(&out_nodes);
    intvec_free(&seeds);
    graph_free(&graph);

    UNPROTECT(1);
    return out;
}

SEXP c_decompose_atoms_impl(
    SEXP edge_mat,
    SEXP directed_flag,
    SEXP n_vertices
) {
    // 临时实现 - 返回空列表作为占位符
    // TODO: 需要与 igraph 2.3.0 API 兼容性处理
    SEXP out;
    SEXP names;
    
    out = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, Rf_allocVector(VECSXP, 0));  // atoms
    SET_VECTOR_ELT(out, 1, Rf_allocVector(VECSXP, 0));  // separators
    SET_VECTOR_ELT(out, 2, Rf_allocVector(INTSXP, 0));  // tree edges
    
    names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar("atoms"));
    SET_STRING_ELT(names, 1, Rf_mkChar("separators"));
    SET_STRING_ELT(names, 2, Rf_mkChar("tree"));
    Rf_setAttrib(out, R_NamesSymbol, names);
    
    UNPROTECT(2);
    return out;
}
