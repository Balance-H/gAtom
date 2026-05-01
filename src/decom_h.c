//#define _CRTDBG_MAP_ALLOC //用于内存泄露检测
//#include <crtdbg.h>

#include <igraph.h>
#include <stdlib.h>
#include <string.h>  // for memset
#include <stdbool.h>
#include <omp.h>
#include "uthash.h"
#include <stdio.h>
#include <time.h>
#include "decom_h.h"

/**
 * @brief 计算图中带有禁忌节点约束的弱连通分量及其边界节点集合。
 * 
 * 本函数在给定图中查找所有的弱连通分量，遍历过程中跳过指定的禁忌节点。
 * 每个连通分量由不包含禁忌节点的顶点集合构成，同时返回该连通分量与禁忌节点相邻的边界节点集合。
 * 禁忌节点既不会作为分量成员，也不会参与遍历扩展，但会被记录为相应连通分量的边界节点。
 * 
 * @param graph 输入的图结构，类型为 igraph_t*，要求图已初始化且合法。
 * @param components 已初始化的 igraph_vector_ptr_t* 指针，用于输出各连通分量，每个元素为指向 igraph_vector_int_t 的指针，存储该分量的节点集合。
 * @param boundaries 已初始化的 igraph_vector_ptr_t* 指针，用于输出对应连通分量的边界节点集合，每个元素为指向 igraph_vector_int_t 的指针。
 * @param forbidden_vertices 禁忌节点集合，类型为 igraph_vector_int_t*，指定的节点不会被访问或包含于分量。
 * 
 * @return 返回状态码，IGRAPH_SUCCESS 表示成功，其他值表示对应错误（如内存不足、无效顶点等）。
 * 
 * @note
 * - 禁忌节点在遍历过程中被标记以避免重复访问。
 * - 函数内部通过广度优先搜索(BFS)实现连通分量检测。
 * - 调用者需确保传入的 components 和 boundaries 已正确初始化，且调用后负责释放其中分配的内存。
 */


igraph_error_t components_forbidden(
    const igraph_t *graph,
    igraph_vector_ptr_t *components,
    igraph_vector_ptr_t *boundaries,
    const igraph_vector_int_t *forbidden_vertices)
{
    igraph_integer_t n = igraph_vcount(graph);

    // 访问标记数组：0未访问，1已访问, 2是禁忌节点
    char *visited = (char*)calloc(n, sizeof(char));
    if (!visited) {
        return IGRAPH_ENOMEM;
    }

    // 标记禁忌节点为已访问，防止遍历
    if (forbidden_vertices) {
        for (igraph_integer_t i = 0; i < igraph_vector_int_size(forbidden_vertices); i++) {
            igraph_integer_t v = VECTOR(*forbidden_vertices)[i];
            visited[v] = 2;
        }
    }

     // 边界节点辅助标记，避免重复加入
    char *boundary_marked = (char*)calloc(n, sizeof(char));
    if (!boundary_marked) {
        free(visited);
        return IGRAPH_ENOMEM;
    }

    // 用于 BFS 的队列（动态数组）
    igraph_integer_t *queue = (igraph_integer_t*)malloc(n * sizeof(igraph_integer_t));
     if (!queue) {
        free(visited);
        free(boundary_marked);
        return IGRAPH_ENOMEM;
    }

    igraph_vector_int_t neighbors;
    igraph_vector_int_t parents;
    igraph_vector_int_init(&neighbors, 0);
    igraph_vector_int_init(&parents, 0);

    for (igraph_integer_t v = 0; v < n; v++) {
        if (visited[v]) continue;  // 已访问或禁忌跳过

        // 新连通分量和边界初始化
        igraph_vector_int_t *comp = (igraph_vector_int_t*)malloc(sizeof(igraph_vector_int_t));
        igraph_vector_int_t *bound = (igraph_vector_int_t*)malloc(sizeof(igraph_vector_int_t));
        if (!comp || !bound) {
            free(visited);
            free(boundary_marked);
            free(queue);
            igraph_vector_int_destroy(&neighbors);
            if (comp) free(comp);
            if (bound) free(bound);
            return IGRAPH_ENOMEM;
        }
        
        igraph_vector_int_init(comp, 0);
        igraph_vector_int_init(bound, 0);

        // BFS 初始化
        igraph_integer_t head = 0, tail = 0;
        queue[tail++] = v;
        visited[v] = 1;

        while (head < tail) {
            igraph_integer_t cur = queue[head++];

            igraph_vector_int_push_back(comp, cur);

            IGRAPH_CHECK(igraph_neighbors(graph, &neighbors, cur, IGRAPH_ALL));
            igraph_integer_t nei_count = igraph_vector_int_size(&neighbors);

            for (igraph_integer_t i = 0; i < nei_count; i++) {
                igraph_integer_t w = VECTOR(neighbors)[i];
                if (visited[w] == 2) {
                    // 禁忌节点，加入边界（避免重复）
                    if (!boundary_marked[w]) {
                        igraph_vector_int_push_back(bound, w);
                        boundary_marked[w] = 1;
                    }
                    continue; // 不加入队列
                }

                if (!visited[w]) {
                    visited[w] = 1;
                    queue[tail++] = w;
                }
            }
        }

        igraph_vector_ptr_push_back(components, comp);
        igraph_vector_ptr_push_back(boundaries, bound);

        // 清空边界辅助标记，为下一个连通分量做准备
        igraph_integer_t bsize = igraph_vector_int_size(bound);
        for (igraph_integer_t i = 0; i < bsize; i++) {
            boundary_marked[VECTOR(*bound)[i]] = 0;
        }

        
    }

    free(visited);
    free(boundary_marked);
    free(queue);
    igraph_vector_int_destroy(&neighbors);

    return IGRAPH_SUCCESS;
}

// 工具函数：unique + sort, 去重
igraph_integer_t vector_int_unique(igraph_vector_int_t *v) {
    if (igraph_vector_int_size(v) == 0) return IGRAPH_SUCCESS;
    igraph_vector_int_sort(v);
    igraph_integer_t  write_pos = 1;
    for (igraph_integer_t  i = 1; i < igraph_vector_int_size(v); i++) {
        if (VECTOR(*v)[i] != VECTOR(*v)[i - 1]) {
            VECTOR(*v)[write_pos++] = VECTOR(*v)[i];
        }
    }
    igraph_vector_int_resize(v, write_pos);
    return IGRAPH_SUCCESS;
}

// 返回 1 如果 S ⊆ U (按元素比较)，否则 0
static igraph_bool_t vector_int_is_subset(const igraph_vector_int_t *S, const igraph_vector_int_t *U) {
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(S); i++) {
        igraph_integer_t sv = VECTOR(*S)[i];
        igraph_bool_t found = 0;
        for (igraph_integer_t j = 0; j < igraph_vector_int_size(U); j++) {
            if (VECTOR(*U)[j] == sv) { found = 1; break; }
        }
        if (!found) return 0;
    }
    return 1;
}

/**
 * @brief 计算从指定起点出发、忽略禁忌节点的单点连通分量及其与禁忌节点相邻的边界节点。
 * 
 * 函数以指定顶点 `vertex` 作为起点，执行广度优先搜索（BFS）遍历图中除禁忌节点外的可达节点，
 * 构成该起点所在的连通分量。同时收集所有与该连通分量邻接但被标记为禁忌的节点，作为边界节点返回。
 * 
 * @param graph 输入图，类型为 `const igraph_t*`，要求已初始化且合法。
 * @param vertex 起始节点编号，类型为 `igraph_integer_t`，必须在图节点范围内。
 * @param forbidden_vertices 禁忌节点集合，类型为 `const igraph_vector_int_t*`，遍历时忽略这些节点。
 * @param bound_b 输出向量，类型为 `igraph_vector_int_t*`，返回遇到的禁忌节点集合（即边界节点）。
 * 
 * @return 返回状态码，`IGRAPH_SUCCESS` 表示执行成功，其他值代表错误状态（如起点无效、内存不足等）。
 * 
 * @note
 * - 禁忌节点在遍历过程中被标记为不可访问，且不包含于连通分量。
 * - 边界节点即为直接邻接于连通分量但为禁忌节点的集合。
 * - 函数内部使用队列实现 BFS，效率较高，适合中大型图遍历。
 * - 调用者需确保输入参数有效，且负责管理 `bound_b` 的内存。
 */


igraph_error_t close_separator(
    const igraph_t *graph,
    igraph_integer_t vertex,  //节点b,算法会序
    const igraph_vector_int_t *forbidden_vertices,
    igraph_vector_int_t *bound_b  // 新增的参数，用于存储遇到的禁忌节点
) {
    igraph_integer_t n = igraph_vcount(graph);
    if (vertex < 0 || vertex >= n) {
        return IGRAPH_EINVVID;
    }

    // 标记数组：0未访问，1已访问或禁忌
    char *visited = calloc(n, sizeof(char));
    if (!visited) {
        return IGRAPH_ENOMEM;
    }

    // 标记禁忌节点为已访问，防止遍历
    if (forbidden_vertices) {
        for (igraph_integer_t i = 0; i < igraph_vector_int_size(forbidden_vertices); i++) {
            igraph_integer_t v = VECTOR(*forbidden_vertices)[i];
            if (v >= 0 && v < n) {
                visited[v] = 2;
            }
        }
    }

    igraph_vector_int_t res;
    igraph_vector_int_init(&res, 0);
    igraph_vector_int_clear(bound_b);

    // 用动态数组做队列
    igraph_integer_t *queue = malloc(n * sizeof(igraph_integer_t));
    if (!queue) {
        free(visited);
        return IGRAPH_ENOMEM;
    }
    igraph_integer_t head = 0, tail = 0;

    // 初始化队列
    queue[tail++] = vertex;
    visited[vertex] = 1;

    igraph_vector_int_t neighbors;
    igraph_vector_int_t parents;
    igraph_vector_int_init(&neighbors, 0);
    igraph_vector_int_init(&parents, 0);

    while (head < tail) {
        igraph_integer_t cur = queue[head++];

        igraph_vector_int_push_back(&res, cur);

        IGRAPH_CHECK(igraph_neighbors(graph, &neighbors, cur, IGRAPH_ALL));
        igraph_integer_t nei_count = igraph_vector_int_size(&neighbors);
        for (igraph_integer_t i = 0; i < nei_count; i++) {
            igraph_integer_t w = VECTOR(neighbors)[i];

            // 如果是禁忌节点 (visited[w] == 2)，将其加入 bound_b，并标记为已处理
            if (visited[w] == 2) {
                igraph_vector_int_push_back(bound_b, w);  // 将禁忌节点 w 加入 bound_b
                visited[w] = 1;  // 标记禁忌节点为已处理
                continue;  // 跳过继续处理该节点，直接进入下一个邻居
            }

            if (!visited[w]) {
                visited[w] = 1;
                queue[tail++] = w;
            }
        }
    }

    igraph_vector_int_destroy(&neighbors);
    igraph_vector_int_destroy(&res);
    free(queue);
    free(visited);

    return IGRAPH_SUCCESS;
}


/**
 * get_minimal_collapsible - 基于节点集 r_nodes 对图 graph 进行 CMSA（Clique Minimal Separator Algorithm）凸包扩展，
 *              并输出扩展后的节点集 H_out。
 *
 * 输入参数：
 *   - graph: 指向已初始化的 igraph_t 图结构的指针，表示待处理的图。
 *   - r_nodes: 指向 igraph_vector_int_t 类型的向量，包含初始的节点集合（即起始的禁忌节点集合）。
 *
 * 输出参数：
 *   - H_out: 指向 igraph_vector_int_t 类型的向量，用于返回扩展后的节点集合结果。
 *
 * 功能说明：
 *   本函数实现对输入图的递归处理，通过计算禁忌节点集 H 中节点的弱连通分量及其边界，
 *   在边界节点中寻找非邻接节点对，并通过调用 close_separator 函数扩展禁忌节点集 H。
 *   该过程重复进行，直到不能再找到新的非邻接节点对，停止迭代。
 *
 * 具体步骤：
 *   1. 初始化节点集合 H 为 r_nodes 的拷贝。
 *   2. 在迭代中，调用 components_forbidden 函数计算带有禁忌节点的连通分量及边界节点。
 *   3. 遍历每个连通分量的边界节点集合，寻找非邻接节点对。
 *   4. 对于每对非邻接节点，调用 close_separator 函数计算最小分割点集合，并将结果加入节点集合 H。
 *   5. 当所有连通分量的边界节点均无非邻接节点对时，迭代终止。
 *
 * 返回值：
 *   - IGRAPH_SUCCESS 表示函数执行成功。
 *   - 其他返回码表示执行过程中发生错误。
 *
 * 注意事项：
 *   - 本函数在迭代过程中会频繁分配和释放内存，可能影响性能，建议结合预分配和向量重用以优化。
 *   - 输入图应已初始化且有效，r_nodes 中的节点应均在图顶点范围内。
 *   - 输出向量 H_out 会被清空并重新初始化。
 */


igraph_error_t get_minimal_collapsible(
    const igraph_t *graph,
    const igraph_vector_int_t *r_nodes,
    igraph_vector_int_t *H_out
) {
    igraph_vector_int_t H;                 // 当前节点集合 H
    igraph_vector_ptr_t components;       // 连通分量列表
    igraph_vector_ptr_t boundaries;       // 边界节点列表
    igraph_bool_t s = 1;                   // 是否继续迭代标志
    igraph_error_t ret = IGRAPH_SUCCESS;

    igraph_integer_t n = igraph_vcount(graph);

    // 初始化 H，预分配空间
    //igraph_vector_int_init(&H,0);
    //igraph_vector_int_reserve(&H, n);
    //igraph_vector_int_copy(&H, r_nodes);

    igraph_vector_int_init_copy(&H, r_nodes);

    

    igraph_vector_ptr_init(&components, 0);
    igraph_vector_ptr_init(&boundaries, 0);

    // 预分配循环内常用向量，单线程分配
    igraph_vector_int_t sub_nodes;
    igraph_vector_int_t neighbors_a, neighbors_b;
    igraph_vector_int_t sep_a_local, sep_b_local;

    igraph_vector_int_init(&sub_nodes,0);
    igraph_vector_int_reserve(&sub_nodes, n);

    igraph_vector_int_init(&neighbors_a,0);
    igraph_vector_int_reserve(&neighbors_a, n);

    igraph_vector_int_init(&neighbors_b,0);
    igraph_vector_int_reserve(&neighbors_b, n);

    igraph_vector_int_init(&sep_a_local,0);
    igraph_vector_int_reserve(&sep_a_local, n);

    igraph_vector_int_init(&sep_b_local,0);
    igraph_vector_int_reserve(&sep_b_local, n);

    igraph_vector_int_t map, invmap;
    igraph_vector_int_init(&map, 0);
    igraph_vector_int_init(&invmap, 0);

    while (s) {
        s = 0;

        // 清理上次结果
        for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&components); i++) {
            igraph_vector_int_t *comp = (igraph_vector_int_t*) VECTOR(components)[i];
            igraph_vector_int_destroy(comp);
            free(comp);
        }
        for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&boundaries); i++) {
            igraph_vector_int_t *bound = (igraph_vector_int_t*) VECTOR(boundaries)[i];
            igraph_vector_int_destroy(bound);
            free(bound);
        }

        igraph_vector_ptr_clear(&components);
        igraph_vector_ptr_clear(&boundaries);

        components_forbidden(graph, &components, &boundaries, &H);

        igraph_integer_t ncomp = igraph_vector_ptr_size(&components);

        // 开启并行循环，每个线程维护局部H_local，循环结束后合并
        #pragma omp parallel
        {
            igraph_vector_int_t H_local;
            igraph_vector_int_init(&H_local,0);
            igraph_vector_int_reserve(&H_local, n);

            #pragma omp for schedule(dynamic)
            for (igraph_integer_t cidx = 0; cidx < ncomp; cidx++) {
                igraph_vector_int_t *component = (igraph_vector_int_t*) VECTOR(components)[cidx];
                igraph_vector_int_t *boundary = (igraph_vector_int_t*) VECTOR(boundaries)[cidx];

                if (igraph_vector_int_size(boundary) < 2) {
                    continue;
                }

                igraph_vector_int_sort(boundary);
                igraph_integer_t bsize = igraph_vector_int_size(boundary);

                igraph_bool_t found = 0;

                for (igraph_integer_t j = 0; j < bsize && !found; j++) {
                    igraph_integer_t hj = VECTOR(*boundary)[j];
                    for (igraph_integer_t k = j + 1; k < bsize && !found; k++) {
                        igraph_integer_t hk = VECTOR(*boundary)[k];

                        igraph_bool_t are_adj = 0;
                        igraph_are_adjacent(graph, hj, hk, &are_adj);
                        if (!are_adj) {
                            found = 1;
                            s = 1; // 标记继续迭代
                            // 清空并复用sub_nodes
                            igraph_vector_int_clear(&sub_nodes);
                            for (igraph_integer_t ci = 0; ci < igraph_vector_int_size(component); ci++) {
                                igraph_vector_int_push_back(&sub_nodes, VECTOR(*component)[ci]);
                            }
                            igraph_vector_int_push_back(&sub_nodes, hj);
                            igraph_vector_int_push_back(&sub_nodes, hk);
                            //igraph_vector_int_sort(&sub_nodes);

                            igraph_t subgraph; 
                            igraph_vector_int_init(&map, 0);
                            igraph_vector_int_init(&invmap, 0);

                            // 创建子图
                            igraph_vs_t vs;
                            igraph_vs_vector(&vs, &sub_nodes);
                            IGRAPH_CHECK(igraph_induced_subgraph_map(graph, &subgraph, vs,
                                                                        IGRAPH_SUBGRAPH_AUTO, &map, &invmap));
                            igraph_vs_destroy(&vs);
                            // 通过 map 获取局部编号
                            igraph_integer_t a = VECTOR(map)[hj] - 1;
                            igraph_integer_t b = VECTOR(map)[hk] - 1;

                            // 固定采用 CMSA 路线：调用两次 close_separator
                            igraph_vector_int_clear(&neighbors_a);
                            igraph_vector_int_clear(&neighbors_b);
                            igraph_neighbors(&subgraph, &neighbors_a, a, IGRAPH_ALL);
                            igraph_neighbors(&subgraph, &neighbors_b, b, IGRAPH_ALL);


                            igraph_vector_int_clear(&sep_a_local);
                            igraph_vector_int_clear(&sep_b_local);

                            close_separator(&subgraph, b, &neighbors_a, &sep_a_local);
                            close_separator(&subgraph, a, &neighbors_b, &sep_b_local);


                            for (igraph_integer_t si = 0; si < igraph_vector_int_size(&sep_a_local); si++) {
                                igraph_integer_t gv = VECTOR(invmap)[VECTOR(sep_a_local)[si]];
                                igraph_vector_int_push_back(&H_local, gv);
                            }
                            for (igraph_integer_t si = 0; si < igraph_vector_int_size(&sep_b_local); si++) {
                                igraph_integer_t gv = VECTOR(invmap)[VECTOR(sep_b_local)[si]];
                                igraph_vector_int_push_back(&H_local, gv);
                            }
                            igraph_destroy(&subgraph);
                            igraph_vector_int_destroy(&map);
                            igraph_vector_int_destroy(&invmap);
                            
                        }
                    }
                }
            }
            // 并行区域结束前合并线程局部结果到全局H，使用临界区保护
            #pragma omp critical
            {
                for (igraph_integer_t i = 0; i < igraph_vector_int_size(&H_local); i++) {
                    igraph_vector_int_push_back(&H, VECTOR(H_local)[i]);
                }
            }

            igraph_vector_int_destroy(&H_local);
        } // omp parallel end

        // 归一化H，去重
        vector_int_unique(&H);
    }

    igraph_vector_int_clear(H_out);
    igraph_vector_int_init_copy(H_out, &H);

    // 释放组件和边界内存
    for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&components); i++) {
        igraph_vector_int_t *comp = (igraph_vector_int_t*) VECTOR(components)[i];
        igraph_vector_int_destroy(comp);
        free(comp);
    }
    for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&boundaries); i++) {
        igraph_vector_int_t *bound = (igraph_vector_int_t*) VECTOR(boundaries)[i];
        igraph_vector_int_destroy(bound);
        free(bound);
    }

    igraph_vector_ptr_destroy(&components);
    igraph_vector_ptr_destroy(&boundaries);

    igraph_vector_int_destroy(&H);
    igraph_vector_int_destroy(&map);
    igraph_vector_int_destroy(&invmap);

    igraph_vector_int_destroy(&sub_nodes);
    igraph_vector_int_destroy(&neighbors_a);
    igraph_vector_int_destroy(&neighbors_b);
    igraph_vector_int_destroy(&sep_a_local);
    igraph_vector_int_destroy(&sep_b_local);

    return ret;
}


/*
 * decompose_atoms - 对无向图进行原子分解，并构建原子分解树。
 *
 * 本函数首先对输入无向图 g 使用最大基数搜索（MCS）计算节点排序，
 * 然后基于该排序迭代进行递归分解。每次从当前子图中选取最小 alpha 值节点 v，
 * 对其闭包 N[G'][v] 进行 atom 扩展，生成一个 atom U 并将其加入 atoms 列表。
 * 对于当前任务队列中的每个pair (S, U_p)，若 S ⊆ U，则在 tree_out 中
 * 添加以 U 为节点的新边，建立 atom tree 的父子关系。
 *
 * 接着，计算 G' - U 的连通分量，每个分量 M 对应的分离集 S_M = N_{G'}(M)，
 * 将 S_M 加入 separators，同时构建该分量的新任务，继续递归分解。
 *
 * 参数:
 *   - g: 输入无向图，必须是无向图，否则返回 IGRAPH_EINVAL。
 *   - atoms: 输出原子集合向量指针，每个元素为一个原子顶点集合（igraph_vector_int_t *）。
 *   - separators: 输出分离子集向量指针，每个元素为一个分离集顶点集合。
 *   - tree_out: 输出 atom tree 的无向边组合，节点对应 atoms 中的 atom 索引，边表示连接关系。
 *
 * 返回:
 *   - IGRAPH_SUCCESS 表示成功。
 *   - 其他错误码表示内部计算或内存分配失败。
 */

igraph_error_t decompose_atoms(
    const igraph_t *g,
    igraph_vector_ptr_t *atoms, // 存储 Atoms (A)
    igraph_vector_ptr_t *separators, // 存储 Separators (C)
    igraph_t *tree_out // 存储 atom tree
) {

        // Step 1: MCS
    igraph_vector_int_t mcs_alpha, mcs_order;
    igraph_vector_int_init(&mcs_alpha, 0);
    igraph_vector_int_init(&mcs_order, 0);
    igraph_maximum_cardinality_search(g, &mcs_alpha, &mcs_order);


    igraph_vector_int_destroy(&mcs_alpha);

    igraph_vector_int_t mcs_index;
    igraph_vector_int_init(&mcs_index, igraph_vcount(g));
    for (igraph_integer_t i = 0; i < igraph_vcount(g); i++) {
        VECTOR(mcs_index)[VECTOR(mcs_order)[i]] = i;
    }
    igraph_vector_int_destroy(&mcs_order);


    igraph_integer_t n = igraph_vcount(g);

    igraph_vector_ptr_t A_list;
    igraph_vector_ptr_init(&A_list, 0);

    igraph_vector_int_t *Vfull = (igraph_vector_int_t *)malloc(sizeof(igraph_vector_int_t));
    igraph_vector_int_init_range(Vfull, 0, n); 

    /* initial empty W for the whole graph */
    igraph_vector_ptr_t *initial_W = (igraph_vector_ptr_t*)malloc(sizeof(igraph_vector_ptr_t));
    igraph_vector_ptr_init(initial_W, 0);

    /* wrap into a task_t and push onto A_list */
    task_t *initial_task = (task_t*)malloc(sizeof(task_t));
    initial_task->V_nodes = Vfull;
    initial_task->W = initial_W;
    igraph_vector_ptr_push_back(&A_list, initial_task);


    igraph_vector_int_t map, invmap;
    igraph_vector_int_init(&map, 0);
    igraph_vector_int_init(&invmap, 0);

    while (!igraph_vector_ptr_empty(&A_list)) {   
        task_t *task = (task_t*)VECTOR(A_list)[igraph_vector_ptr_size(&A_list) - 1];
        igraph_vector_ptr_pop_back(&A_list);

        igraph_vector_int_t *A_global = task->V_nodes;
        igraph_vector_ptr_t *W_current = task->W;


        igraph_t sub_g; 
        igraph_vector_int_clear(&map); 
        igraph_vector_int_clear(&invmap);

        igraph_vs_t vs;
        igraph_vs_vector(&vs, A_global); 
        IGRAPH_CHECK(igraph_induced_subgraph_map(g, &sub_g, vs, IGRAPH_SUBGRAPH_AUTO, &map, &invmap));
        igraph_vs_destroy(&vs);


        igraph_integer_t min_v_global = VECTOR(*A_global)[0];
        for (igraph_integer_t j = 1; j < igraph_vector_int_size(A_global); j++) {
            igraph_integer_t cand = VECTOR(*A_global)[j];
            if (VECTOR(mcs_index)[cand] < VECTOR(mcs_index)[min_v_global]) {
                min_v_global = cand; 
            }
        }
        igraph_integer_t min_v_local = VECTOR(map)[min_v_global]-1;


        igraph_vector_int_t N_v;
        igraph_vector_int_init(&N_v, 0);
        igraph_neighbors(&sub_g, &N_v, min_v_local, IGRAPH_ALL);
        igraph_vector_int_push_back(&N_v, min_v_local);


        igraph_bool_t is_clique;
        igraph_vs_vector(&vs, &N_v);
        igraph_is_clique(&sub_g, vs, 0, &is_clique);
        igraph_vs_destroy(&vs);


        igraph_vector_int_t H;
        igraph_vector_int_init(&H, 0);

        if (is_clique) {
            igraph_vector_int_append(&H, &N_v);
        } else {
            get_minimal_collapsible(&sub_g, &N_v, &H);
        }
        igraph_vector_int_destroy(&N_v);


        igraph_vector_int_t *H_global = malloc(sizeof(igraph_vector_int_t));
        igraph_vector_int_init(H_global, igraph_vector_int_size(&H));
        for (igraph_integer_t i = 0; i < igraph_vector_int_size(&H); i++) {
            VECTOR(*H_global)[i] = VECTOR(invmap)[VECTOR(H)[i]];
        }
        igraph_vector_ptr_push_back(atoms, H_global);


        /* 新产生的 atom 在 atom 集合中的索引（用于在 tree 中连接） */
        igraph_integer_t u_index = igraph_vector_ptr_size(atoms) - 1;
        /* 为 atom tree 增加对应顶点 */
        IGRAPH_CHECK(igraph_add_vertices(tree_out, 1, 0));

        /* 将 W_current 中 S ⊆ U 的 pairs 连接到当前新 atom（u_index），并从 W_current 中移除 */
        if (W_current) {
            igraph_integer_t wsize = igraph_vector_ptr_size(W_current);
            for (igraph_integer_t wi = 0; wi < wsize; wi++) {
                pair_t *pp = (pair_t*) VECTOR(*W_current)[wi];
                if (!pp) continue;
                if (vector_int_is_subset(&pp->S, H_global)) {
                    /* 添加边 u -> pp->u_p_index */
                    IGRAPH_CHECK(igraph_add_edge(tree_out, u_index, pp->u_p_index));
                    /* 释放该 pair 结构并从 W_current 中删除 */
                    igraph_vector_ptr_remove(W_current, wi);
                    igraph_vector_int_destroy(&pp->S);
                    free(pp);
                    wi--; wsize--;
                }
            }
        }


        igraph_vector_ptr_t comps, bounds;
        igraph_vector_ptr_init(&comps, 0);
        igraph_vector_ptr_init(&bounds, 0);

        components_forbidden(&sub_g, &comps, &bounds, &H);
        igraph_vector_int_destroy(&H);


        igraph_integer_t ncomps = igraph_vector_ptr_size(&comps);
        for (igraph_integer_t i = 0; i < ncomps; i++) {
            igraph_vector_int_t *comp = (igraph_vector_int_t *)VECTOR(comps)[i];
            igraph_vector_int_t *bound = (igraph_vector_int_t *)VECTOR(bounds)[i];

            // ************************************************
            // ** 核心修改：将 N_{G'}(M) 映射回全局索引并存储到 separators **
            // ************************************************
            igraph_vector_int_t *C_global = (igraph_vector_int_t *)malloc(sizeof(igraph_vector_int_t));
            if (!C_global) { /* 错误处理 */ }
            IGRAPH_CHECK(igraph_vector_int_init(C_global, igraph_vector_int_size(bound)));
            for (igraph_integer_t j = 0; j < igraph_vector_int_size(bound); j++) {
                VECTOR(*C_global)[j] = VECTOR(invmap)[VECTOR(*bound)[j]];
            }
            IGRAPH_CHECK(igraph_vector_ptr_push_back(separators, C_global)); // 存储 Separator C


            igraph_vector_int_t A;
            igraph_vector_int_init(&A, 0);
            igraph_vector_int_append(&A, comp);
            igraph_vector_int_append(&A, bound);

            igraph_vector_int_t *A_global_ptr = (igraph_vector_int_t *)malloc(sizeof(igraph_vector_int_t));
            igraph_vector_int_init(A_global_ptr, igraph_vector_int_size(&A));
            for (igraph_integer_t j = 0; j < igraph_vector_int_size(&A); j++) {
                VECTOR(*A_global_ptr)[j] = VECTOR(invmap)[VECTOR(A)[j]];
            }


            /* 为本组件 M 构建对应的 W_M 并创建 task 推入 A_list */
            igraph_vector_ptr_t *W_M = (igraph_vector_ptr_t*)malloc(sizeof(igraph_vector_ptr_t));
            igraph_vector_ptr_init(W_M, 0);

            /* 将 (S_M, u_index) 放入 W_M，S_M 使用 C_global（已映射为全局索引） */
            pair_t *p0 = (pair_t*)malloc(sizeof(pair_t));
            igraph_vector_int_init(&p0->S, igraph_vector_int_size(C_global));
            for (igraph_integer_t jj = 0; jj < igraph_vector_int_size(C_global); jj++) VECTOR(p0->S)[jj] = VECTOR(*C_global)[jj];
            p0->u_p_index = u_index;
            igraph_vector_ptr_push_back(W_M, p0);

            /* 将 W_current 中 S ⊆ (M ∪ S_M) 的 pairs 复制到 W_M（按伪代码分发） */
            if (W_current) {
                igraph_integer_t wsize2 = igraph_vector_ptr_size(W_current);
                for (igraph_integer_t wi2 = 0; wi2 < wsize2; wi2++) {
                    pair_t *pp2 = (pair_t*) VECTOR(*W_current)[wi2];
                    if (!pp2) continue;
                    if (vector_int_is_subset(&pp2->S, A_global_ptr)) {
                        pair_t *pp_copy = (pair_t*)malloc(sizeof(pair_t));
                        igraph_vector_int_init(&pp_copy->S, igraph_vector_int_size(&pp2->S));
                        for (igraph_integer_t kk = 0; kk < igraph_vector_int_size(&pp2->S); kk++) {
                            VECTOR(pp_copy->S)[kk] = VECTOR(pp2->S)[kk];
                        }
                        pp_copy->u_p_index = pp2->u_p_index;
                        igraph_vector_ptr_push_back(W_M, pp_copy);
                    }
                }
            }

            /* 将 W_M 包装为任务并推入 A_list */
            task_t *task_M = (task_t*)malloc(sizeof(task_t));
            task_M->V_nodes = A_global_ptr;
            task_M->W = W_M;
            igraph_vector_ptr_push_back(&A_list, task_M);
            igraph_vector_int_destroy(&A);                     
        }

        for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&comps); i++) {
            igraph_vector_int_destroy((igraph_vector_int_t *)VECTOR(comps)[i]);
            free(VECTOR(comps)[i]);
            igraph_vector_int_destroy((igraph_vector_int_t *)VECTOR(bounds)[i]);
            free(VECTOR(bounds)[i]);
        }
        igraph_vector_ptr_destroy(&comps);
        igraph_vector_ptr_destroy(&bounds);
        igraph_destroy(&sub_g);
        // 清理当前 A_global
        igraph_vector_int_destroy(A_global);
        free(A_global);

        /* 清理当前任务的 W_current */
        if (W_current) {
            for (igraph_integer_t wi3 = 0; wi3 < igraph_vector_ptr_size(W_current); wi3++) {
                pair_t *pp3 = (pair_t*) VECTOR(*W_current)[wi3];
                if (!pp3) continue;
                igraph_vector_int_destroy(&pp3->S);
                free(pp3);
            }
            igraph_vector_ptr_destroy(W_current);
            free(W_current);
        }

        /* 释放当前任务结构 */
        free(task);

    }
    // 循环结束，A_list 已空，map / invmap / mcs_index 清理
    igraph_vector_int_destroy(&mcs_index);
    igraph_vector_int_destroy(&map);
    igraph_vector_int_destroy(&invmap);

    return IGRAPH_SUCCESS;
}






igraph_error_t SAHR(igraph_t *g, igraph_integer_t *r, igraph_integer_t r_size, igraph_integer_t **local2global, igraph_integer_t *result_size) {

    igraph_t g_copy;
    IGRAPH_CHECK(igraph_copy(&g_copy, g));

    igraph_integer_t n = igraph_vcount(&g_copy);
    *local2global = (igraph_integer_t *) malloc(n * sizeof(int));
    if (!*local2global) { igraph_destroy(&g_copy); return IGRAPH_ENOMEM; }
    for (igraph_integer_t i = 0; i < n; i++) (*local2global)[i] = i;

    igraph_integer_t local2global_size = n;  // 维护当前 local2global 的有效长度

    // 初始化 M：不在 r 中的节点
    igraph_integer_t *M = (igraph_integer_t *) malloc(n * sizeof(int));
    if (!M) { free(*local2global); igraph_destroy(&g_copy); return IGRAPH_ENOMEM; }
    igraph_integer_t M_size = 0;
    for (igraph_integer_t i = 0; i < n; i++) {
        igraph_integer_t found = 0;
        for (igraph_integer_t j = 0; j < r_size; j++) if (i == r[j]) { found = 1; break; }
        if (!found) M[M_size++] = i;
    }

    igraph_vector_int_t N;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&N, 0);
    igraph_vector_int_t neighbors_vec;
    igraph_vector_int_init(&neighbors_vec, 0);
    igraph_vs_t vs;

    igraph_integer_t changed = 1;
    while (changed) {
        changed = 0;

        igraph_integer_t N_size = igraph_vector_int_size(&N);

        // 处理上一轮删除节点邻集
        for (igraph_integer_t idx = 0; idx < N_size; idx++) {
            igraph_integer_t node = VECTOR(N)[idx];

            if (node >= igraph_vcount(&g_copy)) continue;

            IGRAPH_CHECK(igraph_neighbors(&g_copy, &neighbors_vec, node, IGRAPH_ALL));

            igraph_t subgraph;
            igraph_vs_vector(&vs, &neighbors_vec);
            igraph_error_t ret = igraph_induced_subgraph(&g_copy, &subgraph, vs, IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH);
            if (ret != IGRAPH_SUCCESS) continue;

            igraph_integer_t n_sub = igraph_vcount(&subgraph);
            igraph_integer_t m_sub = igraph_ecount(&subgraph);
            igraph_integer_t complete_edges = n_sub * (n_sub - 1) / 2;

            if (m_sub == complete_edges) {
                igraph_vs_t vs_delete;
                igraph_vs_1(&vs_delete, node);
                IGRAPH_CHECK(igraph_delete_vertices(&g_copy, vs_delete));
                igraph_vs_destroy(&vs_delete);

                // 更新 M
                igraph_integer_t new_M_size = 0;
                for (igraph_integer_t j = 0; j < M_size; j++) {
                    if (M[j] == node) continue;
                    M[new_M_size++] = (M[j] > node) ? M[j] - 1 : M[j];
                }
                M_size = new_M_size;

                // 更新 local2global
                for (igraph_integer_t j = node; j < local2global_size - 1; j++)
                    (*local2global)[j] = (*local2global)[j + 1];
                local2global_size--;
                *result_size = local2global_size;

                // 更新 N
                igraph_vector_int_t new_N;
                igraph_vector_int_init(&new_N, 0);
                for (igraph_integer_t j = 0; j < igraph_vector_int_size(&neighbors_vec); j++) {
                    igraph_integer_t neighbor = VECTOR(neighbors_vec)[j];
                    if (neighbor > node) neighbor--;
                    for (igraph_integer_t k = 0; k < M_size; k++) {
                        if (M[k] == neighbor) { igraph_vector_int_push_back(&new_N, neighbor); break; }
                    }
                }
                igraph_vector_int_init(&N, igraph_vector_int_size(&new_N)); // 初始化 N 大小
                for (long i = 0; i < igraph_vector_int_size(&new_N); i++) {
                    VECTOR(N)[i] = VECTOR(new_N)[i];
                }
                igraph_vector_int_destroy(&new_N);

                igraph_destroy(&subgraph);
                changed = 1;
                break;
            }
            igraph_destroy(&subgraph);
        }
        if (changed) continue;

        // 遍历剩余 M
        for (igraph_integer_t i = 0; i < M_size; i++) {
            igraph_integer_t node = M[i];

            igraph_integer_t in_N = 0;
            for (igraph_integer_t k = 0; k < igraph_vector_int_size(&N); k++) {
                if (VECTOR(N)[k] == node) {
                    in_N = 1;
                    break;
                }
            }
            if (in_N) continue;

            if (node >= igraph_vcount(&g_copy)) continue;

            IGRAPH_CHECK(igraph_neighbors(&g_copy, &neighbors_vec, node, IGRAPH_ALL));
            igraph_t subgraph;
            igraph_vs_vector(&vs, &neighbors_vec);
            IGRAPH_CHECK(igraph_induced_subgraph(&g_copy, &subgraph, vs, IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH));

            igraph_integer_t n_sub = igraph_vcount(&subgraph);
            igraph_integer_t m_sub = igraph_ecount(&subgraph);
            igraph_integer_t complete_edges = n_sub * (n_sub - 1) / 2;

            if (m_sub == complete_edges) {
                igraph_vs_t vs_delete;
                igraph_vs_1(&vs_delete, node);
                IGRAPH_CHECK(igraph_delete_vertices(&g_copy, vs_delete));
                igraph_vs_destroy(&vs_delete);

                igraph_integer_t new_M_size = 0;
                for (igraph_integer_t j = 0; j < M_size; j++) {
                    if (M[j] == node) continue;
                    M[new_M_size++] = (M[j] > node) ? M[j] - 1 : M[j];
                }
                M_size = new_M_size;

                for (igraph_integer_t j = node; j < local2global_size - 1; j++)
                    (*local2global)[j] = (*local2global)[j + 1];
                local2global_size--;
                *result_size = local2global_size;

                igraph_vector_int_t new_N;
                igraph_vector_int_init(&new_N, 0);
                for (igraph_integer_t j = 0; j < igraph_vector_int_size(&neighbors_vec); j++) {
                    igraph_integer_t neighbor = VECTOR(neighbors_vec)[j];
                    if (neighbor > node) neighbor--;
                    for (igraph_integer_t k = 0; k < M_size; k++) {
                        if (M[k] == neighbor) { igraph_vector_int_push_back(&new_N, neighbor); break; }
                    }
                }
                igraph_vector_int_init(&N, igraph_vector_int_size(&new_N)); // 初始化 N 大小
                for (long i = 0; i < igraph_vector_int_size(&new_N); i++) {
                    VECTOR(N)[i] = VECTOR(new_N)[i];
                }
                igraph_vector_int_destroy(&new_N);

                igraph_destroy(&subgraph);
                changed = 1;
                break;
            }
            igraph_destroy(&subgraph);
        }
    }

    *result_size = local2global_size;

    free(M);
    igraph_vector_int_destroy(&neighbors_vec);
    igraph_vector_int_destroy(&N);
    igraph_destroy(&g_copy);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}





/*
 * 实现最大基数搜索（MCS）算法，返回节点的 MCS 序列和对应的编号。
 * BLAIR, J. R. and PEYTON, B. (1993). An introduction to chordal graphs and clique trees. In Graph theory
and sparse matrix computation 1–29. Springer.
 * 参数:
 *  - graph: 输入图的常量指针
 *  - alpha: 输出参数，存储每个节点在 MCS 序列中的编号
 *  - alpham1: 输出参数，存储 MCS 序列中每个位置对应的节点编号
 * 返回:
 *  - igraph_error_t 类型，标示函数执行状态
 */
igraph_error_t mcs_with_cliques(const igraph_t *graph,
                                      igraph_vector_int_t *alpha,
                                      igraph_vector_int_t *alpham1,
                                      igraph_vector_ptr_t *cliques) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t size;
    igraph_vector_int_t head, next, prev; /* doubly linked list with head */
    igraph_integer_t i, j, v, x, k, len, w, ws, nw, pw;
    igraph_adjlist_t adjlist;
    igraph_vector_int_t *neis;
    igraph_integer_t pre_card = 0, new_card = 0;
    igraph_vector_int_t current_clique, L_set;
    igraph_vector_int_t *clique_copy;
    igraph_bool_t cliques_enabled;

    cliques_enabled = (cliques != NULL);

    if (cliques_enabled) {
        igraph_vector_ptr_clear(cliques);
        IGRAPH_CHECK(igraph_vector_int_init(&L_set, 0));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &L_set);
        IGRAPH_CHECK(igraph_vector_int_init(&current_clique, 0));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &current_clique);
    }

    if (no_of_nodes == 0) {
        igraph_vector_int_clear(alpha);
        if (alpham1) {
            igraph_vector_int_clear(alpham1);
        }
        if (cliques_enabled) {
            IGRAPH_FINALLY_CLEAN(2);
            igraph_vector_int_destroy(&current_clique);
            igraph_vector_int_destroy(&L_set);
        }
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&size, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&head, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&next, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&prev, no_of_nodes);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    IGRAPH_CHECK(igraph_vector_int_resize(alpha, no_of_nodes));
    if (alpham1) {
        IGRAPH_CHECK(igraph_vector_int_resize(alpham1, no_of_nodes));
    }

    /***********************************************/
    /* for i in [0,n-1] -> set(i) := emptyset rof; */
    /***********************************************/

    /* nothing to do, 'head' contains all zeros */

    /*********************************************************/
    /* for v in vertices -> size(v):=0; add v to set(0) rof; */
    /*********************************************************/

    VECTOR(head)[0] = 1;
    for (v = 0; v < no_of_nodes; v++) {
        VECTOR(next)[v] = v + 2;
        VECTOR(prev)[v] = v;
    }
    VECTOR(next)[no_of_nodes - 1] = 0;
    /* size is already all zero */

    /***************/
    /* i:=n; j:=0; */
    /***************/

    i = no_of_nodes; j = 0;
    pre_card = 0;

    /**************/
    /* do i>=1 -> */
    /**************/

    while (i >= 1) {
        /********************************/
        /* v :=  delete any from set(j) */
        /********************************/

        v = VECTOR(head)[j] - 1;
        x = VECTOR(next)[v];
        VECTOR(head)[j] = x;
        if (x != 0) {
            VECTOR(prev)[x - 1] = 0;
        }

        /*************************************************/
        /* alpha(v) := i; alpham1(i) := v; size(v) := -1 */
        /*************************************************/

        VECTOR(*alpha)[v] = i - 1;
        if (alpham1) {
            VECTOR(*alpham1)[i - 1] = v;
        }
        VECTOR(size)[v] = -1;

        /* Calculate new_card = |ne(v) ∩ L_{i+1}| and manage cliques */
        if (cliques_enabled) {
            new_card = 0;
            neis = igraph_adjlist_get(&adjlist, v);
            len = igraph_vector_int_size(neis);
            for (k = 0; k < len; k++) {
                w = VECTOR(*neis)[k];
                if (VECTOR(size)[w] == -1) {
                    /* w is in L_set (already processed) */
                    new_card++;
                }
            }
            
            /* if new_card <= pre_card, start a new clique */
            if (new_card <= pre_card) {
                /* Save the previous clique if it's not empty */
                if (igraph_vector_int_size(&current_clique) > 0) {
                    /* Remove duplicates before saving */
                    igraph_vector_int_sort(&current_clique);
                    igraph_integer_t write_idx = 0;
                    for (int ciq_i = 0; ciq_i < igraph_vector_int_size(&current_clique); ciq_i++) {
                        if (ciq_i == 0 || VECTOR(current_clique)[ciq_i] != VECTOR(current_clique)[ciq_i - 1]) {
                            VECTOR(current_clique)[write_idx++] = VECTOR(current_clique)[ciq_i];
                        }
                    }
                    igraph_vector_int_resize(&current_clique, write_idx);
                    
                    clique_copy = (igraph_vector_int_t *)malloc(sizeof(igraph_vector_int_t));
                    if (!clique_copy) {
                        IGRAPH_ERROR("Memory allocation failed", IGRAPH_ENOMEM);
                    }
                    igraph_vector_int_init_copy(clique_copy, &current_clique);
                    igraph_vector_ptr_push_back(cliques, clique_copy);
                }
                /* Start new clique: clear and add neighbors in L_set */
                igraph_vector_int_clear(&current_clique);
            }
            
            /* Add neighbors in L_set to current clique */
            neis = igraph_adjlist_get(&adjlist, v);
            len = igraph_vector_int_size(neis);
            for (k = 0; k < len; k++) {
                w = VECTOR(*neis)[k];
                if (VECTOR(size)[w] == -1) {
                    /* w is in L_set */
                    igraph_vector_int_push_back(&current_clique, w);
                }
            }
            
            /* Add v itself to the clique */
            igraph_vector_int_push_back(&current_clique, v);
            
            /* Update pre_card */
            pre_card = new_card;
            
            /* Add v to L_set */
            igraph_vector_int_push_back(&L_set, v);
        }

        /********************************************/
        /* for {v,w} in E such that size(w) >= 0 -> */
        /********************************************/

        neis = igraph_adjlist_get(&adjlist, v);
        len = igraph_vector_int_size(neis);
        for (k = 0; k < len; k++) {
            w = VECTOR(*neis)[k];
            ws = VECTOR(size)[w];
            if (ws >= 0) {

                /******************************/
                /* delete w from set(size(w)) */
                /******************************/

                nw = VECTOR(next)[w];
                pw = VECTOR(prev)[w];
                if (nw != 0) {
                    VECTOR(prev)[nw - 1] = pw;
                }
                if (pw != 0) {
                    VECTOR(next)[pw - 1] = nw;
                } else {
                    VECTOR(head)[ws] = nw;
                }

                /******************************/
                /* size(w) := size(w)+1       */
                /******************************/

                VECTOR(size)[w] += 1;

                /******************************/
                /* add w to set(size(w))      */
                /******************************/

                ws = VECTOR(size)[w];
                nw = VECTOR(head)[ws];
                VECTOR(next)[w] = nw;
                VECTOR(prev)[w] = 0;
                if (nw != 0) {
                    VECTOR(prev)[nw - 1] = w + 1;
                }
                VECTOR(head)[ws] = w + 1;

            }
        }

        /***********************/
        /* i := i-1; j := j+1; */
        /***********************/

        i -= 1;
        j += 1;

        /*********************************************/
        /* do j>=0 and set(j)=emptyset -> j:=j-1; od */
        /*********************************************/

        if (j < no_of_nodes) {
            while (j >= 0 && VECTOR(head)[j] == 0) {
                j--;
            }
        }
    }

    /* Save the last clique if cliques output is enabled */
    if (cliques_enabled && igraph_vector_int_size(&current_clique) > 0) {
        /* Remove duplicates before saving final clique */
        igraph_vector_int_sort(&current_clique);
        igraph_integer_t write_idx = 0;
        for (int ciq_i = 0; ciq_i < igraph_vector_int_size(&current_clique); ciq_i++) {
            if (ciq_i == 0 || VECTOR(current_clique)[ciq_i] != VECTOR(current_clique)[ciq_i - 1]) {
                VECTOR(current_clique)[write_idx++] = VECTOR(current_clique)[ciq_i];
            }
        }
        igraph_vector_int_resize(&current_clique, write_idx);
        
        clique_copy = (igraph_vector_int_t *)malloc(sizeof(igraph_vector_int_t));
        if (!clique_copy) {
            IGRAPH_ERROR("Memory allocation failed", IGRAPH_ENOMEM);
        }
        igraph_vector_int_init_copy(clique_copy, &current_clique);
        igraph_vector_ptr_push_back(cliques, clique_copy);
    }

    igraph_adjlist_destroy(&adjlist);
    if (cliques_enabled) {
        igraph_vector_int_destroy(&current_clique);
        igraph_vector_int_destroy(&L_set);
        IGRAPH_FINALLY_CLEAN(2);
    }
    igraph_vector_int_destroy(&prev);
    igraph_vector_int_destroy(&next);
    igraph_vector_int_destroy(&head);
    igraph_vector_int_destroy(&size);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}




igraph_error_t dag_components_forbidden(
    const igraph_t *graph,
    igraph_vector_ptr_t *components,
    igraph_vector_ptr_t *boundaries,
    const igraph_vector_int_t *forbidden_vertices)
{
    igraph_integer_t n = igraph_vcount(graph);

    // 访问标记数组：0未访问，1已访问, 2是禁忌节点
    char *visited = (char*)calloc(n, sizeof(char));
    if (!visited) {
        return IGRAPH_ENOMEM;
    }

    // 标记禁忌节点为已访问，防止遍历
    if (forbidden_vertices) {
        for (igraph_integer_t i = 0; i < igraph_vector_int_size(forbidden_vertices); i++) {
            igraph_integer_t v = VECTOR(*forbidden_vertices)[i];
            visited[v] = 2;
        }
    }

     // 边界节点辅助标记，避免重复加入
    char *boundary_marked = (char*)calloc(n, sizeof(char));
    if (!boundary_marked) {
        free(visited);
        return IGRAPH_ENOMEM;
    }

    // 仅对“出边命中的禁忌节点 w”展开其父节点一次
    char *forbidden_out_expanded = (char*)calloc(n, sizeof(char));
    if (!forbidden_out_expanded) {
        free(visited);
        free(boundary_marked);
        return IGRAPH_ENOMEM;
    }

    // 用于 BFS 的队列（动态数组）
    igraph_integer_t *queue = (igraph_integer_t*)malloc(n * sizeof(igraph_integer_t));
     if (!queue) {
        free(visited);
        free(boundary_marked);
        free(forbidden_out_expanded);
        return IGRAPH_ENOMEM;
    }

    igraph_vector_int_t neighbors;
    igraph_vector_int_t parents;
    igraph_vector_int_init(&neighbors, 0);
    igraph_vector_int_init(&parents, 0);

    for (igraph_integer_t v = 0; v < n; v++) {
        if (visited[v]) continue;  // 已访问或禁忌跳过

        // 新连通分量和边界初始化
        igraph_vector_int_t *comp = (igraph_vector_int_t*)malloc(sizeof(igraph_vector_int_t));
        igraph_vector_int_t *bound = (igraph_vector_int_t*)malloc(sizeof(igraph_vector_int_t));
        if (!comp || !bound) {
            free(visited);
            free(boundary_marked);
            free(forbidden_out_expanded);
            free(queue);
            igraph_vector_int_destroy(&neighbors);
            igraph_vector_int_destroy(&parents);
            if (comp) free(comp);
            if (bound) free(bound);
            return IGRAPH_ENOMEM;
        }
        
        igraph_vector_int_init(comp, 0);
        igraph_vector_int_init(bound, 0);

        // BFS 初始化
        igraph_integer_t head = 0, tail = 0;
        queue[tail++] = v;
        visited[v] = 1;

        while (head < tail) {
            igraph_integer_t cur = queue[head++];

            igraph_vector_int_push_back(comp, cur);

            IGRAPH_CHECK(igraph_neighbors(graph, &neighbors, cur, IGRAPH_IN));
            igraph_integer_t nei_count = igraph_vector_int_size(&neighbors);

            for (igraph_integer_t i = 0; i < nei_count; i++) {
                igraph_integer_t w = VECTOR(neighbors)[i];
                if (visited[w] == 2) {
                    // 禁忌节点，加入边界（避免重复）
                    if (!boundary_marked[w]) {
                        igraph_vector_int_push_back(bound, w);
                        boundary_marked[w] = 1;
                    }
                    continue; // 不加入队列
                }

                if (!visited[w]) {
                    visited[w] = 1;
                    queue[tail++] = w;
                }
            }

            IGRAPH_CHECK(igraph_neighbors(graph, &neighbors, cur, IGRAPH_OUT));
            igraph_integer_t out_count = igraph_vector_int_size(&neighbors);

            for (igraph_integer_t i = 0; i < out_count; i++) {
                igraph_integer_t w = VECTOR(neighbors)[i];
                if (visited[w] == 2) {
                    // 出边遇到禁忌节点：w加入边界；w的禁忌父节点加入边界，非禁忌父节点并入当前分量
                    if (!boundary_marked[w]) {
                        igraph_vector_int_push_back(bound, w);
                        boundary_marked[w] = 1;
                    }

                    if (!forbidden_out_expanded[w]) {
                        forbidden_out_expanded[w] = 1;

                        IGRAPH_CHECK(igraph_neighbors(graph, &parents, w, IGRAPH_IN));
                        igraph_integer_t pcount = igraph_vector_int_size(&parents);
                        for (igraph_integer_t pi = 0; pi < pcount; ++pi) {
                            igraph_integer_t p = VECTOR(parents)[pi];
                            if (visited[p] == 2) {
                                if (!boundary_marked[p]) {
                                    igraph_vector_int_push_back(bound, p);
                                    boundary_marked[p] = 1;
                                }
                            } else if (!visited[p]) {
                                visited[p] = 1;
                                queue[tail++] = p;
                            }
                        }
                    }


                    continue; // 不加入队列
                }

                if (!visited[w]) {
                    visited[w] = 1;
                    queue[tail++] = w;
                }
            }

        }

        igraph_vector_ptr_push_back(components, comp);
        igraph_vector_ptr_push_back(boundaries, bound);

        // 清空边界辅助标记，为下一个连通分量做准备
        igraph_integer_t bsize = igraph_vector_int_size(bound);
        for (igraph_integer_t i = 0; i < bsize; i++) {
            boundary_marked[VECTOR(*bound)[i]] = 0;
        }

        
    }

    free(visited);
    free(boundary_marked);
    free(forbidden_out_expanded);
    free(queue);
    igraph_vector_int_destroy(&neighbors);
    igraph_vector_int_destroy(&parents);

    return IGRAPH_SUCCESS;
}









/*
 * 下面是有向图凸包找寻算法
 * 下面是一些辅助数据结构和函数，用于递归图分解过程中的工作空间管理和节点对标记。
 * 这些结构和函数主要用于存储和操作在分解过程中需要频繁访问的节点状态、掩码和队列等信息，
 * 以提高算法的效率和可读性。
 */
static size_t dag_pair_bitmap_index(igraph_integer_t a, igraph_integer_t b, igraph_integer_t n) {
    igraph_integer_t lo = a < b ? a : b;
    igraph_integer_t hi = a < b ? b : a;
    size_t nsz = (size_t) n;
    size_t lo_sz = (size_t) lo;
    return (lo_sz * (2 * nsz - lo_sz - 1)) / 2 + ((size_t) hi - lo_sz - 1);
}

static igraph_error_t dag_pair_bitmap_mark_seen(
    unsigned char *pair_seen,
    igraph_integer_t n,
    igraph_integer_t a,
    igraph_integer_t b,
    igraph_bool_t *inserted
) {
    if (a == b) {
        if (inserted) {
            *inserted = 0;
        }
        return IGRAPH_SUCCESS;
    }

    if (!pair_seen || n < 2) {
        if (inserted) {
            *inserted = 1;
        }
        return IGRAPH_SUCCESS;
    }

    size_t idx = dag_pair_bitmap_index(a, b, n);
    size_t byte_idx = idx >> 3;
    unsigned char bit = (unsigned char) (1u << (idx & 7u));

    if (pair_seen[byte_idx] & bit) {
        if (inserted) {
            *inserted = 0;
        }
        return IGRAPH_SUCCESS;
    }

    pair_seen[byte_idx] |= bit;
    if (inserted) {
        *inserted = 1;
    }
    return IGRAPH_SUCCESS;
}

typedef struct dag_workspace_t {
    igraph_integer_t n;
    char *seed_mask;
    char *anc_uv_mask;
    char *mb_u_mask;
    char *reach_v_mask;
    char *an_re_mask;
    char *visited_down;
    char *visited_up;
    igraph_integer_t *stack;
    igraph_integer_t *queue_node;
    int *queue_dir;
    igraph_vector_int_t preds;
    igraph_vector_int_t neigh;
    igraph_vector_int_t child_preds;
    igraph_bool_t vectors_initialized;
} dag_workspace_t;


static igraph_error_t dag_workspace_init(dag_workspace_t *ws, igraph_integer_t n) {
    igraph_error_t ret = IGRAPH_SUCCESS;
    igraph_bool_t preds_initialized = 0;
    igraph_bool_t neigh_initialized = 0;
    igraph_bool_t child_preds_initialized = 0;

    memset(ws, 0, sizeof(*ws));
    ws->n = n;

    ws->seed_mask = (char *) calloc((size_t) n, sizeof(char));
    ws->anc_uv_mask = (char *) calloc((size_t) n, sizeof(char));
    ws->mb_u_mask = (char *) calloc((size_t) n, sizeof(char));
    ws->reach_v_mask = (char *) calloc((size_t) n, sizeof(char));
    ws->an_re_mask = (char *) calloc((size_t) n, sizeof(char));
    ws->visited_down = (char *) calloc((size_t) n, sizeof(char));
    ws->visited_up = (char *) calloc((size_t) n, sizeof(char));
    ws->stack = (igraph_integer_t *) malloc((size_t) n * sizeof(igraph_integer_t));
    ws->queue_node = (igraph_integer_t *) malloc((size_t) (2 * n) * sizeof(igraph_integer_t));
    ws->queue_dir = (int *) malloc((size_t) (2 * n) * sizeof(int));

    if (!ws->seed_mask || !ws->anc_uv_mask || !ws->mb_u_mask || !ws->reach_v_mask ||
        !ws->an_re_mask || !ws->visited_down || !ws->visited_up || !ws->stack ||
        !ws->queue_node || !ws->queue_dir) {
        ret = IGRAPH_ENOMEM;
        goto fail;
    }

    ret = igraph_vector_int_init(&ws->preds, 0);
    if (ret != IGRAPH_SUCCESS) {
        goto fail;
    }
    preds_initialized = 1;

    ret = igraph_vector_int_init(&ws->neigh, 0);
    if (ret != IGRAPH_SUCCESS) {
        goto fail;
    }
    neigh_initialized = 1;

    ret = igraph_vector_int_init(&ws->child_preds, 0);
    if (ret != IGRAPH_SUCCESS) {
        goto fail;
    }
    child_preds_initialized = 1;

    ws->vectors_initialized = 1;
    return IGRAPH_SUCCESS;

fail:
    if (child_preds_initialized) {
        igraph_vector_int_destroy(&ws->child_preds);
    }
    if (neigh_initialized) {
        igraph_vector_int_destroy(&ws->neigh);
    }
    if (preds_initialized) {
        igraph_vector_int_destroy(&ws->preds);
    }

    free(ws->seed_mask);
    free(ws->anc_uv_mask);
    free(ws->mb_u_mask);
    free(ws->reach_v_mask);
    free(ws->an_re_mask);
    free(ws->visited_down);
    free(ws->visited_up);
    free(ws->stack);
    free(ws->queue_node);
    free(ws->queue_dir);
    memset(ws, 0, sizeof(*ws));
    return ret;
}

static void dag_workspace_destroy(dag_workspace_t *ws) {
    if (ws->vectors_initialized) {
        igraph_vector_int_destroy(&ws->preds);
        igraph_vector_int_destroy(&ws->neigh);
        igraph_vector_int_destroy(&ws->child_preds);
    }

    free(ws->seed_mask);
    free(ws->anc_uv_mask);
    free(ws->mb_u_mask);
    free(ws->reach_v_mask);
    free(ws->an_re_mask);
    free(ws->visited_down);
    free(ws->visited_up);
    free(ws->stack);
    free(ws->queue_node);
    free(ws->queue_dir);
    memset(ws, 0, sizeof(*ws));
}

static igraph_error_t dag_ancestors_from_mask_ws(
    const igraph_t *graph,
    const char *seed_mask,
    const char *allowed_mask,
    char *out_mask,
    igraph_integer_t *stack,
    igraph_vector_int_t *preds
) {
    igraph_integer_t n = igraph_vcount(graph);
    memset(out_mask, 0, (size_t) n);

    igraph_integer_t top = 0;
    for (igraph_integer_t i = 0; i < n; ++i) {
        if (seed_mask[i] && (!allowed_mask || allowed_mask[i])) {
            out_mask[i] = 1;
            stack[top++] = i;
        }
    }

    igraph_error_t ret = IGRAPH_SUCCESS;
    while (top > 0) {
        igraph_integer_t cur = stack[--top];
        ret = igraph_neighbors(graph, preds, cur, IGRAPH_IN);
        if (ret != IGRAPH_SUCCESS) {
            break;
        }
        igraph_integer_t pn = igraph_vector_int_size(preds);
        for (igraph_integer_t i = 0; i < pn; ++i) {
            igraph_integer_t pa = VECTOR(*preds)[i];
            if ((!allowed_mask || allowed_mask[pa]) && !out_mask[pa]) {
                out_mask[pa] = 1;
                stack[top++] = pa;
            }
        }
    }

    return ret;
}

static igraph_error_t dag_ancestors_from_mask(
    const igraph_t *graph,
    const char *seed_mask,
    const char *allowed_mask,
    char *out_mask
) {
    igraph_integer_t n = igraph_vcount(graph);
    igraph_integer_t *stack = (igraph_integer_t *) malloc((size_t) n * sizeof(igraph_integer_t));
    if (!stack) {
        return IGRAPH_ENOMEM;
    }

    igraph_vector_int_t preds;
    igraph_vector_int_init(&preds, 0);

    igraph_error_t ret = dag_ancestors_from_mask_ws(graph, seed_mask, allowed_mask, out_mask, stack, &preds);

    igraph_vector_int_destroy(&preds);
    free(stack);
    return ret;
}

static igraph_error_t dag_ancestors_from_vector(
    const igraph_t *graph,
    const igraph_vector_int_t *seeds,
    const char *allowed_mask,
    char *out_mask
) {
    igraph_integer_t n = igraph_vcount(graph);
    char *seed_mask = (char *) calloc((size_t) n, sizeof(char));
    if (!seed_mask) {
        return IGRAPH_ENOMEM;
    }

    igraph_integer_t k = igraph_vector_int_size(seeds);
    for (igraph_integer_t i = 0; i < k; ++i) {
        igraph_integer_t v = VECTOR(*seeds)[i];
        if (v >= 0 && v < n) {
            seed_mask[v] = 1;
        }
    }

    igraph_error_t ret = dag_ancestors_from_mask(graph, seed_mask, allowed_mask, out_mask);
    free(seed_mask);
    return ret;
}

igraph_error_t dag_get_ancestors(
    const igraph_t *graph,
    const igraph_vector_int_t *r_nodes,
    igraph_vector_int_t *ancestors_out
) {
    if (!graph || !r_nodes || !ancestors_out) {
        return IGRAPH_EINVAL;
    }

    igraph_integer_t n = igraph_vcount(graph);
    char *anc_mask = (char *) calloc((size_t) n, sizeof(char));
    if (!anc_mask) {
        return IGRAPH_ENOMEM;
    }

    igraph_error_t ret = dag_ancestors_from_vector(graph, r_nodes, NULL, anc_mask);
    if (ret == IGRAPH_SUCCESS) {
        igraph_vector_int_clear(ancestors_out);
        for (igraph_integer_t i = 0; i < n; ++i) {
            if (anc_mask[i]) {
                igraph_vector_int_push_back(ancestors_out, i);
            }
        }
    }

    free(anc_mask);
    return ret;
}

static int dag_pass_rule(int direction_in, igraph_integer_t node, int direction_out, const char *an_re_mask, const char *re_mask) {
    if (!an_re_mask[node]) {
        return !(direction_in == 0 && direction_out == 1);
    }
    if (re_mask[node]) {
        return (direction_in == 0 && direction_out == 1);
    }
    return 1;
}

static igraph_error_t dag_rech(
    const igraph_t *graph,
    igraph_integer_t source,
    const char *re_mask,
    const char *allowed_mask,
    char *reach_mask,
    dag_workspace_t *ws
) {
    igraph_integer_t n = igraph_vcount(graph);

    igraph_error_t ret = dag_ancestors_from_mask_ws(
        graph,
        re_mask,
        allowed_mask,
        ws->an_re_mask,
        ws->stack,
        &ws->preds
    );
    if (ret != IGRAPH_SUCCESS) {
        return ret;
    }

    igraph_integer_t head = 0, tail = 0;
    memset(ws->visited_down, 0, (size_t) n);
    memset(ws->visited_up, 0, (size_t) n);
    memset(reach_mask, 0, (size_t) n);

    if (allowed_mask && !allowed_mask[source]) {
        ret = IGRAPH_SUCCESS;
        return ret;
    }

    ret = igraph_neighbors(graph, &ws->neigh, source, IGRAPH_OUT);
    if (ret != IGRAPH_SUCCESS) {
        return ret;
    }
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(&ws->neigh); ++i) {
        igraph_integer_t ch = VECTOR(ws->neigh)[i];
        if ((!allowed_mask || allowed_mask[ch]) && !ws->visited_down[ch]) {
            ws->visited_down[ch] = 1;
            ws->queue_node[tail] = ch;
            ws->queue_dir[tail] = 0;
            tail++;
        }
    }

    ret = igraph_neighbors(graph, &ws->neigh, source, IGRAPH_IN);
    if (ret != IGRAPH_SUCCESS) {
        return ret;
    }
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(&ws->neigh); ++i) {
        igraph_integer_t pa = VECTOR(ws->neigh)[i];
        if ((!allowed_mask || allowed_mask[pa]) && !ws->visited_up[pa]) {
            ws->visited_up[pa] = 1;
            ws->queue_node[tail] = pa;
            ws->queue_dir[tail] = 1;
            tail++;
        }
    }

    while (head < tail) {
        igraph_integer_t node = ws->queue_node[head];
        int dir_in = ws->queue_dir[head];
        head++;

        ret = igraph_neighbors(graph, &ws->neigh, node, IGRAPH_OUT);
        if (ret != IGRAPH_SUCCESS) {
            return ret;
        }
        for (igraph_integer_t i = 0; i < igraph_vector_int_size(&ws->neigh); ++i) {
            igraph_integer_t nb = VECTOR(ws->neigh)[i];
            if ((!allowed_mask || allowed_mask[nb]) && !ws->visited_down[nb] && dag_pass_rule(dir_in, node, 0, ws->an_re_mask, re_mask)) {
                ws->visited_down[nb] = 1;
                ws->queue_node[tail] = nb;
                ws->queue_dir[tail] = 0;
                tail++;
            }
        }

        ret = igraph_neighbors(graph, &ws->neigh, node, IGRAPH_IN);
        if (ret != IGRAPH_SUCCESS) {
            return ret;
        }
        for (igraph_integer_t i = 0; i < igraph_vector_int_size(&ws->neigh); ++i) {
            igraph_integer_t nb = VECTOR(ws->neigh)[i];
            if ((!allowed_mask || allowed_mask[nb]) && !ws->visited_up[nb] && dag_pass_rule(dir_in, node, 1, ws->an_re_mask, re_mask)) {
                ws->visited_up[nb] = 1;
                ws->queue_node[tail] = nb;
                ws->queue_dir[tail] = 1;
                tail++;
            }
        }
    }

    for (igraph_integer_t i = 0; i < n; ++i) {
        reach_mask[i] = (char) (ws->visited_down[i] || ws->visited_up[i]);
    }

    return ret;
}

static igraph_error_t dag_fcms(
    const igraph_t *graph,
    igraph_integer_t u,
    igraph_integer_t v,
    const char *allowed_mask,
    igraph_vector_int_t *out_nodes,
    dag_workspace_t *ws
) {
    igraph_integer_t n = igraph_vcount(graph);

    memset(ws->seed_mask, 0, (size_t) n);
    memset(ws->mb_u_mask, 0, (size_t) n);
    ws->seed_mask[u] = 1;
    ws->seed_mask[v] = 1;

    igraph_error_t ret = dag_ancestors_from_mask_ws(
        graph,
        ws->seed_mask,
        allowed_mask,
        ws->anc_uv_mask,
        ws->stack,
        &ws->preds
    );
    if (ret != IGRAPH_SUCCESS) {
        return ret;
    }

    ret = igraph_neighbors(graph, &ws->neigh, u, IGRAPH_IN);
    if (ret != IGRAPH_SUCCESS) {
        return ret;
    }
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(&ws->neigh); ++i) {
        igraph_integer_t p = VECTOR(ws->neigh)[i];
        if ((!allowed_mask || allowed_mask[p]) && ws->anc_uv_mask[p]) {
            ws->mb_u_mask[p] = 1;
        }
    }

    ret = igraph_neighbors(graph, &ws->neigh, u, IGRAPH_OUT);
    if (ret != IGRAPH_SUCCESS) {
        return ret;
    }
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(&ws->neigh); ++i) {
        igraph_integer_t c = VECTOR(ws->neigh)[i];
        if ((allowed_mask && !allowed_mask[c]) || !ws->anc_uv_mask[c]) {
            continue;
        }
        ws->mb_u_mask[c] = 1;

        ret = igraph_neighbors(graph, &ws->child_preds, c, IGRAPH_IN);
        if (ret != IGRAPH_SUCCESS) {
            return ret;
        }
        for (igraph_integer_t j = 0; j < igraph_vector_int_size(&ws->child_preds); ++j) {
            igraph_integer_t cp = VECTOR(ws->child_preds)[j];
            if ((!allowed_mask || allowed_mask[cp]) && ws->anc_uv_mask[cp]) {
                ws->mb_u_mask[cp] = 1;
            }
        }
    }
    ws->mb_u_mask[u] = 0;

    ret = dag_rech(graph, v, ws->mb_u_mask, allowed_mask, ws->reach_v_mask, ws);
    if (ret != IGRAPH_SUCCESS) {
        return ret;
    }

    igraph_vector_int_clear(out_nodes);
    for (igraph_integer_t i = 0; i < n; ++i) {
        if ((!allowed_mask || allowed_mask[i]) && ws->mb_u_mask[i] && ws->reach_v_mask[i]) {
            igraph_vector_int_push_back(out_nodes, i);
        }
    }

    return ret;
}

igraph_error_t cmdsa(
    const igraph_t *graph,
    const igraph_vector_int_t *r_nodes,
    igraph_vector_int_t *H_out
) {
    igraph_bool_t directed = igraph_is_directed(graph);
    if (!directed) {
        return IGRAPH_EINVAL;
    }

    igraph_error_t ret = IGRAPH_SUCCESS;
    igraph_integer_t rsize = igraph_vector_int_size(r_nodes);


    igraph_vector_int_t anc_nodes;
    igraph_vector_int_t tmp_anc;
    igraph_vector_int_t topo_order;
    igraph_vector_int_t h_nodes;
    igraph_vector_int_t s_a;
    igraph_vector_int_t s_b;
    igraph_vector_int_t map;
    igraph_vector_int_t invmap;
    igraph_vector_int_init(&anc_nodes, 0);
    igraph_vector_int_init(&tmp_anc, 0);
    igraph_vector_int_init(&topo_order, 0);
    igraph_vector_int_init(&h_nodes, 0);
    igraph_vector_int_init(&s_a, 0);
    igraph_vector_int_init(&s_b, 0);
    igraph_vector_int_init(&map, 0);
    igraph_vector_int_init(&invmap, 0);

    unsigned char *in_h = NULL;
    unsigned char *pair_seen = NULL;

    igraph_vector_ptr_t components;
    igraph_vector_ptr_t boundaries;
    igraph_vector_ptr_init(&components, 0);
    igraph_vector_ptr_init(&boundaries, 0);

    igraph_t an_graph;
    igraph_bool_t an_graph_initialized = 0;
    dag_workspace_t ws;
    igraph_bool_t ws_initialized = 0;
    igraph_integer_t *topo_rank = NULL;
    igraph_integer_t an = 0;

    // 直接收集所有根节点的祖先到向量中
    for (igraph_integer_t i = 0; i < rsize; ++i) {
        igraph_integer_t rv = VECTOR(*r_nodes)[i];
        ret = igraph_subcomponent(graph, &tmp_anc, rv, IGRAPH_IN);
        if (ret != IGRAPH_SUCCESS) {
            goto cleanup;
        }
        for (igraph_integer_t j = 0; j < igraph_vector_int_size(&tmp_anc); ++j) {
            ret = igraph_vector_int_push_back(&anc_nodes, VECTOR(tmp_anc)[j]);
            if (ret != IGRAPH_SUCCESS) {
                goto cleanup;
            }
        }
    }

    // 排序后去重
    igraph_vector_int_sort(&anc_nodes);
    vector_int_unique(&anc_nodes);

    igraph_vs_t vs;
    igraph_vs_vector(&vs, &anc_nodes);
    ret = igraph_induced_subgraph_map(graph, &an_graph, vs, IGRAPH_SUBGRAPH_AUTO, &map, &invmap);
    igraph_vs_destroy(&vs);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }
    an_graph_initialized = 1;

    an = igraph_vcount(&an_graph);

    in_h = (unsigned char *) calloc((size_t) an, sizeof(unsigned char));
    if (!in_h) {
        ret = IGRAPH_ENOMEM;
        goto cleanup;
    }

    {
        size_t an_sz = (size_t) an;
        size_t pair_bits = an_sz > 1 ? (an_sz * (an_sz - 1)) / 2 : 0;
        size_t pair_bytes = (pair_bits + 7) / 8;
        if (pair_bytes > 0) {
            pair_seen = (unsigned char *) calloc(pair_bytes, sizeof(unsigned char));
            if (!pair_seen) {
                ret = IGRAPH_ENOMEM;
                goto cleanup;
            }
        }
    }

    for (igraph_integer_t i = 0; i < rsize; ++i) {
        igraph_integer_t rv = VECTOR(*r_nodes)[i];
        igraph_integer_t lv = VECTOR(map)[rv] - 1;
        if (lv < 0 || lv >= an) {
            continue;
        }
        if (!in_h[lv]) {
            in_h[lv] = 1;
            ret = igraph_vector_int_push_back(&h_nodes, lv);
            if (ret != IGRAPH_SUCCESS) {
                goto cleanup;
            }
        }
    }

    ret = igraph_topological_sorting(&an_graph, &topo_order, IGRAPH_OUT);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    topo_rank = (igraph_integer_t *) malloc((size_t) an * sizeof(igraph_integer_t));
    if (!topo_rank) {
        ret = IGRAPH_ENOMEM;
        goto cleanup;
    }
    for (igraph_integer_t i = 0; i < an; ++i) {
        topo_rank[VECTOR(topo_order)[i]] = i;
    }

    if (dag_workspace_init(&ws, an) != IGRAPH_SUCCESS) {
        ret = IGRAPH_ENOMEM;
        goto cleanup;
    }
    ws_initialized = 1;

    igraph_bool_t changed = 1;
    while (changed) {
        changed = 0;

        for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&components); i++) {
            igraph_vector_int_t *comp = (igraph_vector_int_t *) VECTOR(components)[i];
            igraph_vector_int_destroy(comp);
            free(comp);
        }
        for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&boundaries); i++) {
            igraph_vector_int_t *bound = (igraph_vector_int_t *) VECTOR(boundaries)[i];
            igraph_vector_int_destroy(bound);
            free(bound);
        }
        igraph_vector_ptr_clear(&components);
        igraph_vector_ptr_clear(&boundaries);

        ret = dag_components_forbidden(&an_graph, &components, &boundaries, &h_nodes);
        if (ret != IGRAPH_SUCCESS) {
            break;
        }

        igraph_integer_t ngroups = igraph_vector_ptr_size(&boundaries);
        for (igraph_integer_t cidx = 0; cidx < ngroups && !changed; ++cidx) {
            igraph_vector_int_t *boundary = (igraph_vector_int_t *) VECTOR(boundaries)[cidx];
            if (igraph_vector_int_size(boundary) < 2) {
                continue;
            }

            igraph_vector_int_sort(boundary);
            igraph_integer_t bsize = igraph_vector_int_size(boundary);

            for (igraph_integer_t i = 0; i < bsize && !changed; ++i) {
                igraph_integer_t a = VECTOR(*boundary)[i];
                for (igraph_integer_t j = i + 1; j < bsize && !changed; ++j) {
                    igraph_integer_t b = VECTOR(*boundary)[j];
                    igraph_integer_t aa = a;
                    igraph_integer_t bb = b;

                    if (topo_rank[aa] > topo_rank[bb]) {
                        igraph_integer_t t = aa;
                        aa = bb;
                        bb = t;
                    }
                    if (topo_rank[aa] >= topo_rank[bb]) {
                        continue;
                    }

                    igraph_bool_t pair_is_new = 0;
                    ret = dag_pair_bitmap_mark_seen(pair_seen, an, aa, bb, &pair_is_new);
                    if (ret != IGRAPH_SUCCESS) {
                        break;
                    }
                    if (!pair_is_new) {
                        continue;
                    }

                    igraph_bool_t ab = 0;
                    ret = igraph_are_adjacent(&an_graph, aa, bb, &ab);
                    if (ret != IGRAPH_SUCCESS) {
                        break;
                    }
                    if (ab) {
                        continue;
                    }

                    igraph_vector_int_clear(&s_a);
                    igraph_vector_int_clear(&s_b);

                    ret = dag_fcms(&an_graph, aa, bb, NULL, &s_a, &ws);
                    if (ret == IGRAPH_SUCCESS) {
                        ret = dag_fcms(&an_graph, bb, aa, NULL, &s_b, &ws);
                    }
                    if (ret != IGRAPH_SUCCESS) {
                        break;
                    }

                    igraph_bool_t local_changed = 0;
                    for (igraph_integer_t si = 0; si < igraph_vector_int_size(&s_a); ++si) {
                        igraph_integer_t node = VECTOR(s_a)[si];
                        if (node >= 0 && node < an && !in_h[node]) {
                            in_h[node] = 1;
                            ret = igraph_vector_int_push_back(&h_nodes, node);
                            if (ret != IGRAPH_SUCCESS) {
                                break;
                            }
                            local_changed = 1;
                        }
                    }
                    if (ret != IGRAPH_SUCCESS) {
                        break;
                    }

                    for (igraph_integer_t si = 0; si < igraph_vector_int_size(&s_b); ++si) {
                        igraph_integer_t node = VECTOR(s_b)[si];
                        if (node >= 0 && node < an && !in_h[node]) {
                            in_h[node] = 1;
                            ret = igraph_vector_int_push_back(&h_nodes, node);
                            if (ret != IGRAPH_SUCCESS) {
                                break;
                            }
                            local_changed = 1;
                        }
                    }
                    if (ret != IGRAPH_SUCCESS) {
                        break;
                    }

                    if (local_changed) {
                        changed = 1;
                    }
                }
            }
        }

        if (ret != IGRAPH_SUCCESS) {
            break;
        }
    }

    if (ret == IGRAPH_SUCCESS) {
        igraph_vector_int_sort(&h_nodes);
        igraph_vector_int_clear(H_out);
        for (igraph_integer_t i = 0; i < igraph_vector_int_size(&h_nodes); ++i) {
            igraph_integer_t lv = VECTOR(h_nodes)[i];
            ret = igraph_vector_int_push_back(H_out, VECTOR(invmap)[lv]);
            if (ret != IGRAPH_SUCCESS) {
                break;
            }
        }
    }

cleanup:
    if (ws_initialized) {
        dag_workspace_destroy(&ws);
    }
    if (an_graph_initialized) {
        igraph_destroy(&an_graph);
    }

    free(in_h);
    free(pair_seen);
    free(topo_rank);

    for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&components); i++) {
        igraph_vector_int_t *comp = (igraph_vector_int_t *) VECTOR(components)[i];
        igraph_vector_int_destroy(comp);
        free(comp);
    }
    for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&boundaries); i++) {
        igraph_vector_int_t *bound = (igraph_vector_int_t *) VECTOR(boundaries)[i];
        igraph_vector_int_destroy(bound);
        free(bound);
    }
    igraph_vector_ptr_destroy(&components);
    igraph_vector_ptr_destroy(&boundaries);

    igraph_vector_int_destroy(&anc_nodes);
    igraph_vector_int_destroy(&tmp_anc);
    igraph_vector_int_destroy(&topo_order);
    igraph_vector_int_destroy(&h_nodes);
    igraph_vector_int_destroy(&s_a);
    igraph_vector_int_destroy(&s_b);
    igraph_vector_int_destroy(&map);
    igraph_vector_int_destroy(&invmap);

    return ret;
}
