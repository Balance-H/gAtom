import os
import sys
import time

import igraph
import netdecom as nd

os.environ['PATH'] = r"D:\vcpkg-master\installed\x64-windows\bin;" + os.environ['PATH']
sys.path.append(r"D:\Jupyter\A_final_algorithms\netdecom_c\build\python_modules")
import decom_h


def normalize_sets(items):
    return {frozenset(item) for item in items}


def check_tree_and_separators(atoms, separators, tree_edges):
    atom_sets = [set(atom) for atom in atoms]
    separator_sets = [set(sep) for sep in separators]

    tree_graph = igraph.Graph(n=len(atoms), edges=tree_edges, directed=False)
    tree_is_connected = tree_graph.is_connected()
    tree_is_tree = tree_is_connected and tree_graph.ecount() == max(0, len(atoms) - 1)

    edge_intersections = []
    for u, v in tree_edges:
        edge_intersections.append(atom_sets[u] & atom_sets[v])

    intersection_set = normalize_sets(edge_intersections)
    separator_set = normalize_sets(separator_sets)

    return {
        "tree_is_connected": tree_is_connected,
        "tree_is_tree": tree_is_tree,
        "edge_intersections": edge_intersections,
        "intersection_set": intersection_set,
        "separator_set": separator_set,
        "separators_match": intersection_set == separator_set,
    }


filenames = ["Animal-Network.txt", "bio-CE-GT.txt", "CA-HepTh.txt"]
repeat = 1
for f in filenames:
    G_ig = nd.get_example(f, class_type="ig")
    components = G_ig.components()
    k = len(components)
    times = []
    for i in range(k - 1):
        u = components[i][0]
        v = components[i + 1][0]
        G_ig.add_edge(u, v)

    for i in range(repeat):
        start = time.perf_counter()
        atoms, separators, tree = decom_h.decompose_atoms(G_ig)
        end = time.perf_counter()
        times.append(end - start)

    avg_time = sum(times) / repeat
    checks = check_tree_and_separators(atoms, separators, tree)

    print(
        f"CMSA平均耗时: {avg_time:.6f} 秒",
        f"原子数量: {len(atoms)}",
        f"分离子数量: {len(separators)}",
        f"树宽: {max(len(atom) for atom in atoms) - 1}",
        f"tree连通: {checks['tree_is_connected']}",
        f"tree是树: {checks['tree_is_tree']}",
        f"边交集与separators一致: {checks['separators_match']}"
    )

    if not checks["separators_match"]:
        print("  边交集集合:", checks["intersection_set"])
        print("  separators集合:", checks["separator_set"])
