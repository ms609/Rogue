#ifndef ROGUE_GEODESIC_H
#define ROGUE_GEODESIC_H

/* Shared declarations so the C++ translation unit (tip_instability.cpp) can
 * reuse the geodesic worker and log lookup table defined in graph_geodesic.c.
 * Wrapped in extern "C" when included from C++ so the symbols match the
 * (unmangled) C definitions. */

#ifdef __cplusplus
extern "C" {
#endif

/* Lookup table of natural logs, lg[k] == log(k); populated at load time by a
 * constructor in graph_geodesic.c. */
extern double lg[];

/* Fill `ret` (an all_nodes * all_nodes int matrix) with pairwise node graph
 * geodesics for one tree.  See graph_geodesic.c for the algorithm. */
void graph_geodesic_phylo(const int *n_tip, const int *n_node,
                          const int *parent, const int *child,
                          const int *n_edge, const int *all_nodes,
                          int *ret);

#ifdef __cplusplus
}
#endif

#endif  /* ROGUE_GEODESIC_H */
