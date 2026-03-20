#include <R.h>
#include <Rinternals.h>
#include <assert.h>

#define GEOD_MAX 32768

double lg[GEOD_MAX];
__attribute__((constructor))
  void compute_logs(void) {
    for (int i = GEOD_MAX; --i; ) {
      lg[i] = log(i);
    }
  }

#define RET(i, j) ret[(i) + *all_nodes * (j)]
#define GET(i, j) ((i) > (j) ? RET((i), (j)) : RET((j), (i)))
#define SETBOTH(i, j, k) RET((i), (j)) = RET((j), (i)) = (k)
#define SETNODE(i, j, k) if ((i) > (j)) RET((i), (j)) = (k); else RET((j), (i)) = (k)

// Algorithm following ape::dist_nodes
void graph_geodesic_phylo(const int *n_tip, const int *n_node,
                          const int *parent, const int *child,
                          const int *n_edge, const int *all_nodes,
                          int *ret) {
  const int
    root_node = parent[0],
    root_child = child[0]
  ;
  for (int i = *all_nodes; i--; ) {
    RET(i, i) = 0;
  }

#ifdef MYDEBUG
  Rprintf("Start: SetN %i > %i = %i\n", root_node, root_child, 1);
#endif
  SETNODE(root_node, root_child, 1);
  for (int i = 1; i != *n_edge; ++i) {
    const int parent_i = parent[i], child_i = child[i];
#ifdef MYDEBUG
    Rprintf("\n edge_i = %i: SetN %i > %i = %i\n", i, parent_i, child_i, 1);
#endif
    SETNODE(parent_i, child_i, 1);

    for (int j = i - 1; j >= 0; --j) {
      const int child_j = child[j];
      if (child_j == parent_i) {
#ifdef MYDEBUG
        Rprintf("j = %i: Skip %i = %i\n", j, parent_i, child_j);
#endif
        continue;
      }
#ifdef MYDEBUG
      Rprintf("j = %i: SetB [%i]-%i = 1+ (%i - %i) =  %i\n", j, child_j, child_i,
              parent_i, child_j, GET(parent_i, child_j) + 1);
#endif
      SETBOTH(child_j, child_i, GET(parent_i, child_j) + 1);
    }
    if (parent_i != root_node) {
#ifdef MYDEBUG
      Rprintf("Finish edge %i: SetN %i-%i = 1+ %i-%i = %i\n", i, root_node,
              child_i, root_node, parent_i, RET(parent_i, root_node) + 1);
      if (parent_i <= root_node) {
        REprintf("Something went wrong! parent_i = %i <= root_node = %i",
                 parent_i, root_node);
      }
#endif
      assert(parent_i > root_node);
      if (child_i > root_node) {
        SETNODE(child_i, root_node, RET(parent_i, root_node) + 1);
      } else {
        SETNODE(root_node, child_i, RET(parent_i, root_node) + 1);
      }
    }
  }
}

SEXP GRAPH_GEODESIC(SEXP n_tip, SEXP n_node, SEXP parent, SEXP child,
                    SEXP n_edge) {
  const int all_nodes = INTEGER(n_tip)[0] + INTEGER(n_node)[0];
  SEXP RESULT = PROTECT(allocVector(INTSXP, all_nodes * all_nodes));
  int *result = INTEGER(RESULT);

  graph_geodesic_phylo(INTEGER(n_tip), INTEGER(n_node), INTEGER(parent),
                       INTEGER(child), INTEGER(n_edge), &all_nodes, result);
  UNPROTECT(1);
  return(RESULT);
}

SEXP LOG_GRAPH_GEODESIC(SEXP n_tip, SEXP n_node, SEXP parent, SEXP child,
                        SEXP n_edge) {
  const int
    n_tips = INTEGER(n_tip)[0],
    all_nodes = n_tips + INTEGER(n_node)[0]
  ;
  SEXP RESULT = PROTECT(allocVector(REALSXP, n_tips * n_tips));
  SEXP INTERIM = PROTECT(allocVector(INTSXP, all_nodes * all_nodes));
  int *interim = INTEGER(INTERIM);

  graph_geodesic_phylo(INTEGER(n_tip), INTEGER(n_node), INTEGER(parent),
                       INTEGER(child), INTEGER(n_edge), &all_nodes, interim);

  double *result = REAL(RESULT);
  for (int i = INTEGER(n_tip)[0]; i--; ) {
    for (int j = 0; j < i; ++j) {
      result[i + (n_tips * j)] =
        result[j + (n_tips * i)] =
        lg[interim[j + (all_nodes * i)]];
    }
    result[i + (n_tips * i)] = R_NegInf;
  }
  UNPROTECT(2);
  return(RESULT);
}

// Batch version: compute log-geodesics for multiple trees in one call.
// Returns lower-triangle entries only (nPairs x nTree matrix).
// Reuses a single interim buffer across all trees.
//
// parent_all, child_all: concatenated edge arrays for all trees
//   (each tree contributes n_edge consecutive entries)
// n_tree: number of trees
SEXP LOG_GRAPH_GEODESIC_MULTI(SEXP n_tip, SEXP n_node, SEXP parent_all,
                              SEXP child_all, SEXP n_edge, SEXP n_tree) {
  const int
    n_tips = INTEGER(n_tip)[0],
    n_nodes = INTEGER(n_node)[0],
    all_nodes = n_tips + n_nodes,
    n_edges = INTEGER(n_edge)[0],
    n_trees = INTEGER(n_tree)[0],
    n_pairs = n_tips * (n_tips - 1) / 2
  ;

  SEXP RESULT = PROTECT(allocVector(REALSXP, (R_xlen_t)n_pairs * n_trees));
  SEXP INTERIM = PROTECT(allocVector(INTSXP, all_nodes * all_nodes));
  double *result = REAL(RESULT);
  int *interim = INTEGER(INTERIM);
  const int *par_all = INTEGER(parent_all);
  const int *ch_all = INTEGER(child_all);

  for (int t = 0; t < n_trees; ++t) {
    const int *par = par_all + (R_xlen_t)t * n_edges;
    const int *ch = ch_all + (R_xlen_t)t * n_edges;

    graph_geodesic_phylo(&n_tips, &n_nodes, par, ch, &n_edges,
                         &all_nodes, interim);

    // Extract lower triangle (row > col in column-major order)
    // Read interim[i + all_nodes * j] (stride-1) instead of
    // interim[j + all_nodes * i] (stride all_nodes); matrix is symmetric.
    double *res_col = result + (R_xlen_t)t * n_pairs;
    int pair_idx = 0;
    for (int j = 0; j < n_tips - 1; ++j) {
      const int col_offset = all_nodes * j;
      for (int i = j + 1; i < n_tips; ++i) {
        res_col[pair_idx++] = lg[interim[i + col_offset]];
      }
    }
  }

  UNPROTECT(2);
  return(RESULT);
}
