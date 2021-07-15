
SEXP COPHENETIC(SEXP n_tip, SEXP n_node, SEXP parent, SEXP child, SEXP n_edge) {
  const int all_nodes = INTEGER(n_tip)[0] + INTEGER(n_node)[0];
  SEXP RESULT = PROTECT(allocVector(INTSXP, all_nodes * all_nodes));
  int *result = INTEGER(RESULT);

  cophenetic_phylo(INTEGER(n_tip), INTEGER(n_node), INTEGER(parent),
                   INTEGER(child), INTEGER(n_edge), &all_nodes, result);
  UNPROTECT(1);
  return(RESULT);
}
