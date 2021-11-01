#define USE_RINTERNALS

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> /* for NULL */
#include <R_ext/Rdynload.h>

extern SEXP RogueNaRok(SEXP, SEXP, SEXP, SEXP, SEXP,
                       SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GRAPH_GEODESIC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LOG_GRAPH_GEODESIC(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef callMethods[] = {
  {"GRAPH_GEODESIC", (DL_FUNC) &GRAPH_GEODESIC, 5},
  {"LOG_GRAPH_GEODESIC", (DL_FUNC) &LOG_GRAPH_GEODESIC, 5},
  {"RogueNaRok", (DL_FUNC) &RogueNaRok, 10},
  {NULL, NULL, 0}
};

void R_init_Rogue(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
