#define USE_RINTERNALS /* TODO RESTORE: Eventually faster, but may make bugs harder to spot*/

#include <Rmath.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> /* for NULL */
#include <R_ext/Rdynload.h>

extern SEXP RogueNaRok(SEXP, SEXP, SEXP, SEXP, SEXP,
                       SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef callMethods[] = {
  {"RogueNaRok", (DL_FUNC) &RogueNaRok, 10},
  {NULL, NULL, 0}
};

void R_init_Rogue(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
