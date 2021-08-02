#define USE_RINTERNALS

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> /* for NULL */
#include <R_ext/Rdynload.h>

extern SEXP RogueNaRok(SEXP, SEXP, SEXP, SEXP, SEXP,
                       SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP COPHENETIC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP COPHENETIC_LOG(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef callMethods[] = {
  {"COPHENETIC", (DL_FUNC) &COPHENETIC, 5},
  {"COPHENETIC_LOG", (DL_FUNC) &COPHENETIC_LOG, 5},
  {"RogueNaRok", (DL_FUNC) &RogueNaRok, 10},
  {NULL, NULL, 0}
};

void R_init_Rogue(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
