#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP PhylogeneticEM_log_likelihood_BM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PhylogeneticEM_log_likelihood_OU(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PhylogeneticEM_upward_downward_BM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PhylogeneticEM_upward_downward_OU(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"PhylogeneticEM_log_likelihood_BM",  (DL_FUNC) &PhylogeneticEM_log_likelihood_BM,  6},
  {"PhylogeneticEM_log_likelihood_OU",  (DL_FUNC) &PhylogeneticEM_log_likelihood_OU,  7},
  {"PhylogeneticEM_upward_downward_BM", (DL_FUNC) &PhylogeneticEM_upward_downward_BM, 6},
  {"PhylogeneticEM_upward_downward_OU", (DL_FUNC) &PhylogeneticEM_upward_downward_OU, 7},
  {NULL, NULL, 0}
};

void R_init_PhylogeneticEM(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}