#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _JACA_jacaCpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _JACA_transformYCpp(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_JACA_jacaCpp",       (DL_FUNC) &_JACA_jacaCpp,       8},
  {"_JACA_transformYCpp", (DL_FUNC) &_JACA_transformYCpp, 1},
  {NULL, NULL, 0}
};

void R_init_JACA(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
