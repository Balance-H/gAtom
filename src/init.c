#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP c_get_minimal_collapsible_impl(SEXP, SEXP, SEXP, SEXP);
extern SEXP c_decompose_atoms_impl(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_gAtom_get_minimal_collapsible_impl", (DL_FUNC) &c_get_minimal_collapsible_impl, 4},
    {"_gAtom_decompose_atoms_impl", (DL_FUNC) &c_decompose_atoms_impl, 3},
    {NULL, NULL, 0}
};

void R_init_gAtom(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
