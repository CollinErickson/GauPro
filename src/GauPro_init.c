#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP GauPro_cholC(SEXP);
extern SEXP GauPro_corr_gauss_matrix_sym_armaC(SEXP, SEXP);
extern SEXP GauPro_corr_gauss_matrix_symC(SEXP, SEXP);
extern SEXP GauPro_corr_gauss_matrixC(SEXP, SEXP, SEXP);
extern SEXP GauPro_corr_gauss_matrixvecC(SEXP, SEXP, SEXP);
extern SEXP GauPro_deviance_fngr_joint(SEXP, SEXP, SEXP);
extern SEXP GauPro_deviance_fngr_nug(SEXP, SEXP, SEXP);
extern SEXP GauPro_deviance_fngr_theta(SEXP, SEXP, SEXP);
extern SEXP GauPro_deviance_grad_joint(SEXP, SEXP, SEXP, SEXP);
extern SEXP GauPro_deviance_grad_nug(SEXP, SEXP, SEXP, SEXP);
extern SEXP GauPro_deviance_grad_theta(SEXP, SEXP, SEXP, SEXP);
extern SEXP GauPro_deviance_part(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GauPro_devianceC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GauPro_Gaussian_deviance_part(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GauPro_Gaussian_devianceC(SEXP, SEXP, SEXP, SEXP);
extern SEXP GauPro_Gaussian_hessianCC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GauPro_pred_cov(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GauPro_pred_meanC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GauPro_pred_var(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GauPro_rcpp_hello_world();
extern SEXP GauPro_solveC(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"GauPro_cholC",                       (DL_FUNC) &GauPro_cholC,                       1},
  {"GauPro_corr_gauss_matrix_sym_armaC", (DL_FUNC) &GauPro_corr_gauss_matrix_sym_armaC, 2},
  {"GauPro_corr_gauss_matrix_symC",      (DL_FUNC) &GauPro_corr_gauss_matrix_symC,      2},
  {"GauPro_corr_gauss_matrixC",          (DL_FUNC) &GauPro_corr_gauss_matrixC,          3},
  {"GauPro_corr_gauss_matrixvecC",       (DL_FUNC) &GauPro_corr_gauss_matrixvecC,       3},
  {"GauPro_deviance_fngr_joint",         (DL_FUNC) &GauPro_deviance_fngr_joint,         3},
  {"GauPro_deviance_fngr_nug",           (DL_FUNC) &GauPro_deviance_fngr_nug,           3},
  {"GauPro_deviance_fngr_theta",         (DL_FUNC) &GauPro_deviance_fngr_theta,         3},
  {"GauPro_deviance_grad_joint",         (DL_FUNC) &GauPro_deviance_grad_joint,         4},
  {"GauPro_deviance_grad_nug",           (DL_FUNC) &GauPro_deviance_grad_nug,           4},
  {"GauPro_deviance_grad_theta",         (DL_FUNC) &GauPro_deviance_grad_theta,         4},
  {"GauPro_deviance_part",               (DL_FUNC) &GauPro_deviance_part,               5},
  {"GauPro_devianceC",                   (DL_FUNC) &GauPro_devianceC,                   5},
  {"GauPro_Gaussian_deviance_part",      (DL_FUNC) &GauPro_Gaussian_deviance_part,      5},
  {"GauPro_Gaussian_devianceC",          (DL_FUNC) &GauPro_Gaussian_devianceC,          4},
  {"GauPro_Gaussian_hessianCC",          (DL_FUNC) &GauPro_Gaussian_hessianCC,          6},
  {"GauPro_pred_cov",                    (DL_FUNC) &GauPro_pred_cov,                    6},
  {"GauPro_pred_meanC",                  (DL_FUNC) &GauPro_pred_meanC,                  5},
  {"GauPro_pred_var",                    (DL_FUNC) &GauPro_pred_var,                    6},
  {"GauPro_rcpp_hello_world",            (DL_FUNC) &GauPro_rcpp_hello_world,            0},
  {"GauPro_solveC",                      (DL_FUNC) &GauPro_solveC,                      2},
  {NULL, NULL, 0}
};

void R_init_GauPro(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  //R_registerRoutines(dll, NULL, NULL, NULL, NULL);
  //R_useDynamicSymbols(dll, TRUE);
}
