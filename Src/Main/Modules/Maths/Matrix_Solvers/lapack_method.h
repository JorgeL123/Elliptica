#include "core_lib.h"
#include "error_handling_lib.h"
#include "maths_linear_algebra_lib.h"
#include "maths_matrix_solvers_lib.h"

/* LAPACK (Fortran interface, LP64), provided by the linked BLAS library. */
void dgetrf_(const int *m,const int *n,double *a,const int *lda,
             int *ipiv,int *info);
void dgetrs_(const char *trans,const int *n,const int *nrhs,
             const double *a,const int *lda,const int *ipiv,
             double *b,const int *ldb,int *info,size_t trans_len);
void dgecon_(const char *norm,const int *n,const double *a,const int *lda,
             const double *anorm,double *rcond,double *work,int *iwork,
             int *info,size_t norm_len);

int direct_solver_series_lapack_dense(void *vp);
Lapack_Dense_T *init_lapack_dense(void);
void free_lapack_dense(Lapack_Dense_T *lpk);
