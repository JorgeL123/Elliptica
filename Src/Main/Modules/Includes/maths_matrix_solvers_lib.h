#ifndef maths_matrix_solvers
#define maths_matrix_solvers
#include "elliptica_system_lib.h"

struct MATRIX_T;

/* umfpack direct solver */
typedef struct UMFPACK_T
{
  char description[9000];
  struct MATRIX_T *a;/* a in a.x = b */
  double *b;/* in ax=b */
  double *x;/* in ax=b */
  double **xs;/* x for series solving Ax[i]=b[i], i runs 0,...,ns */
  double **bs;/* b for series solving Ax[i]=b[i], i runs 0,...,ns */
  Uint ns;/* number of series */
  double refinement;/* max iter. refinement step, if DBL_MAX is set then it uses the default of UMFPACK which is 2 */
}Umfpack_T;

/* dense LAPACK direct solver for a series of b's (dgetrf + dgetrs) */
typedef struct LAPACK_DENSE_T
{
  char description[9000];
  struct MATRIX_T *a;/* a in a.x = b, square and in regular format */
  double **xs;/* x for series solving Ax[i]=b[i], i runs 0,...,ns */
  double **bs;/* b for series solving Ax[i]=b[i], i runs 0,...,ns */
  Uint ns;/* number of series */
}Lapack_Dense_T;

int direct_solver_umfpack_di(void *vp);
int direct_solver_umfpack_dl(void *vp);
int direct_solver_series_umfpack_di(void *vp);
int direct_solver_series_umfpack_dl(void *vp);
Umfpack_T *init_umfpack(void);
void free_umfpack(Umfpack_T *umf);
int direct_solver_series_lapack_dense(void *vp);
int direct_solver_umfpack_dl_reuse(void *vp,const char *const key,
                                   const double tol,const Uint max_iter);
void free_umfpack_reuse_cache(void);
Lapack_Dense_T *init_lapack_dense(void);
void free_lapack_dense(Lapack_Dense_T *lpk);

void solver_tests(void);

#endif


