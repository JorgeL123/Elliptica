#include "core_lib.h"
#include "error_handling_lib.h"
#include <umfpack.h>
#include "maths_linear_algebra_lib.h"
#include "maths_matrix_solvers_lib.h"

/* max Krylov dimension of one GMRES cycle */
#define UMF_REUSE_MAX_KRYLOV (100)

/* a cached UMFPACK factorization, one per key (e.g. field name) */
typedef struct UMF_REUSE_CACHE_T
{
  char key[256];
  long n;/* dimension of the factorized matrix */
  void *Numeric;/* UMFPACK numeric factorization */
  Uint nreuse;/* number of solves since the factorization */
}Umf_Reuse_Cache_T;

int direct_solver_umfpack_dl_reuse(void *vp,const char *const key,
                                   const double tol,const Uint max_iter);
void free_umfpack_reuse_cache(void);
void umfpack_error_dl(const double *const Control,const long status,
                      const char *const file,const int line);
static Umf_Reuse_Cache_T *find_cache(const char *const key);
static void ccs_long_matvec(const Matrix_T *const a,const double *const x,
                            double *const y);
static int gmres_preconditioned(const Matrix_T *const a,void *const Numeric,
                                double *const Control,
                                const double *const b,double *const x,
                                const double tol,const Uint max_iter,
                                Uint *const niter,double *const res);
static double norm2(const double *const v,const long n);
static Uint gmres_cycle(const Matrix_T *const a,void *const Numeric,
                        double *const Control,
                        const double *const r,const double nr,
                        const double stop,const Uint m,double *const dx);
