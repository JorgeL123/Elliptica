/*
// Oliver Steppohn
// October 2026
*/

#include "lapack_method.h"

#define MSG_LEN (99)

/* dense LAPACK direct solver for a series of a.x[i] = b[i] with the same a.
// a must be square and in regular (reg) format, i.e. a->reg->A[row][col].
// a is factorized once (dgetrf) and all b's are solved in one call (dgetrs),
// so the solve phase runs as matrix-matrix (BLAS-3) operations instead of
// one triangular solve per b. a itself is not modified.
// ->return value: EXIT_SUCCESS
*/
int direct_solver_series_lapack_dense(void *vp)
{
  Lapack_Dense_T *const lpk = vp;
  const int n       = (int)lpk->a->row;
  const int ns      = (int)lpk->ns;/* number of series to be solved */
  double **const bs = lpk->bs;/* b's */
  double **const xs = lpk->xs;/* x's */
  double **const A  = lpk->a->reg->A;
  const double time1 = get_time_sec();
  double *LU = 0,*X = 0,*work = 0;
  double anorm = 0,rcond = 0;
  int *ipiv = 0,*iwork = 0;
  int info = 0;
  char err[MSG_LEN];

  if (!lpk->a->reg_f)
    Error0("The matrix a in a.x = b must be in regular format.\n");
  if (lpk->a->row != lpk->a->col)
    Error0("The matrix a in a.x = b must be square.\n");

  /* row major a == column major a^T: copy rows contiguously, factorize
  // a^T and solve with trans = 'T', i.e. (a^T)^T.x = a.x = b.
  // ||a^T||_1 = ||a||_inf is needed for the condition number. */
  LU = alloc_double((Uint)n*(Uint)n);
  for (int r = 0; r < n; ++r)
  {
    double rsum = 0;
    memcpy(&LU[(size_t)r*(size_t)n],A[r],(size_t)n*sizeof(*LU));
    for (int c = 0; c < n; ++c)
      rsum += fabs(A[r][c]);
    if (rsum > anorm)
      anorm = rsum;
  }

  ipiv = calloc((size_t)n,sizeof(*ipiv));
  IsNull(ipiv);
  dgetrf_(&n,&n,LU,&n,ipiv,&info);
  if (info < 0)
  {
    sprintf(err,"argument %d had an illegal value",-info);
    Errors("dgetrf: %s.\n",err);
  }
  if (info > 0)
  {
    sprintf(err,"U(%d,%d) = 0",info,info);
    Errors("dgetrf: matrix is exactly singular, %s.\n",err);
  }

  /* condition number estimate, O(n^2) */
  work  = alloc_double(4*(Uint)n);
  iwork = calloc((size_t)n,sizeof(*iwork));
  IsNull(iwork);
  dgecon_("1",&n,LU,&n,&anorm,&rcond,work,iwork,&info,1);
  free(work);
  free(iwork);

  /* all b's as one column major n x ns matrix */
  X = alloc_double((Uint)n*(Uint)ns);
  for (int i = 0; i < ns; ++i)
    memcpy(&X[(size_t)i*(size_t)n],bs[i],(size_t)n*sizeof(*X));

  dgetrs_("T",&n,&ns,LU,&n,ipiv,X,&n,&info,1);
  if (info < 0)
  {
    sprintf(err,"argument %d had an illegal value",-info);
    Errors("dgetrs: %s.\n",err);
  }

  for (int i = 0; i < ns; ++i)
    memcpy(xs[i],&X[(size_t)i*(size_t)n],(size_t)n*sizeof(*X));

  free(X);
  free(LU);
  free(ipiv);

  lpk->description[0] = '\0';
  sprintf(lpk->description,
            "o.  Series Solver     = LAPACK dgetrf/dgetrs (dense)\n"
            "o.  Matrix Dimension  = %dx%d\n"
            "o.  Number of RHS     = %d\n"
            "o.  Condition Number  = %e (1-norm estimate of a^T)\n"
            "o.  Solver Wall-Clock = %g(s)\n",
                                           n,n,ns,
                                           rcond > 0 ? 1/rcond : DBL_MAX,
                                           get_time_sec()-time1);

  return EXIT_SUCCESS;
}

/* initialize Lapack_Dense_T struct.
// ->return value: Lapack_Dense_T struct
*/
Lapack_Dense_T *init_lapack_dense(void)
{
  Lapack_Dense_T *lpk = calloc(1,sizeof(*lpk));
  IsNull(lpk);

  return lpk;
}

/* free Lapack_Dense_T struct; a, bs and xs belong to the caller. */
void free_lapack_dense(Lapack_Dense_T *lpk)
{
  if (!lpk)
    return;

  free(lpk);
}
