/*
// Oliver Steppohn
// October 2026
*/

#include "umfpack_reuse.h"

/* the cache lives across calls (and across outer iterations, where the
// grid is rebuilt), keyed by a name. */
#define UMF_REUSE_MAX_CACHE (32)
static Umf_Reuse_Cache_T *umf_reuse_cache[UMF_REUSE_MAX_CACHE];

/* solve a.x = b with UMFPACK, reusing an earlier factorization when
// possible. a is in CCS long format.
//
// if a factorization with the same key and dimension exists, x is found by
// GMRES on a, preconditioned with that (old) factorization. this converges
// in a few iterations when a changed only a little, e.g. between outer
// iterations. if GMRES does not reach the relative residual tol within
// max_iter iterations, or there is no usable factorization, a is factorized
// (as direct_solver_umfpack_dl does) and the factorization is cached.
// convergence is always checked against the current a, so an outdated
// factorization only costs time, never accuracy.
// ->return value: EXIT_SUCCESS
*/
int direct_solver_umfpack_dl_reuse(void *vp,const char *const key,
                                   const double tol,const Uint max_iter)
{
  Umfpack_T *const umf = vp;
  const long n     = (long)umf->a->row;
  long *const Ap   = umf->a->ccs_long->Ap;
  long *const Ai   = umf->a->ccs_long->Ai;
  double *const Ax = umf->a->ccs_long->Ax;
  double *const b  = umf->b;
  double *const x  = umf->x;
  Umf_Reuse_Cache_T *cache = 0;
  double Control[UMFPACK_CONTROL];
  double Info[UMFPACK_INFO];
  double res = 0;
  Uint niter = 0;
  long status;
  const double time1 = get_time_sec();
  const double Gb = 1.0E9;

  if (!umf->a->ccs_l_f)
    Error0("The matrix a in a.x = b in umfpack must be in CCS long format.\n");
  if (umf->a->row != umf->a->col)
    Error0("The matrix a in a.x = b must be square.\n");

  umfpack_dl_defaults(Control);
  if (umf->refinement != DBL_MAX)/* refinment step */
    Control[UMFPACK_IRSTEP] = umf->refinement;

  cache = find_cache(key);

  /* try the old factorization as preconditioner */
  if (cache->Numeric && cache->n == n)
  {
    if (gmres_preconditioned(umf->a,cache->Numeric,Control,b,x,
                             tol,max_iter,&niter,&res))
    {
      cache->nreuse++;
      umf->description[0] = '\0';
      sprintf(umf->description,
                "o.  Solver            = UMFPACK_dl, reused LU + GMRES\n"
                "o.  Matrix Dimension  = %ldx%ld\n"
                "o.  Reuse of LU       = %u (solves since factorization)\n"
                "o.  GMRES Iterations  = %u\n"
                "o.  Relative Residual = %e (tolerance %e)\n"
                "o.  Solver Wall-Clock = %g(s)\n",
                 n,n,cache->nreuse,niter,res,tol,get_time_sec()-time1);
      return EXIT_SUCCESS;
    }

    printf(Pretty0"%s: reused LU did not converge (%u iterations, "
           "rel. residual %e), refactorizing.\n",key,niter,res);
    fflush(stdout);
  }

  /* (re)factorize */
  if (cache->Numeric)
    umfpack_dl_free_numeric(&cache->Numeric);
  cache->Numeric = 0;

  void *Symbolic = 0;
  status = umfpack_dl_symbolic(n,n,Ap,Ai,Ax,&Symbolic,Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_dl(Control,status,__FILE__,__LINE__);

  status = umfpack_dl_numeric(Ap,Ai,Ax,Symbolic,&cache->Numeric,Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_dl(Control,status,__FILE__,__LINE__);

  umfpack_dl_free_symbolic(&Symbolic);
  cache->n = n;
  cache->nreuse = 0;

  status = umfpack_dl_solve(UMFPACK_A,Ap,Ai,Ax,x,b,cache->Numeric,
                            Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_dl(Control,status,__FILE__,__LINE__);

  umf->description[0] = '\0';
  sprintf(umf->description,
            "o.  Solver            = UMFPACK_dl, new LU (kept for reuse)\n"
            "o.  Matrix Dimension  = %ldx%ld\n"
            "o.  Condition Number  = %e\n"
            "o.  Sparsity          = %0.2f(%%)\n"
            "o.  Refinement Step   = %d\n"
            "o.  Peak Memory       = %g(Gbs)\n"
            "o.  Number of Flops   = %g\n"
            "o.  Number of Realloc = %g\n"
            "o.  Solver Wall-Clock = %g(s)\n",
             n,n,
             1/Info[UMFPACK_RCOND],
             (1-(double)Ap[n]/((double)n*(double)n))*100.,
             (int)Control[UMFPACK_IRSTEP],
             Info[UMFPACK_PEAK_MEMORY]*Info[UMFPACK_SIZE_OF_UNIT]/Gb,
             Info[UMFPACK_FLOPS],
             Info[UMFPACK_NUMERIC_REALLOC],
             get_time_sec()-time1);

  return EXIT_SUCCESS;
}

/* free all cached factorizations. */
void free_umfpack_reuse_cache(void)
{
  for (Uint i = 0; i < UMF_REUSE_MAX_CACHE && umf_reuse_cache[i]; ++i)
  {
    if (umf_reuse_cache[i]->Numeric)
      umfpack_dl_free_numeric(&umf_reuse_cache[i]->Numeric);
    free(umf_reuse_cache[i]);
    umf_reuse_cache[i] = 0;
  }
}

/* ->return value: the cache entry of key, a new empty one if none. */
static Umf_Reuse_Cache_T *find_cache(const char *const key)
{
  Uint i;

  for (i = 0; i < UMF_REUSE_MAX_CACHE && umf_reuse_cache[i]; ++i)
    if (!strcmp(umf_reuse_cache[i]->key,key))
      return umf_reuse_cache[i];

  if (i == UMF_REUSE_MAX_CACHE)
    Error0("Too many keys for the UMFPACK reuse cache.\n");

  umf_reuse_cache[i] = calloc(1,sizeof(*umf_reuse_cache[i]));
  IsNull(umf_reuse_cache[i]);
  snprintf(umf_reuse_cache[i]->key,sizeof(umf_reuse_cache[i]->key),
           "%s",key);

  return umf_reuse_cache[i];
}

/* y = a.x for a in CCS long format, threads work on column blocks and
// sum their partial results. */
static void ccs_long_matvec(const Matrix_T *const a,const double *const x,
                            double *const y)
{
  const long n = a->col;
  const long nr = a->row;
  const long *const Ap = a->ccs_long->Ap;
  const long *const Ai = a->ccs_long->Ai;
  const double *const Ax = a->ccs_long->Ax;

  memset(y,0,(size_t)nr*sizeof(*y));

  #pragma omp parallel
  {
    double *yp = calloc((size_t)nr,sizeof(*yp));
    IsNull(yp);

    #pragma omp for schedule(static)
    for (long c = 0; c < n; ++c)
    {
      const double xc = x[c];
      for (long k = Ap[c]; k < Ap[c+1]; ++k)
        yp[Ai[k]] += Ax[k]*xc;
    }

    #pragma omp critical (ccs_long_matvec_sum)
    {
      for (long r = 0; r < nr; ++r)
        y[r] += yp[r];
    }
    free(yp);
  }
}

static double norm2(const double *const v,const long n)
{
  double s = 0;
  for (long i = 0; i < n; ++i)
    s += v[i]*v[i];
  return sqrt(s);
}

/* one GMRES cycle on a.dx = r, right preconditioned with the UMFPACK
// factorization Numeric (M^-1 = Numeric solve), starting from dx = 0.
// it stops when the estimated residual is <= stop*nb or after m
// iterations. ->return value: number of iterations done. */
static Uint gmres_cycle(const Matrix_T *const a,void *const Numeric,
                        double *const Control,
                        const double *const r,const double nr,
                        const double stop,const Uint m,double *const dx)
{
  const long n = a->row;
  const long *const Ap = a->ccs_long->Ap;
  const long *const Ai = a->ccs_long->Ai;
  const double *const Ax = a->ccs_long->Ax;
  double Info[UMFPACK_INFO];
  double *V   = alloc_double((Uint)n*(m+1));/* Krylov basis */
  double *H   = alloc_double((m+1)*m);/* Hessenberg, H[row*m+col] */
  double *cs  = alloc_double(m);
  double *sn  = alloc_double(m);
  double *g   = alloc_double(m+1);
  double *yk  = alloc_double(m);
  double *w   = alloc_double((Uint)n);
  double *z   = alloc_double((Uint)n);
  Uint j = 0;
  long i;
  
  for (i = 0; i < n; ++i)
    V[i] = r[i]/nr;
  g[0] = nr;
  
  for (j = 0; j < m; ++j)
  {
    double h;
    
    /* w = a.M^-1.v_j */
    if (umfpack_dl_solve(UMFPACK_A,Ap,Ai,Ax,z,&V[(long)j*n],Numeric,
                         Control,Info) != UMFPACK_OK)
      break;
    ccs_long_matvec(a,z,w);
    
    /* modified Gram-Schmidt */
    for (Uint k = 0; k <= j; ++k)
    {
      const double *const vk = &V[(long)k*n];
      h = 0;
      for (i = 0; i < n; ++i)
        h += w[i]*vk[i];
      H[k*m+j] = h;
      for (i = 0; i < n; ++i)
        w[i] -= h*vk[i];
    }
    h = norm2(w,n);
    H[(j+1)*m+j] = h;
    if (h != 0)
      for (i = 0; i < n; ++i)
        V[(long)(j+1)*n+i] = w[i]/h;
    
    /* apply previous Givens rotations and make a new one */
    for (Uint k = 0; k < j; ++k)
    {
      const double h1 = H[k*m+j],h2 = H[(k+1)*m+j];
      H[k*m+j]     =  cs[k]*h1+sn[k]*h2;
      H[(k+1)*m+j] = -sn[k]*h1+cs[k]*h2;
    }
    {
      const double h1 = H[j*m+j],h2 = H[(j+1)*m+j];
      const double d  = sqrt(h1*h1+h2*h2);
      cs[j] = d == 0 ? 1 : h1/d;
      sn[j] = d == 0 ? 0 : h2/d;
      H[j*m+j]     = d;
      H[(j+1)*m+j] = 0;
      g[j+1] = -sn[j]*g[j];
      g[j]   =  cs[j]*g[j];
    }
    
    if (fabs(g[j+1]) <= stop || h == 0)
    {
      ++j;
      break;
    }
  }
  
  /* dx = M^-1.V.y, H.y = g */
  for (long k = (long)j-1; k >= 0; --k)
  {
    yk[k] = g[k];
    for (Uint l = (Uint)k+1; l < j; ++l)
      yk[k] -= H[k*m+l]*yk[l];
    yk[k] /= H[k*m+k];
  }
  memset(w,0,(size_t)n*sizeof(*w));
  for (Uint k = 0; k < j; ++k)
    for (i = 0; i < n; ++i)
      w[i] += yk[k]*V[(long)k*n+i];
  umfpack_dl_solve(UMFPACK_A,Ap,Ai,Ax,dx,w,Numeric,Control,Info);
  
  free(V);free(H);free(cs);free(sn);free(g);free(yk);free(w);free(z);
  
  return j;
}

/* restarted GMRES on a.x = b, right preconditioned with the UMFPACK
// factorization Numeric, starting from x = 0. after each cycle the true
// residual b - a.x is computed; if it is still above tol, the next cycle
// solves for a correction (like iterative refinement). at most max_iter
// iterations in total.
// ->return value: 1 if ||b - a.x||/||b|| <= tol, 0 otherwise. */
static int gmres_preconditioned(const Matrix_T *const a,void *const Numeric,
                                double *const Control,
                                const double *const b,double *const x,
                                const double tol,const Uint max_iter,
                                Uint *const niter,double *const res)
{
  const long n = a->row;
  double *r  = alloc_double((Uint)n);
  double *dx = alloc_double((Uint)n);
  const double nb = norm2(b,n);
  double nr = nb;
  long i;
  
  memset(x,0,(size_t)n*sizeof(*x));
  memcpy(r,b,(size_t)n*sizeof(*r));
  *niter = 0;
  *res   = 0;
  
  while (nb > 0 && nr > tol*nb && *niter < max_iter)
  {
    Uint left = max_iter - *niter;
    Uint m    = left < UMF_REUSE_MAX_KRYLOV ? left : UMF_REUSE_MAX_KRYLOV;
    
    /* stop a bit below tol since the estimate can be optimistic */
    Uint j = gmres_cycle(a,Numeric,Control,r,nr,0.1*tol*nb,m,dx);
    
    *niter += j;
    for (i = 0; i < n; ++i)
      x[i] += dx[i];
    
    /* true residual */
    ccs_long_matvec(a,x,r);
    for (i = 0; i < n; ++i)
      r[i] = b[i]-r[i];
    nr = norm2(r,n);
    
    if (!j)
      break;
  }
  *res = nb > 0 ? nr/nb : 0;
  
  free(r);
  free(dx);
  
  return *res <= tol;
}
