Welcome to

```ascii
EEEEEEEEEEEE  ll  ll
EE            ll  ll
EE            ll  ll  O                       O
EE            ll  ll                   tt
EEEEEEEEEEEE  ll  ll  ii  pppppppp. tttttttt  ii    .ccccc   .aaaaa.
EE            ll  ll  ii  pp      p.   tt     ii   .c       .a     a.
EE            ll  ll  ii  pp       p.  tt     ii  .c       .a       a.
EE            ll  ll  ii  pp      p.   tt     ii   .c       .a     .a
EEEEEEEEEEEE  ll  ll  ii  pppppppp.    tt     ii    .ccccc   .aaaaa.aa
                          pp                                          aa
                          pp
                          pp
                          pp

```

# Installation

Requirements:
- GNU make.
- Git version control system.
- A C programming compiler, for instance, `gcc`.
- A compatible `openmp` library with the compiler.
- The `UMFPACK` library (`suitesparse`).
- A `LAPACK` library with 32-bit integers (LP64), for instance `OpenBLAS`.

```bash
git clone --depth=1 git@github.com:rashti-alireza/Elliptica.git
cd Elliptica
make MyConfig
make git_clone
make -j4
```

You can see the list of available projects after issuing `make git_clone` in `Src/Projects`.

## Linking BLAS/LAPACK

This mirror of Elliptica modifies the Schur complement solver factorization of the (dense) matrix of each patch with
LAPACK (`dgetrf`/`dgetrs`), which results in a tremendous speed-up compared to the publicly available version of Elliptica. So `MyConfig` must link a LAPACK, e.g.

```make
SPECIAL_LIBS += -lopenblas
```

See `Doc/Make/MyConfig.example` for an example. The example also lists alternatives (reference LAPACK, MKL).
Things to watch out for:

- Link an LP64 build (32-bit integers). ILP64 builds such as
  `libopenblas64.so` or `mkl_intel_ilp64` compile and link fine but
  give wrong results or crashes.
- Link only one BLAS/LAPACK provider. If several are linked (e.g. OpenBLAS
  and MKL), the link order decides which one is used.
- OpenBLAS comes in serial, pthreads and OpenMP flavors. Prefer the OpenMP
  (`libopenblaso`) or serial one; the pthreads one starts its own threads
  inside Elliptica's OpenMP loops.
- Older OpenBLAS versions (e.g. 0.3.15) misdetect newer CPUs (e.g. Ice Lake)
  and run several times slower; `export OPENBLAS_CORETYPE=SkylakeX` fixes it.
- Memory: each patch needs about `8*(NS^2 + NS*NI)` bytes for the dense
  solve (NS, NI: subdomain and interface points of the patch); all patches
  are solved concurrently. `solve_ddm_schur_thread_cap` limits this.

Relevant solver parameters (defaults first):

- `solve_ddm_schur_patch_solver = lapack_dense | umfpack`
- `solve_ddm_schur_S_reuse = yes | no` reuses the LU of the Schur
  complement S of the previous outer iteration as GMRES preconditioner
  (keeps one LU per field in memory); see also
  `solve_ddm_schur_S_reuse_tolerance` (1E-12) and
  `solve_ddm_schur_S_reuse_max_iter` (15), above which S is refactorized.
- `solve_ddm_schur_jacobian_loop = stencil | full` and
  `solve_ddm_schur_jacobian_check = no | yes` (compares both loops).

When restarting from a checkpoint, parameters are taken from the checkpoint;
use `modify:<parameter> = <value>` to change them.

# Run

In each project there is a directory named `XXXX_ParFiles`, where `XXXX` denotes the project name, 
e.g., `NSNS_ParFiles`.
One can find various parameter files for different gravitational systems. 
Having set the parameters of interest, one invokes the code:

`$ ./Exe/elliptica /path/to/the/parameter/file`


# Citation

inspirehep.net citation key = `{Rashti:2021ihv,Rashti:2024drr}`


