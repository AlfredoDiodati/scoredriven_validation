# Prompt for the et_al session

Copy everything below the line into a new session started in
/media/alfredo/programs3/_py/et_al.

---

Work on et_al only. Do not touch /media/alfredo/programs3/_py/ABM_collab/_C
beyond reading it; its changes come after this and in their own session.

Note on paths: my global CLAUDE.md gives the et_al development folder as
/home/dioda/Documents/_py/et_al., which does not exist on this machine. The
real one is /media/alfredo/programs3/_py/et_al. and the installed copy is
/usr/local/include/et_al./. Read the README policies and docs/ in full before
writing anything, per that same file.

## Why

Fitting t-QVARMA to a large batch of simulated series is about to go from
21,600 fits to 500,000, and two things in et_al make that far more expensive
than the arithmetic requires. Both were measured on 2026-08-29 in the ABM
project. Read the setup and the numbers there rather than re-running them:

    /media/alfredo/programs3/_py/ABM_collab/_C/out/fit_speedup_options.txt

The benchmarks behind it are that project's tests/qvarma_*.c and
tests/small_call_scaling.c, run by `make bench`. Every one of them raises
M_MMAP_THRESHOLD and M_TRIM_THRESHOLD through mallopt before timing, so the
numbers already exclude a separate 1.52x that was nothing but the kernel taking
back the tape's blocks on every tape_free.

The two findings that drive this work:

1. A 5x5 by 5x1 cblas_dgemm costs 158 ns on this machine, which is 50 flops at
   0.32 GFLOP/s, so the call overhead is the whole cost at this size. Under
   four concurrent threads the same call costs 1370 ns, because OpenBLAS's
   buffer table is one shared structure per process. malloc at the same rate
   scales 3.59x and pure arithmetic 3.74x, so neither the allocator nor the
   machine is the problem. The qvarma filter issues about 5,600 BLAS calls per
   value-and-gradient evaluation at T = 400, r = 2, which is 63 percent of its
   1.41 ms single-threaded cost and 78 percent of its 9.82 ms cost at four
   threads. Consequence: an OpenMP loop over independent fits runs slower than
   a serial one, 29.5 s against 52.9 s for the same 8 fits.

2. Hand-writing the filter's forward pass with no tape and no BLAS gives the
   same log-likelihood to every printed digit (-2599.66967) in 0.0362 ms
   against the taped 0.6867 ms, a 19.0x speedup, and it scales 3.51x on four
   threads where the taped path scales 0.58x.

## What to do, in this order

### 1. Small-matrix kernels in ad.h

ad_matmul calls mat_mul, hence cblas_dgemm, and its backward issues two more.
ad_chol_quadform calls _trtrs forward, and vec_chol_solve plus a malloc and
free of an n-length buffer backward, once per period. Below some dimension
threshold these should be inline kernels rather than BLAS calls, and the
per-call malloc in ad_chol_quadform_backward should be gone.

This is general: every score-driven model in et_al issues tiny gemms through
the tape, so this is not a qvarma fix. Measure the threshold rather than
picking one, and put the benchmark under tests/performance/ following the
naming policy. Report the crossover and the four-thread scaling before and
after.

### 2. A fused value-and-gradient filter in sd/qvarma.h

The forward recursion and its adjoint written by hand as one allocation-free
loop, with no tape.

Decisions already taken, so do not re-open them:

- It goes in sd/qvarma.h, beside _qvarma_filter. The traced variant stays;
  et_al's own model policy already requires a traced and an untraced variant
  that agree, with a test that checks it.
- Every dimension is read from QvarmaParams at runtime. The ABM pipeline fits
  only p = 1, q = 1, r = 2, K = 5, K_star = 3, R = 1, but nothing about that
  shape belongs in the code.
- The forward path the adjoint reads back (v_t, z_t, s_t, u_t, mu_star_t,
  mu_dag_t, so 5K + 1 doubles per period, about 83 KB at T = 400) is one
  allocation made when fit starts and reused across every evaluation. No
  allocation inside the loop.
- The link adjoints read the same name/transform/inverse/derivative table
  _qvarma_link and _qvarma_unlink already share, so the three cannot drift.
- qvarma_fit_cached and every other public signature stay as they are.

Derive the adjoint yourself. One part of it is where a hand derivation
normally goes wrong, so check it against the paper's own form: the quadratic
form solves Omega_inv z_t = v_t by forward substitution, and its adjoint needs
a *back* substitution against Omega_inv transposed, after which Omega_inv's
gradient takes a rank-one update that must be masked to the lower triangle.

### 3. Tests, all three required

- Fused value against _qvarma_filter's value, over random theta and random
  shapes, not only the ABM pipeline's shape.
- Fused gradient against the taped gradient componentwise, same sweep.
- Both against a central finite difference at a few theta. Without this a
  shared error in the two passes the second test.

Name them for the question they answer, per the naming policy.

## What to report

Per-evaluation cost of value-only and value-and-gradient, taped and fused, at
one and four threads, with the setup stated. I expect roughly 10x to 13x on the
value-and-gradient from a hand-written reverse pass costing two to four times
its forward, but that is an estimate and the measurement is the point.

Do not install over /usr/local until the tests pass.
