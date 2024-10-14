#include <complex.h>
#include <ekf/ekf.hpp>
#include <time.h>

/* 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise */
static int is_sym(cs_ci *A) {
  int is_upper, is_lower, j, p, n = A->n, m = A->m, *Ap = A->p, *Ai = A->i;
  if (m != n)
    return (0);
  is_upper = 1;
  is_lower = 1;
  for (j = 0; j < n; j++) {
    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      if (Ai[p] > j)
        is_upper = 0;
      if (Ai[p] < j)
        is_lower = 0;
    }
  }
  return (is_upper ? 1 : (is_lower ? -1 : 0));
}

/* true for off-diagonal entries */
static int dropdiag(int i, int j, cs_complex_t aij, void *other) {
  return (i != j);
}

/* C = A + triu(A,1)' */
static cs_ci *make_sym(cs_ci *A) {
  cs_ci *AT, *C;
  AT = cs_ci_transpose(A, 1);       /* AT = A' */
  cs_ci_fkeep(AT, &dropdiag, NULL); /* drop diagonal entries from AT */
  C = cs_ci_add(A, AT, 1, 1);       /* C = A+AT */
  cs_ci_spfree(AT);
  return (C);
}

/* create a right-hand side */
static void rhs(cs_complex_t *x, cs_complex_t *b, int m) {
  int i;
  for (i = 0; i < m; i++)
    b[i] = 1 + ((double)i) / m;
  for (i = 0; i < m; i++)
    x[i] = b[i];
}

/* infinity-norm of x */
static double norm(cs_complex_t *x, int n) {
  int i;
  double normx = 0;
  for (i = 0; i < n; i++)
    normx = CS_MAX(normx, cabs(x[i]));
  return (normx);
}

/* compute residual, norm(A*x-b,inf) / (norm(A,1)*norm(x,inf) + norm(b,inf)) */
static void print_resid(int ok, cs_ci *A, cs_complex_t *x, cs_complex_t *b,
                        cs_complex_t *resid) {
  int i, m, n;
  if (!ok) {
    printf("    (failed)\n");
    return;
  }
  m = A->m;
  n = A->n;
  for (i = 0; i < m; i++)
    resid[i] = -b[i];       /* resid = -b */
  cs_ci_gaxpy(A, x, resid); /* resid = resid + A*x  */
  printf("resid: %8.2e\n",
         norm(resid, m) /
             ((n == 0) ? 1 : (cs_ci_norm(A) * norm(x, n) + norm(b, m))));
}

static double tic(void) { return (clock() / (double)CLOCKS_PER_SEC); }

static double toc(double t) {
  double s = tic();
  return (CS_MAX(0, s - t));
}

static void print_order(int order) {
  switch (order) {
  case 0:
    printf("natural    ");
    break;
  case 1:
    printf("amd(A+A')  ");
    break;
  case 2:
    printf("amd(S'*S)  ");
    break;
  case 3:
    printf("amd(A'*A)  ");
    break;
  }
}

int solve(Problem *problem) {
  cs_ci *A, *C;
  cs_complex_t *b, *x, *resid;
  double t, tol;
  int k, m, n, ok, order, nb, ns, *r, *s, *rr, sprank;
  cs_cid *D;
  if (!problem)
    return (0);
  A = problem->A;
  C = problem->C;
  b = problem->b;
  x = problem->x;
  resid = problem->resid;
  m = A->m;
  n = A->n;
  tol = problem->sym ? 0.001 : 1; /* partial pivoting tolerance */
  D = cs_ci_dmperm(C, 1);         /* randomized dmperm analysis */
  if (!D)
    return (0);
  nb = D->nb;
  r = D->r;
  s = D->s;
  rr = D->rr;
  sprank = rr[3];
  for (ns = 0, k = 0; k < nb; k++) {
    ns += ((r[k + 1] == r[k] + 1) && (s[k + 1] == s[k] + 1));
  }
  printf("blocks: %g singletons: %g structural rank: %g\n", (double)nb,
         (double)ns, (double)sprank);
  cs_ci_dfree(D);
  for (order = 0; order <= 3; order += 3) /* natural and amd(A'*A) */
  {
    if (!order && m > 1000)
      continue;
    printf("QR   ");
    print_order(order);
    rhs(x, b, m); /* compute right-hand side */
    t = tic();
    ok = cs_ci_qrsol(order, C, x); /* min norm(Ax-b) with QR */
    printf("time: %8.2f ", toc(t));
    print_resid(ok, C, x, b, resid); /* print residual */
  }
  if (m != n || sprank < n)
    return (1);                        /* return if rect. or singular*/
  for (order = 0; order <= 3; order++) /* try all orderings */
  {
    if (!order && m > 1000)
      continue;
    printf("LU   ");
    print_order(order);
    rhs(x, b, m); /* compute right-hand side */
    t = tic();
    ok = cs_ci_lusol(order, C, x, tol); /* solve Ax=b with LU */
    printf("time: %8.2f ", toc(t));
    print_resid(ok, C, x, b, resid); /* print residual */
  }
  if (!problem->sym)
    return (1);
  for (order = 0; order <= 1; order++) /* natural and amd(A+A') */
  {
    if (!order && m > 1000)
      continue;
    printf("Chol ");
    print_order(order);
    rhs(x, b, m); /* compute right-hand side */
    t = tic();
    ok = cs_ci_cholsol(order, C, x); /* solve Ax=b with Cholesky */
    printf("time: %8.2f ", toc(t));
    print_resid(ok, C, x, b, resid); /* print residual */
  }
  return (1);
}
