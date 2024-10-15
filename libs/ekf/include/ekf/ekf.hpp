#include <cs.h>
#include <stdio.h>

struct Problem {
  Problem();
  ~Problem();
  cs_ci *A;
  cs_ci *C;
  int sym;
  cs_complex_t *x;
  cs_complex_t *b;
  cs_complex_t *resid;
};

Problem *getProblem(FILE *file, float tol);
Problem *freeProblem(Problem *problem);
int solve(Problem *problem);
