#include <cs.h>

struct Problem {
  Problem() {};
  ~Problem() {
    cs_ci_spfree(A);
    cs_ci_spfree(C);
    cs_ci_free(x);
    cs_ci_free(b);
    cs_ci_free(resid);
  };

  cs_ci *A;
  cs_ci *C;
  int sym;
  cs_complex_t *x;
  cs_complex_t *b;
  cs_complex_t *resid;
};

int solve(Problem *problem);
