#include "test_utils.h"
#include "gtest/gtest.h"

#include "qoco.h"

TEST(kkt_test, qoco_kkt_multiply_matches_kkt_multiply)
{
  // Small problem: minimize (1/2) x' I x  s.t.  G x <= h (non-negative cone)
  QOCOInt n = 2;
  QOCOInt p = 0;
  QOCOInt m = 1;
  QOCOInt l = 1;      // non-neg orthant dimension
  QOCOInt nsoc = 0;   // no SOCs

  // P = I_2 in CSC
  QOCOFloat Px[] = {1.0, 1.0};
  QOCOInt Pnnz = 2;
  QOCOInt Pp[] = {0, 1, 2};
  QOCOInt Pi[] = {0, 1};

  // No equality constraints (p = 0)

  // G = [1 0]
  QOCOFloat Gx[] = {1.0};
  QOCOInt Gnnz = 1;
  QOCOInt Gp[] = {0, 1, 1};
  QOCOInt Gi[] = {0};

  QOCOFloat c[] = {0.0, 0.0};
  QOCOFloat h[] = {1.0};
  QOCOInt* q = nullptr;

  QOCOCscMatrix* P = (QOCOCscMatrix*)malloc(sizeof(QOCOCscMatrix));
  QOCOCscMatrix* G = (QOCOCscMatrix*)malloc(sizeof(QOCOCscMatrix));

  qoco_set_csc(P, n, n, Pnnz, Px, Pp, Pi);
  qoco_set_csc(G, m, n, Gnnz, Gx, Gp, Gi);

  QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
  QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

  set_default_settings(settings);

  QOCOInt exit = qoco_setup(solver, n, m, p, P, c, nullptr, nullptr, G, h,
                            l, nsoc, q, settings);
  ASSERT_EQ(exit, QOCO_NO_ERROR);

  // Solve once to ensure NT scaling (Wfull) is computed.
  exit = qoco_solve(solver);
  ASSERT_TRUE(exit == QOCO_SOLVED || exit == QOCO_SOLVED_INACCURATE);

  // Test vector (length n + p + m)
  QOCOInt len = n + p + m;
  QOCOFloat x[] = {1.0, 2.0, 3.0};
  QOCOFloat y_api[3] = {0.0, 0.0, 0.0};
  QOCOFloat y_ref[3] = {0.0, 0.0, 0.0};

  // Call new public API.
  exit = qoco_kkt_multiply(solver, x, y_api);
  ASSERT_EQ(exit, QOCO_NO_ERROR);

  // Compute reference result using internal kkt_multiply directly.
  QOCOWorkspace* work = solver->work;
  kkt_multiply(x, y_ref, work->data, get_data_vectorf(work->Wfull),
               get_data_vectorf(work->xbuff), get_data_vectorf(work->ubuff1),
               get_data_vectorf(work->ubuff2));

  // Compare results.
  expect_eq_vectorf(y_api, y_ref, len, 1e-8);

  // Cleanup.
  qoco_cleanup(solver);
  free(settings);
  free(P);
  free(G);
}
