#include <stdlib.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/least_squares.h"

// Random number generator.
#define RNG_SIZE 256
static char rng_state[RNG_SIZE];

// LAPACK prototypes.
void dgetrf(int *N, int *NRHS, double *A, int *LDA, int *IPIV, int *INFO);
void dgetrs(char *TRANS, int *N, int *NRHS, double *A, 
            int *LDA, int *IPIV, double *B, int *LDB, int *INFO);

void test_poly_fit(int p, point_t* x0, point_t* points, int num_points, double* coeffs)
{
  int dim = poly_ls_basis_size(p);

  // Create scatter data.
  double data[num_points], basis[dim];
  memset(data, 0, sizeof(double)*num_points);
  for (int i = 0; i < num_points; ++i)
  {
    compute_poly_ls_basis_vector(p, &points[i], basis);
    for (int k = 0; k < dim; ++k)
      data[i] += coeffs[k]*basis[k];
  }

  // Assemble the linear system.
  double A[dim*dim], b[dim];
  compute_poly_ls_system(p, x0, points, num_points, data, A, b);

  // Solve the system.
  int lda = dim, ldb = dim, pivot[dim], info, one = 1;
  char trans = 'N';
  dgetrf(&dim, &dim, A, &lda, pivot, &info);
  assert_int_equal(0, info);
  dgetrs(&trans, &dim, &one, A, &lda, pivot, b, &ldb, &info);
  assert_int_equal(0, info);

  // Check the coefficients of the fit.
  for (int i = 0; i < dim; ++i)
  {
    assert_true(fabs(b[i] - coeffs[i]) < 1e-12);
  }
}

static void generate_random_points(int num_points, point_t* points)
{
  for (int i = 0; i < num_points; ++i)
  {
    points[i].x = 1.0*random()/RAND_MAX;
    points[i].y = 1.0*random()/RAND_MAX;
    points[i].z = 1.0*random()/RAND_MAX;
  }
}

void test_p0(void** state)
{
  static double coeffs[] = {1.0};
  point_t points[4];
  generate_random_points(4, points);
  test_poly_fit(0, NULL, points, 4, coeffs);
}

void test_p1(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t points[8];
  generate_random_points(8, points);
  test_poly_fit(1, NULL, points, 8, coeffs);
}

void test_p2(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t points[16];
  generate_random_points(16, points);
  test_poly_fit(2, NULL, points, 16, coeffs);
}

void test_p3(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t points[30];
  generate_random_points(30, points);
  test_poly_fit(3, NULL, points, 30, coeffs);
}

int main(int argc, char* argv[]) 
{
  arbi_init(argc, argv);

  // Initialize the random number generator.
  unsigned seed = 1;
  initstate(seed, rng_state, RNG_SIZE);

  const UnitTest tests[] = 
  {
    unit_test(test_p0),
    unit_test(test_p1),
    unit_test(test_p2)
  };
  return run_tests(tests);
}
