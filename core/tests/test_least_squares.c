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

static void generate_random_points(int num_points, point_t* points)
{
  for (int i = 0; i < num_points; ++i)
  {
    points[i].x = 1.0*random()/RAND_MAX;
    points[i].y = 1.0*random()/RAND_MAX;
    points[i].z = 1.0*random()/RAND_MAX;
  }
}

static void average_points(point_t* points, int num_points, point_t* average)
{
  average->x = 0.0;
  average->y = 0.0;
  average->z = 0.0;
  for (int i = 0; i < num_points; ++i)
  {
    average->x += points[i].x;
    average->y += points[i].y;
    average->z += points[i].z;
  }
  average->x /= num_points;
  average->y /= num_points;
  average->z /= num_points;
}

// LAPACK prototypes.
void dgetrf(int *N, int *NRHS, double *A, int *LDA, int *IPIV, int *INFO);
void dgetrs(char *TRANS, int *N, int *NRHS, double *A, 
            int *LDA, int *IPIV, double *B, int *LDB, int *INFO);

// Weighting function: W(x, x0) = 1/(|x-x0|^2 + eps^2)
static const double epsilon = 1e-4;
static void weighting_func(void* context, point_t* x, point_t* x0, double* W, vector_t* gradient)
{
  double D = (x0 != NULL) ? point_distance(x, x0) : sqrt(x->x*x->x+x->y*x->y+x->z*x->z);
  *W = 1.0 / (pow(D, 2) + pow(epsilon, 2));
  double dDdx = x->x / D, dDdy = x->y / D, dDdz = x->z / D;
  double deriv_term = -(*W)*(*W) * 2.0 * D;
  gradient->x = deriv_term * dDdx;
  gradient->y = deriv_term * dDdy;
  gradient->z = deriv_term * dDdz;
}

void test_poly_fit(int p, point_t* x0, point_t* points, int num_points, double* coeffs, bool weighted)
{
  int dim = poly_ls_basis_size(p);

  // Create scatter data.
  double data[num_points], basis[dim];
  memset(data, 0, sizeof(double)*num_points);
  for (int i = 0; i < num_points; ++i)
  {
    if (x0 != NULL)
    {
      point_t y = {.x = points[i].x - x0->x, 
                   .y = points[i].y - x0->y,
                   .z = points[i].z - x0->z};
      compute_poly_ls_basis_vector(p, &y, basis);
    }
    else
      compute_poly_ls_basis_vector(p, &points[i], basis);
    for (int k = 0; k < dim; ++k)
      data[i] += coeffs[k]*basis[k];
  }

  // Assemble the linear system.
  double A[dim*dim], b[dim];
  if (weighted)
    compute_weighted_poly_ls_system(p, weighting_func, x0, points, num_points, data, A, b);
  else
    compute_poly_ls_system(p, x0, points, num_points, data, A, b);

  // Solve the system.
  int lda = dim, ldb = dim, pivot[dim], info, one = 1;
  char trans = 'N';
  dgetrf(&dim, &dim, A, &lda, pivot, &info);
  assert_int_equal(0, info);
  dgetrs(&trans, &dim, &one, A, &lda, pivot, b, &ldb, &info);
  assert_int_equal(0, info);

  // Check the coefficients of the fit.
  static const double tolerances[] = {1e-15, 2e-14, 2e-12, 1e-9};
  for (int i = 0; i < dim; ++i)
  {
    assert_true(fabs(b[i] - coeffs[i]) < tolerances[p]);
  }
}

void test_poly_shape_functions(int p, point_t* x0, point_t* points, int num_points, double* coeffs, bool weighted)
{
  int dim = poly_ls_basis_size(p);

  // Create scatter data.
  double data[num_points], basis[dim];
  memset(data, 0, sizeof(double)*num_points);
  for (int i = 0; i < num_points; ++i)
  {
    point_t y = {.x = points[i].x - x0->x, 
                 .y = points[i].y - x0->y,
                 .z = points[i].z - x0->z};
    compute_poly_ls_basis_vector(p, &y, basis);
    for (int k = 0; k < dim; ++k)
      data[i] += coeffs[k]*basis[k];
  }

  // Compute shape functions for the given data.
  poly_ls_shape_t* N = poly_ls_shape_new(p, false);
  if (weighted)
    poly_ls_shape_set_simple_weighting_func(N, 2.0, 1e-8);
  poly_ls_shape_set_domain(N, x0, points, num_points);
  double Nk[num_points];

  // Make sure the shape functions interpolate the data.
  for (int i = 0; i < num_points; ++i)
  {
    double value = 0.0;
    poly_ls_shape_compute(N, &points[i], Nk);
    for (int k = 0; k < num_points; ++k)
      value += Nk[k]*data[k];
    assert_true(fabs(value - data[i]) < 1e-12);
  }

  // Now make sure that the fit matches the polynomial at another point.
  point_t point;
  generate_random_points(1, &point);
  poly_ls_shape_compute(N, &point, Nk);
  double phi = 0.0, phi_fit = 0.0;
  point_t z = {.x = point.x - x0->x, 
               .y = point.y - x0->y,
               .z = point.z - x0->z};
  compute_poly_ls_basis_vector(p, &z, basis);
  for (int k = 0; k < dim; ++k)
    phi += coeffs[k] * basis[k];
  for (int k = 0; k < num_points; ++k)
    phi_fit += Nk[k] * data[k];
  assert_true(fabs(phi_fit - phi) < 1e-12);

  // Clean up.
  N = NULL;
}

void test_poly_shape_function_gradients(int p, point_t* x0, point_t* points, int num_points, double* coeffs, bool weighted)
{
  int dim = poly_ls_basis_size(p);

  // Create scatter data.
  double data[num_points], basis[dim];
  vector_t data_grads[num_points], basis_grad[dim];
  memset(data, 0, sizeof(double)*num_points);
  memset(data_grads, 0, sizeof(vector_t)*num_points);
  for (int i = 0; i < num_points; ++i)
  {
    point_t y = {.x = points[i].x - x0->x, 
                 .y = points[i].y - x0->y,
                 .z = points[i].z - x0->z};
    compute_poly_ls_basis_vector(p, &y, basis);
    compute_poly_ls_basis_gradient(p, &y, basis_grad);
    for (int k = 0; k < dim; ++k)
    {
      data[i] += coeffs[k]*basis[k];
      data_grads[i].x += coeffs[k]*basis_grad[k].x;
      data_grads[i].y += coeffs[k]*basis_grad[k].y;
      data_grads[i].z += coeffs[k]*basis_grad[k].z;
    }
  }

  // Compute shape functions for the given data.
  poly_ls_shape_t* N = poly_ls_shape_new(p, true);
  if (weighted)
    poly_ls_shape_set_simple_weighting_func(N, 2.0, 1e-8);
  poly_ls_shape_set_domain(N, x0, points, num_points);

  // Make sure the shape functions interpolate the data and their gradients.
  double Nk[num_points];
  vector_t gradNk[num_points];
  for (int i = 0; i < num_points; ++i)
  {
    double value = 0.0;
    vector_t gradient = {.x = 0.0, .y = 0.0, .z = 0.0};
    poly_ls_shape_compute_gradients(N, &points[i], Nk, gradNk);
    for (int k = 0; k < num_points; ++k)
    {
      value += Nk[k]*data[k];
      gradient.x += gradNk[k].x*data[k];
      gradient.y += gradNk[k].y*data[k];
      gradient.z += gradNk[k].z*data[k];
    }
    assert_true(fabs(value - data[i]) < 1e-12);
    assert_true(fabs(gradient.x - data_grads[i].x) < 1e-12);
    assert_true(fabs(gradient.y - data_grads[i].y) < 1e-12);
    assert_true(fabs(gradient.z - data_grads[i].z) < 1e-12);
  }

  // Now make sure that the fit matches the polynomial and its 
  // derivative at another point.
  point_t point;
  generate_random_points(1, &point);
  poly_ls_shape_compute(N, &point, Nk);
  poly_ls_shape_compute_gradients(N, &point, Nk, gradNk);
  double phi = 0.0, phi_fit = 0.0;
  vector_t grad_phi = {.x = 0.0, .y = 0.0, .z = 0.0},
           grad_phi_fit = {.x = 0.0, .y = 0.0, .z = 0.0};
  point_t z = {.x = point.x - x0->x, 
               .y = point.y - x0->y,
               .z = point.z - x0->z};
  compute_poly_ls_basis_vector(p, &z, basis);
  compute_poly_ls_basis_gradient(p, &z, basis_grad);
  for (int k = 0; k < dim; ++k)
  {
    phi += coeffs[k] * basis[k];
    grad_phi.x += coeffs[k] * basis_grad[k].x;
    grad_phi.y += coeffs[k] * basis_grad[k].y;
    grad_phi.z += coeffs[k] * basis_grad[k].z;
  }
  for (int k = 0; k < num_points; ++k)
  {
    phi_fit += Nk[k] * data[k];
    grad_phi_fit.x += gradNk[k].x * data[k];
    grad_phi_fit.y += gradNk[k].y * data[k];
    grad_phi_fit.z += gradNk[k].z * data[k];
  }
  assert_true(fabs(phi_fit - phi) < 1e-12);
  assert_true(fabs(grad_phi_fit.x - grad_phi.x) < 1e-12);
  assert_true(fabs(grad_phi_fit.y - grad_phi.y) < 1e-12);
  assert_true(fabs(grad_phi_fit.z - grad_phi.z) < 1e-12);

  // Clean up.
  N = NULL;
}

void test_p0_fit(void** state)
{
  static double coeffs[] = {1.0};
  point_t x0, points[4];
  generate_random_points(4, points);
  test_poly_fit(0, NULL, points, 4, coeffs, false);
  average_points(points, 4, &x0);
  test_poly_fit(0, &x0, points, 4, coeffs, false);
}

void test_weighted_p0_fit(void** state)
{
  static double coeffs[] = {1.0};
  point_t x0, points[4];
  generate_random_points(4, points);
  test_poly_fit(0, NULL, points, 4, coeffs, true);
  average_points(points, 4, &x0);
  test_poly_fit(0, &x0, points, 4, coeffs, true);
}

void test_p0_shape_funcs(void** state)
{
  static double coeffs[] = {1.0};
  point_t x0, points[4];
  generate_random_points(4, points);
  average_points(points, 4, &x0);
  test_poly_shape_functions(0, &x0, points, 4, coeffs, false);
}

void test_weighted_p0_shape_funcs(void** state)
{
  static double coeffs[] = {1.0};
  point_t x0, points[4];
  generate_random_points(4, points);
  average_points(points, 4, &x0);
  test_poly_shape_functions(0, &x0, points, 4, coeffs, true);
}

void test_p0_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0};
  point_t x0, points[4];
  generate_random_points(4, points);
  average_points(points, 4, &x0);
  test_poly_shape_function_gradients(0, &x0, points, 4, coeffs, false);
}

void test_weighted_p0_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0};
  point_t x0, points[4];
  generate_random_points(4, points);
  average_points(points, 4, &x0);
  test_poly_shape_function_gradients(0, &x0, points, 4, coeffs, true);
}

void test_p1_fit(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0, points[8];
  generate_random_points(8, points);
  test_poly_fit(1, NULL, points, 8, coeffs, false);
  average_points(points, 8, &x0);
  test_poly_fit(1, &x0, points, 8, coeffs, false);
}

void test_weighted_p1_fit(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0, points[8];
  generate_random_points(8, points);
  test_poly_fit(1, NULL, points, 8, coeffs, true);
  average_points(points, 8, &x0);
  test_poly_fit(1, &x0, points, 8, coeffs, true);
}

void test_p1_shape_funcs(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0, points[8];
  generate_random_points(8, points);
  average_points(points, 8, &x0);
  test_poly_shape_functions(1, &x0, points, 8, coeffs, false);
}

void test_weighted_p1_shape_funcs(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0, points[8];
  generate_random_points(8, points);
  average_points(points, 8, &x0);
  test_poly_shape_functions(1, &x0, points, 8, coeffs, true);
}

void test_p1_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0, points[8];
  generate_random_points(8, points);
  average_points(points, 8, &x0);
  test_poly_shape_function_gradients(1, &x0, points, 8, coeffs, false);
}

void test_weighted_p1_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0, points[8];
  generate_random_points(8, points);
  average_points(points, 8, &x0);
  test_poly_shape_function_gradients(1, &x0, points, 8, coeffs, true);
}

void test_p2_fit(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[16];
  generate_random_points(16, points);
  test_poly_fit(2, NULL, points, 16, coeffs, false);
  average_points(points, 16, &x0);
  test_poly_fit(2, &x0, points, 16, coeffs, false);
}

void test_weighted_p2_fit(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[16];
  generate_random_points(16, points);
  test_poly_fit(2, NULL, points, 16, coeffs, true);
  average_points(points, 16, &x0);
  test_poly_fit(2, &x0, points, 16, coeffs, true);
}

void test_p2_shape_funcs(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[16];
  generate_random_points(16, points);
  average_points(points, 16, &x0);
  test_poly_shape_functions(2, &x0, points, 16, coeffs, false);
}

void test_weighted_p2_shape_funcs(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[16];
  generate_random_points(16, points);
  average_points(points, 16, &x0);
  test_poly_shape_functions(2, &x0, points, 16, coeffs, true);
}

void test_p2_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[16];
  generate_random_points(16, points);
  average_points(points, 16, &x0);
  test_poly_shape_function_gradients(2, &x0, points, 16, coeffs, false);
}

void test_weighted_p2_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[16];
  generate_random_points(16, points);
  average_points(points, 16, &x0);
  test_poly_shape_function_gradients(2, &x0, points, 16, coeffs, true);
}

void test_p3_fit(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  test_poly_fit(3, NULL, points, 30, coeffs, false);
  average_points(points, 30, &x0);
  test_poly_fit(3, &x0, points, 30, coeffs, false);
}

void test_weighted_p3_fit(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  test_poly_fit(3, NULL, points, 30, coeffs, true);
  average_points(points, 30, &x0);
  test_poly_fit(3, &x0, points, 30, coeffs, true);
}

void test_p3_shape_funcs(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  average_points(points, 30, &x0);
  test_poly_shape_functions(3, &x0, points, 30, coeffs, false);
}

void test_weighted_p3_shape_funcs(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  average_points(points, 30, &x0);
  test_poly_shape_functions(3, &x0, points, 30, coeffs, true);
}

void test_p3_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  average_points(points, 30, &x0);
  test_poly_shape_function_gradients(3, &x0, points, 30, coeffs, false);
}

void test_weighted_p3_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  average_points(points, 30, &x0);
  test_poly_shape_function_gradients(3, &x0, points, 30, coeffs, true);
}

int main(int argc, char* argv[]) 
{
  arbi_init(argc, argv);

  // Initialize the random number generator.
  unsigned seed = 1;
  initstate(seed, rng_state, RNG_SIZE);

  const UnitTest tests[] = 
  {
    unit_test(test_p0_fit),
    unit_test(test_p1_fit),
    unit_test(test_p2_fit), 
    unit_test(test_p3_fit),
    unit_test(test_weighted_p0_fit),
    unit_test(test_weighted_p1_fit),
    unit_test(test_weighted_p2_fit), 
    unit_test(test_weighted_p3_fit),
    unit_test(test_p0_shape_funcs),
    unit_test(test_p1_shape_funcs),
    unit_test(test_p2_shape_funcs),
    unit_test(test_p3_shape_funcs),
    unit_test(test_weighted_p0_shape_funcs),
    unit_test(test_weighted_p1_shape_funcs),
    unit_test(test_weighted_p2_shape_funcs),
    unit_test(test_weighted_p3_shape_funcs),
    unit_test(test_p0_shape_func_gradients),
    unit_test(test_p1_shape_func_gradients),
    unit_test(test_p2_shape_func_gradients),
    unit_test(test_p3_shape_func_gradients),
    unit_test(test_weighted_p0_shape_func_gradients),
    unit_test(test_weighted_p1_shape_func_gradients),
    unit_test(test_weighted_p2_shape_func_gradients),
    unit_test(test_weighted_p3_shape_func_gradients)
  };
  return run_tests(tests);
}
