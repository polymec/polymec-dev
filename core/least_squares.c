#include <stdlib.h>
#include <string.h>
#include <gc/gc.h>
#include "core/least_squares.h"

#ifdef __cplusplus
extern "C" {
#endif

struct multi_index_t 
{
  int p, x_order, y_order, z_order;
  int offset; // Offset in flattened index space.
};

static void multi_index_free(void* ctx, void* dummy)
{
  multi_index_t* m = (multi_index_t*)ctx;
  free(m);
}

multi_index_t* multi_index_new(int p)
{
  ASSERT(p >= 0);
  ASSERT(p < 4);
  multi_index_t* m = GC_MALLOC(sizeof(multi_index_t));
  m->p = p;
  m->x_order = m->y_order = m->z_order = 0;
  m->offset = 0;
  GC_register_finalizer(m, &multi_index_free, m, NULL, NULL);
  return m;
}

static int multi_index_sizes[] = {1, 4, 10, 18};

bool multi_index_next(multi_index_t* m, int* x_order, int* y_order, int* z_order)
{
  ASSERT(m->p >= 0);
  ASSERT(m->p < 4);
  if (m->offset == -1)
    return false;
  if (m->p == 0)
  {
    *x_order = m->x_order;
    *y_order = m->y_order;
    *z_order = m->z_order;
    m->offset = -1;
  }
  if (m->p == 1)
  {
    static const multi_index_t multi_index[] = 
    {
      {.p = 1, .x_order = 0, .y_order = 0, .z_order = 0, .offset = 0}, 
      {.p = 1, .x_order = 1, .y_order = 0, .z_order = 0, .offset = 1}, 
      {.p = 1, .x_order = 0, .y_order = 1, .z_order = 0, .offset = 2}, 
      {.p = 1, .x_order = 0, .y_order = 0, .z_order = 1, .offset = 3},
      {.p = 1, .offset = -1}
    };
    *x_order = m->x_order;
    *y_order = m->y_order;
    *z_order = m->z_order;
    *m = multi_index[m->offset+1];
  }
  else if (m->p == 2)
  {
    static const multi_index_t multi_index[] = 
    {
      {.p = 2, .x_order = 0, .y_order = 0, .z_order = 0, .offset = 0}, 
      {.p = 2, .x_order = 1, .y_order = 0, .z_order = 0, .offset = 1}, 
      {.p = 2, .x_order = 0, .y_order = 1, .z_order = 0, .offset = 2}, 
      {.p = 2, .x_order = 0, .y_order = 0, .z_order = 1, .offset = 3},
      {.p = 2, .x_order = 2, .y_order = 0, .z_order = 0, .offset = 4},
      {.p = 2, .x_order = 1, .y_order = 1, .z_order = 0, .offset = 5},
      {.p = 2, .x_order = 1, .y_order = 0, .z_order = 1, .offset = 6},
      {.p = 2, .x_order = 0, .y_order = 2, .z_order = 0, .offset = 7},
      {.p = 2, .x_order = 0, .y_order = 1, .z_order = 1, .offset = 8},
      {.p = 2, .x_order = 0, .y_order = 0, .z_order = 2, .offset = 9},
      {.p = 1, .offset = -1}
    };
    *x_order = m->x_order;
    *y_order = m->y_order;
    *z_order = m->z_order;
    *m = multi_index[m->offset+1];
  }
  else if (m->p == 3)
  {
    static const multi_index_t multi_index[] = 
    {
      {.p = 3, .x_order = 0, .y_order = 0, .z_order = 0, .offset = 0}, 
      {.p = 3, .x_order = 1, .y_order = 0, .z_order = 0, .offset = 1}, 
      {.p = 3, .x_order = 0, .y_order = 1, .z_order = 0, .offset = 2}, 
      {.p = 3, .x_order = 0, .y_order = 0, .z_order = 1, .offset = 3},
      {.p = 3, .x_order = 2, .y_order = 0, .z_order = 0, .offset = 4},
      {.p = 3, .x_order = 1, .y_order = 1, .z_order = 0, .offset = 5},
      {.p = 3, .x_order = 1, .y_order = 0, .z_order = 1, .offset = 6},
      {.p = 3, .x_order = 0, .y_order = 2, .z_order = 0, .offset = 7},
      {.p = 3, .x_order = 0, .y_order = 1, .z_order = 1, .offset = 8},
      {.p = 3, .x_order = 0, .y_order = 0, .z_order = 2, .offset = 9},
      {.p = 3, .x_order = 3, .y_order = 0, .z_order = 0, .offset = 10},
      {.p = 3, .x_order = 2, .y_order = 1, .z_order = 0, .offset = 11},
      {.p = 3, .x_order = 2, .y_order = 0, .z_order = 1, .offset = 12},
      {.p = 3, .x_order = 1, .y_order = 2, .z_order = 0, .offset = 13},
      {.p = 3, .x_order = 1, .y_order = 1, .z_order = 1, .offset = 14},
      {.p = 3, .x_order = 0, .y_order = 3, .z_order = 0, .offset = 15},
      {.p = 3, .x_order = 0, .y_order = 2, .z_order = 1, .offset = 16},
      {.p = 3, .x_order = 0, .y_order = 1, .z_order = 2, .offset = 17},
      {.p = 1, .offset = -1}
    };
    *x_order = m->x_order;
    *y_order = m->y_order;
    *z_order = m->z_order;
    *m = multi_index[m->offset+1];
  }
  return true; 
}

void multi_index_reset(multi_index_t* m)
{
  m->offset = 0;
  m->x_order = m->y_order = m->z_order = 0;
}

int multi_index_order(multi_index_t* m)
{
  return m->p;
}

int multi_index_size(multi_index_t* m)
{
  return multi_index_sizes[m->p];
}

int poly_ls_basis_size(int p)
{
  return multi_index_sizes[p];
}

void compute_poly_ls_basis_vector(int p, point_t* point, double* basis)
{
  multi_index_t* m = multi_index_new(p);
  int i = 0, x, y, z;
  while (multi_index_next(m, &x, &y, &z))
    basis[i++] = pow(point->x, x)*pow(point->y, y)*pow(point->z, z);
  m = NULL;
}

void compute_poly_ls_basis_gradient(int p, point_t* point, vector_t* gradients)
{
  multi_index_t* m = multi_index_new(p);
  int i = 0, x, y, z;
  while (multi_index_next(m, &x, &y, &z))
  {
    gradients[i++].x = (x == 0) ? 0.0 : x*pow(point->x, x-1)*pow(point->y, y)*pow(point->z, z);
    gradients[i++].y = (y == 0) ? 0.0 : pow(point->x, x)*y*pow(point->y, y-1)*pow(point->z, z);
    gradients[i++].z = (z == 0) ? 0.0 : pow(point->x, x)*pow(point->y, y)*z*pow(point->z, z-1);
  }
  m = NULL;
}

void compute_poly_ls_system(int p, point_t* x0, point_t* points, int num_points, 
                            double* data, double* moment_matrix, double* rhs)
{
  ASSERT(p >= 0);
  ASSERT(p < 4);
  ASSERT(num_points >= multi_index_sizes[p]);
  ASSERT(moment_matrix != NULL);
  ASSERT(rhs != NULL);
  int size = multi_index_sizes[p];
  double basis[size];

  // Zero the system.
  memset(moment_matrix, 0, sizeof(double)*size*size);
  memset(rhs, 0, sizeof(double)*size);
 
  for (int n = 0; n < num_points; ++n)
  {
    if (x0 != NULL)
    {
      point_t y = {.x = points[n].x - x0->x, 
                   .y = points[n].y - x0->y,
                   .z = points[n].z - x0->z};
      compute_poly_ls_basis_vector(p, &y, basis);
    }
    else
    {
      compute_poly_ls_basis_vector(p, &points[n], basis);
    }
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
        moment_matrix[size*j+i] += basis[i]*basis[j];
      rhs[i] += basis[i]*data[n];
    }
  }
}

void compute_weighted_poly_ls_system(int p, ls_weighting_func_t W, point_t* x0, point_t* points, int num_points, 
                                     double* data, double* moment_matrix, double* rhs)
{
  ASSERT(p >= 0);
  ASSERT(p < 4);
  ASSERT(moment_matrix != NULL);
  ASSERT(rhs != NULL);
  int size = multi_index_sizes[p];
  double basis[size];

  memset(moment_matrix, 0, sizeof(double)*size*size);
  memset(rhs, 0, sizeof(double)*size);
 
  for (int n = 0; n < num_points; ++n)
  {
    double d;
    if (x0 != NULL)
    {
      point_t y = {.x = points[n].x - x0->x, 
                   .y = points[n].y - x0->y,
                   .z = points[n].z - x0->z};
      compute_poly_ls_basis_vector(p, &y, basis);
      d = y.x*y.x + y.y*y.y + y.z*y.z;
    }
    else
    {
      compute_poly_ls_basis_vector(p, &points[n], basis);
      d = points[n].x*points[n].x + points[n].y*points[n].y + points[n].z*points[n].z;
    }
    double Wd;
    vector_t gradWd;
    W(NULL, d, &Wd, &gradWd);
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
        moment_matrix[size*j+i] += Wd*basis[i]*basis[j];
      rhs[i] += basis[i]*data[n];
    }
  }
}

// Some LAPACK prototypes.

// LU factorization.
void dgetrf(int *N, int *NRHS, double *A, int *LDA, int *IPIV, int *INFO);

// Solution of a linear system using an LU factorization.
void dgetrs(char *TRANS, int *N, int *NRHS, double *A, 
            int *LDA, int *IPIV, double *B, int *LDB, int *INFO);

// Matrix-vector multiplication: y := alpha*A*x + beta*y.
void dgemv(const char *trans, int *m, int *n, double *alpha,
           void *a, int *lda, void *x, int *incx,
           double *beta, void *y, int *incy);

// Shape function basis.
struct poly_ls_shape_basis_t 
{
  int p; // Order of basis.
  ls_weighting_func_t weighting_func;
};

static void poly_ls_shape_basis_free(void* context, void* dummy)
{
  poly_ls_shape_basis_t* N = (poly_ls_shape_basis_t*)context;
  free(N);
}

poly_ls_shape_basis_t* poly_ls_shape_basis_new(int p, ls_weighting_func_t weighting_func)
{
  ASSERT(p >= 0);
  ASSERT(p < 4);
  poly_ls_shape_basis_t* N = GC_MALLOC(sizeof(poly_ls_shape_basis_t));
  N->p = p;
  N->weighting_func = weighting_func;
  GC_register_finalizer(N, &poly_ls_shape_basis_free, N, NULL, NULL);
  return N;
}

void poly_ls_shape_basis_compute(poly_ls_shape_basis_t* N, point_t* x0, point_t* points, int num_points, point_t* x, double* values)
{
  int dim = poly_ls_basis_size(N->p);

  // Compute the moment matrix A and the basis matrix B.
  double A[dim*dim], basis[dim], B[dim*num_points], W[num_points];
  memset(A, 0, sizeof(double)*dim*dim);
  memset(B, 0, sizeof(double)*dim*num_points);
  for (int n = 0; n < num_points; ++n)
  {
    // Expand about x0.
    point_t y = {.x = points[n].x - x0->x, 
                 .y = points[n].y - x0->y,
                 .z = points[n].z - x0->z};
    double d = y.x*y.x + y.y*y.y + y.z*y.z;
    W[n] = 1.0;
    vector_t gradWn = {.x = 0, .y = 0, .z = 0};
    if (N->weighting_func != NULL)
      N->weighting_func(NULL, d, &W[n], &gradWn);
    compute_poly_ls_basis_vector(N->p, &y, basis);
    for (int i = 0; i < dim; ++i)
    {
      B[dim*n+i] = W[n]*basis[i];
      for (int j = 0; j < dim; ++j)
        A[dim*j+i] += basis[i]*W[n]*basis[j];
    }
  }
printf("P = [");
for (int i = 0; i < num_points*dim; ++i)
  printf("%g ", B[i]);
printf("]\n");
printf("A = [");
for (int i = 0; i < dim*dim; ++i)
  printf("%g ", A[i]);
printf("]\n");

  // Factor the moment matrix.
  int pivot[dim], info;
  dgetrf(&dim, &dim, A, &dim, pivot, &info);
  ASSERT(info == 0);

  // Now compute the values of the shape function basis at x.

  // Ainv * B -> B.
  char no_trans = 'N';
  dgetrs(&no_trans, &dim, &num_points, A, &dim, pivot, B, &dim, &info);
  ASSERT(info == 0);
printf("B = [");
for (int i = 0; i < num_points*dim; ++i)
  printf("%g ", B[i]);
printf("]\n");

  // values^T = basis^T * Ainv * B (or values = (Ainv * B)^T * basis.)
  {
    double alpha = 1.0, beta = 0.0;
    int one = 1;
    char trans = 'T';
    point_t y = {.x = x->x - x0->x, 
                 .y = x->y - x0->y,
                 .z = x->z - x0->z};
    compute_poly_ls_basis_vector(N->p, &y, basis);
    dgemv(&trans, &dim, &num_points, &alpha, B, &dim, basis, &one, &beta, values, &one);
  }
}

void poly_ls_shape_basis_compute_gradients(poly_ls_shape_basis_t* N, point_t* x0, point_t* points, int num_points, point_t* x, double* values, vector_t* gradients)
{
  int dim = poly_ls_basis_size(N->p);

  // Compute the moment matrix and its gradient.
  double A[dim*dim], basis[dim], P[dim*dim], W[num_points],
         gradAx[dim*dim], gradAy[dim*dim], gradAz[dim*dim];
  vector_t gradW[num_points];
  memset(A, 0, sizeof(double)*dim*dim);
  memset(gradAx, 0, sizeof(double)*dim*dim);
  memset(gradAy, 0, sizeof(double)*dim*dim);
  memset(gradAz, 0, sizeof(double)*dim*dim);
  memset(gradW, 0, sizeof(vector_t)*num_points);
  for (int n = 0; n < num_points; ++n)
  {
    // Expand about x0.
    point_t y = {.x = points[n].x - x0->x, 
                 .y = points[n].y - x0->y,
                 .z = points[n].z - x0->z};

    // Compute the weighting functions.
    double d = y.x*y.x + y.y*y.y + y.z*y.z;
    W[n] = 1.0;
    if (N->weighting_func != NULL)
      N->weighting_func(NULL, d, &W[n], &gradW[n]);

    compute_poly_ls_basis_vector(N->p, &y, basis);
    memcpy(&P[dim*n], basis, dim*sizeof(double));
    for (int i = 0; i < dim; ++i)
    {
      for (int j = 0; j < dim; ++j)
      {
        A[dim*j+i] += basis[i]*W[n]*basis[j];
        gradAx[dim*j+i] += basis[i]*gradW[n].x*basis[j];
        gradAy[dim*j+i] += basis[i]*gradW[n].y*basis[j];
        gradAz[dim*j+i] += basis[i]*gradW[n].z*basis[j];
      }
    }
  }

  // Factor the moment matrix.
  int lda = dim, pivot[dim], info;
  dgetrf(&dim, &dim, A, &lda, pivot, &info);
  ASSERT(info == 0);

  // Now compute the values and gradients of the shape function basis at x.
  char trans = 'N';
  int ldb = dim, one = 1;
  double alpha = 1.0, beta = 0.0;
  double Ainv_p[dim];
  memset(values, 0, sizeof(double)*num_points);
  memset(gradients, 0, sizeof(vector_t)*num_points);
  {
    compute_poly_ls_basis_vector(N->p, x, basis);
    compute_poly_ls_basis_gradient(N->p, x, gradients);
    for (int n = 0; n < num_points; ++n)
    {
      // Compute Ainv * p.
      memcpy(Ainv_p, &P[dim*n], dim*sizeof(double));
      dgetrs(&trans, &dim, &one, A, &lda, pivot, Ainv_p, &ldb, &info);

      // First pass: pT_Ainv_p and dpT_Ainv_p.
      double pT_Ainv_p = 0.0;
      vector_t gradpT_Ainv_p = {.x = 0., .y = 0., .z = 0.};
      for (int i = 0; i < dim; ++i)
      {
        pT_Ainv_p += basis[i]*Ainv_p[i];
        gradpT_Ainv_p.x += gradients[i].x * Ainv_p[i];
        gradpT_Ainv_p.y += gradients[i].y * Ainv_p[i];
        gradpT_Ainv_p.z += gradients[i].z * Ainv_p[i];
      }
      values[n] = pT_Ainv_p;

      // Second pass: Ainv_gradA_Ainv_p.
      double Ainv_gradA_Ainv_px[dim],
             Ainv_gradA_Ainv_py[dim],
             Ainv_gradA_Ainv_pz[dim];
      for (int i = 0; i < dim; ++i)
      {
        // gradA * Ainv_p -> Ainv_gradA_Ainv_p.
        dgemv(&trans, &dim, &dim, &alpha, gradAx, &dim, Ainv_p, &one, &beta,
              Ainv_gradA_Ainv_px, &one);
        dgemv(&trans, &dim, &dim, &alpha, gradAy, &dim, Ainv_p, &one, &beta,
              Ainv_gradA_Ainv_py, &one);
        dgemv(&trans, &dim, &dim, &alpha, gradAz, &dim, Ainv_p, &one, &beta,
              Ainv_gradA_Ainv_pz, &one);

        // Solve A * X = gradA * Ainv_p; X -> Ainv_gradA_Ainv_p.
        dgetrs(&trans, &dim, &one, A, &lda, pivot, Ainv_gradA_Ainv_px, &ldb, &info);
        ASSERT(info == 0);
        dgetrs(&trans, &dim, &one, A, &lda, pivot, Ainv_gradA_Ainv_py, &ldb, &info);
        ASSERT(info == 0);
        dgetrs(&trans, &dim, &one, A, &lda, pivot, Ainv_gradA_Ainv_pz, &ldb, &info);
        ASSERT(info == 0);
      }
      
      // Third pass: pT_Ainv_gradA_Ainv_p.
      vector_t pT_Ainv_gradA_Ainv_p = {.x = 0., .y = 0., .z = 0.};
      for (int i = 0; i < dim; ++i)
      {
        pT_Ainv_gradA_Ainv_p.x += basis[i] * Ainv_gradA_Ainv_px[i];
        pT_Ainv_gradA_Ainv_p.y += basis[i] * Ainv_gradA_Ainv_py[i];
        pT_Ainv_gradA_Ainv_p.y += basis[i] * Ainv_gradA_Ainv_pz[i];
      }

      // Third pass: compute the shape function gradients.
      for (int i = 0; i < dim; ++i)
      {
        gradients[n].x += (gradpT_Ainv_p.x - pT_Ainv_gradA_Ainv_p.x) * W[n] + 
                          pT_Ainv_p * gradW[n].x;
        gradients[n].y += (gradpT_Ainv_p.y - pT_Ainv_gradA_Ainv_p.y) * W[n] + 
                          pT_Ainv_p * gradW[n].y;
        gradients[n].z += (gradpT_Ainv_p.z - pT_Ainv_gradA_Ainv_p.z) * W[n] + 
                          pT_Ainv_p * gradW[n].z;
      }
    }
  }
}

#ifdef __cplusplus
}
#endif

