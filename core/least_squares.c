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
    gradients[i].x = (x == 0) ? 0.0 : x*pow(point->x, x-1)*pow(point->y, y)*pow(point->z, z);
    gradients[i].y = (y == 0) ? 0.0 : pow(point->x, x)*y*pow(point->y, y-1)*pow(point->z, z);
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
    W(NULL, &points[n], x0, &Wd, &gradWd);
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
        moment_matrix[size*j+i] += Wd*basis[i]*basis[j];
      rhs[i] += Wd*basis[i]*data[n];
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
void dgemv(char *trans, int *m, int *n, double *alpha,
           void *a, int *lda, void *x, int *incx,
           double *beta, void *y, int *incy);

// Matrix-matrix multiplication: C := alpha*op(A)*op(B) + beta*C.
void dgemm(char *transa, char* transB, int *m, int *n, int *k, double *alpha,
           double *A, int *lda, double *B, int *ldb, double *beta, double *C, 
           int *ldc);

// Shape function basis.
struct poly_ls_shape_t 
{
  int p; // Order of basis.
  bool compute_gradients; // Compute gradients, or no?
  int dim; // Dimension of basis.
  int num_points; // Number of points in domain.
  point_t* points; // Points in domain.
  point_t x0; // Origin.
  double *A, *dAdx, *dAdy, *dAdz; // moment matrix and derivatives.
  double *AinvB; // Ainv * B.
  double *dBdx, *dBdy, *dBdz; // Derivatives of B.
  double *dAinvBdx, *dAinvBdy, *dAinvBdz; // Derivatives of Ainv*B.
  ls_weighting_func_t weighting_func; // Weighting function.
  void* w_context; // Context pointer for weighting function.
  void (*w_dtor)(void*); // Destructor for weighting function context pointer.
  double* weights; // Weight function values.
  vector_t* gradients; // Weight function gradients.
};

static void poly_ls_shape_free(void* context, void* dummy)
{
  poly_ls_shape_t* N = (poly_ls_shape_t*)context;
  if ((N->w_context != NULL) && (N->w_dtor != NULL))
    (*N->w_dtor)(N->w_context);
  if (N->points != NULL)
    free(N->points);
  free(N->A);
  free(N->dAdx);
  free(N->dAdy);
  free(N->dAdz);
  if (N->AinvB != NULL)
    free(N->AinvB);
  if (N->dBdx != NULL)
  {
    free(N->dBdx);
    free(N->dBdy);
    free(N->dBdz);
  }
  if (N->weights != NULL)
    free(N->weights);
  if (N->gradients != NULL)
    free(N->gradients);
  free(N);
}

static void no_weighting_func(void* context, point_t* x, point_t* x0, double* W, vector_t* gradient)
{
  *W = 1.0;
  gradient->x = gradient->y = gradient->z = 0.0;
//  arbi_error("No weighting function has been set for this LS shape function.");
}

poly_ls_shape_t* poly_ls_shape_new(int p, bool compute_gradients)
{
  ASSERT(p >= 0);
  ASSERT(p < 4);
  poly_ls_shape_t* N = GC_MALLOC(sizeof(poly_ls_shape_t));
  N->p = p;
  N->compute_gradients = compute_gradients;
  N->weighting_func = &no_weighting_func;
  N->w_context = NULL;
  N->w_dtor = NULL;
  N->dim = poly_ls_basis_size(p);
  N->num_points = 0;
  N->points = NULL;
  N->A = malloc(sizeof(double)*N->dim*N->dim);
  N->dAdx = malloc(sizeof(double)*N->dim*N->dim);
  N->dAdy = malloc(sizeof(double)*N->dim*N->dim);
  N->dAdz = malloc(sizeof(double)*N->dim*N->dim);
  N->AinvB = NULL;
  N->dBdx = NULL;
  N->dBdy = NULL;
  N->dBdz = NULL;
  N->dAinvBdx = N->dAinvBdy = N->dAinvBdz = NULL;
  N->weights = NULL;
  N->gradients = NULL;
  GC_register_finalizer(N, &poly_ls_shape_free, N, NULL, NULL);
  return N;
}

void poly_ls_shape_set_domain(poly_ls_shape_t* N, point_t* x0, point_t* points, int num_points)
{
  int dim = N->dim;
  if (num_points != N->num_points)
  {
    N->num_points = num_points;
    N->points = realloc(N->points, sizeof(point_t)*num_points);
    memcpy(N->points, points, sizeof(point_t)*num_points);
    N->AinvB = realloc(N->AinvB, sizeof(double)*dim*num_points);
    N->dBdx = realloc(N->dBdx, sizeof(double)*dim*num_points);
    N->dBdy = realloc(N->dBdy, sizeof(double)*dim*num_points);
    N->dBdz = realloc(N->dBdz, sizeof(double)*dim*num_points);
    N->dAinvBdx = realloc(N->dAinvBdx, sizeof(double)*dim*num_points);
    N->dAinvBdy = realloc(N->dAinvBdy, sizeof(double)*dim*num_points);
    N->dAinvBdz = realloc(N->dAinvBdz, sizeof(double)*dim*num_points);
    N->weights = realloc(N->weights, sizeof(double)*num_points);
    N->gradients = realloc(N->gradients, sizeof(vector_t)*num_points);
  }
  N->x0.x = x0->x;
  N->x0.y = x0->y;
  N->x0.z = x0->z;

  // Compute the moment matrix A and the basis matrix B.
  memset(N->A, 0, sizeof(double)*dim*dim);
  memset(N->dAdx, 0, sizeof(double)*dim*dim);
  memset(N->dAdy, 0, sizeof(double)*dim*dim);
  memset(N->dAdz, 0, sizeof(double)*dim*dim);
  memset(N->dBdx, 0, sizeof(double)*dim*num_points);
  memset(N->dBdy, 0, sizeof(double)*dim*num_points);
  memset(N->dBdz, 0, sizeof(double)*dim*num_points);
  memset(N->AinvB, 0, sizeof(double)*dim*num_points);
  memset(N->dAinvBdx, 0, sizeof(double)*dim*num_points);
  memset(N->dAinvBdy, 0, sizeof(double)*dim*num_points);
  memset(N->dAinvBdz, 0, sizeof(double)*dim*num_points);
  double basis[dim];
  for (int n = 0; n < num_points; ++n)
  {
    // Expand about x0.
    point_t y = {.x = points[n].x - x0->x, 
                 .y = points[n].y - x0->y,
                 .z = points[n].z - x0->z};
    N->weighting_func(N->w_context, &points[n], x0, &N->weights[n], &N->gradients[n]);
    compute_poly_ls_basis_vector(N->p, &y, basis);
    for (int i = 0; i < dim; ++i)
    {
      N->AinvB[dim*n+i] = N->weights[n]*basis[i];
      N->dBdx[dim*n+i] = N->gradients[n].x*basis[i];
      N->dBdy[dim*n+i] = N->gradients[n].y*basis[i];
      N->dBdz[dim*n+i] = N->gradients[n].z*basis[i];
      for (int j = 0; j < dim; ++j)
      {
        N->A[dim*j+i] += basis[i]*N->weights[n]*basis[j];
        N->dAdx[dim*j+i] += basis[i]*N->gradients[n].x*basis[j];
        N->dAdy[dim*j+i] += basis[i]*N->gradients[n].y*basis[j];
        N->dAdz[dim*j+i] += basis[i]*N->gradients[n].z*basis[j];
      }
    }
  }

  // Factor the moment matrix.
  int pivot[dim], info;
  dgetrf(&dim, &dim, N->A, &dim, pivot, &info);
  ASSERT(info == 0);

  // Compute Ainv * B.
  char no_trans = 'N';
  dgetrs(&no_trans, &dim, &num_points, N->A, &dim, pivot, N->AinvB, &dim, &info);
  ASSERT(info == 0);

  // If we are in the business of computing gradients, compute the 
  // partial derivatives of Ainv * B.
  if (N->compute_gradients)
  {
    // The partial derivatives of A inverse are:
    // d(Ainv) = -Ainv * dA * Ainv, so
    // d(Ainv * B) = -Ainv * dA * Ainv * B + Ainv * dB
    //             = Ainv * (-dA * Ainv * B + dB).

    // We left-multiply Ainv*B by the gradient of A, placing the results 
    // in dAinvBdx, dAinvBdy, and dAinvBdz.
    double alpha = 1.0, beta = 0.0;
    dgemm(&no_trans, &no_trans, &N->dim, &N->num_points, &N->dim, &alpha, 
          N->dAdx, &N->dim, N->AinvB, &N->dim, &beta, N->dAinvBdx, &N->dim);
    dgemm(&no_trans, &no_trans, &N->dim, &N->num_points, &N->dim, &alpha, 
          N->dAdy, &N->dim, N->AinvB, &N->dim, &beta, N->dAinvBdy, &N->dim);
    dgemm(&no_trans, &no_trans, &N->dim, &N->num_points, &N->dim, &alpha, 
          N->dAdz, &N->dim, N->AinvB, &N->dim, &beta, N->dAinvBdz, &N->dim);

    // Flip the sign of dA * Ainv * B, and add dB.
    for (int i = 0; i < dim*num_points; ++i)
    {
      N->dAinvBdx[i] = -N->dAinvBdx[i] + N->dBdx[i];
      N->dAinvBdy[i] = -N->dAinvBdy[i] + N->dBdy[i];
      N->dAinvBdz[i] = -N->dAinvBdz[i] + N->dBdz[i];
    }

    // Now "left-multiply by Ainv" by solving the equation (e.g.)
    // A * (dAinvBdx) = (-dA * Ainv * B + dB).
    dgetrs(&no_trans, &dim, &num_points, N->A, &dim, pivot, N->dAinvBdx, &dim, &info);
    ASSERT(info == 0);
    dgetrs(&no_trans, &dim, &num_points, N->A, &dim, pivot, N->dAinvBdy, &dim, &info);
    ASSERT(info == 0);
    dgetrs(&no_trans, &dim, &num_points, N->A, &dim, pivot, N->dAinvBdz, &dim, &info);
    ASSERT(info == 0);
  }
}

void poly_ls_shape_compute(poly_ls_shape_t* N, point_t* x, double* values)
{
  ASSERT(N->AinvB != NULL);
  double basis[N->dim];

  // values^T = basis^T * Ainv * B (or values = (Ainv * B)^T * basis.)
  double alpha = 1.0, beta = 0.0;
  int one = 1;
  char trans = 'T';
  point_t y = {.x = x->x - N->x0.x, 
               .y = x->y - N->x0.y,
               .z = x->z - N->x0.z};
  compute_poly_ls_basis_vector(N->p, &y, basis);
  dgemv(&trans, &N->dim, &N->num_points, &alpha, N->AinvB, &N->dim, basis, &one, &beta, values, &one);
}

void poly_ls_shape_compute_gradients(poly_ls_shape_t* N, point_t* x, double* values, vector_t* gradients)
{
  ASSERT(N->compute_gradients);
  ASSERT(N->AinvB != NULL);

  int dim = N->dim;

  // values^T = basis^T * Ainv * B (or values = (Ainv * B)^T * basis.)
  double alpha = 1.0, beta = 0.0;
  int one = 1;
  char trans = 'T';
  point_t y = {.x = x->x - N->x0.x, 
               .y = x->y - N->x0.y,
               .z = x->z - N->x0.z};
  double basis[dim];
  compute_poly_ls_basis_vector(N->p, &y, basis);
  dgemv(&trans, &N->dim, &N->num_points, &alpha, N->AinvB, &N->dim, basis, &one, &beta, values, &one);

  // Now compute the gradients.

  // 1st term: gradient of basis, dotted with Ainv * B.
  vector_t basis_grads[dim];
  compute_poly_ls_basis_gradient(N->p, &y, basis_grads);
  double dpdx[dim], dpdy[dim], dpdz[dim];
  for (int i = 0; i < dim; ++i)
  {
    dpdx[i] = basis_grads[i].x;
    dpdy[i] = basis_grads[i].y;
    dpdz[i] = basis_grads[i].z;
  }
  double dpdx_AinvB[N->num_points], dpdy_AinvB[N->num_points], dpdz_AinvB[N->num_points];
  dgemv(&trans, &N->dim, &N->num_points, &alpha, N->AinvB, &N->dim, dpdx, &one, &beta, dpdx_AinvB, &one);
  dgemv(&trans, &N->dim, &N->num_points, &alpha, N->AinvB, &N->dim, dpdy, &one, &beta, dpdy_AinvB, &one);
  dgemv(&trans, &N->dim, &N->num_points, &alpha, N->AinvB, &N->dim, dpdz, &one, &beta, dpdz_AinvB, &one);

  // Second term: basis dotted with gradient of Ainv * B.
  double p_dAinvBdx[N->num_points], p_dAinvBdy[N->num_points], p_dAinvBdz[N->num_points];
  dgemv(&trans, &N->dim, &N->num_points, &alpha, N->dAinvBdx, &N->dim, basis, &one, &beta, p_dAinvBdx, &one);
  dgemv(&trans, &N->dim, &N->num_points, &alpha, N->dAinvBdy, &N->dim, basis, &one, &beta, p_dAinvBdy, &one);
  dgemv(&trans, &N->dim, &N->num_points, &alpha, N->dAinvBdz, &N->dim, basis, &one, &beta, p_dAinvBdz, &one);

  // Gradients are the sum of these terms.
  for (int i = 0; i < N->num_points; ++i)
  {
    gradients[i].x = dpdx_AinvB[i] + p_dAinvBdx[i];
    gradients[i].y = dpdy_AinvB[i] + p_dAinvBdy[i];
    gradients[i].z = dpdz_AinvB[i] + p_dAinvBdz[i];
  }
}

void poly_ls_shape_compute_constraint_transform(poly_ls_shape_t* N, int* constraint_indices, int num_constraints,
                                                double* a, double* b, double* c, double* d, double* e,
                                                double* A, double* B)
{
  ASSERT(N->p > 0); // Constraints cannot be applied to constants.
  ASSERT(N->compute_gradients);
  ASSERT(constraint_indices != NULL);
  ASSERT(num_constraints < N->num_points);
  ASSERT(a != NULL);
  ASSERT(b != NULL);
  ASSERT(c != NULL);
  ASSERT(d != NULL);
  ASSERT(e != NULL);
  ASSERT(A != NULL);
  ASSERT(B != NULL);

  // Set up the constraint matrices.
  double amat[num_constraints*num_constraints];
  for (int i = 0; i < num_constraints; ++i)
  {
    // Compute the shape functions at xi.
    double N_vals[N->num_points];
    vector_t N_grads[N->num_points];
    poly_ls_shape_compute_gradients(N, &N->points[constraint_indices[i]], N_vals, N_grads);

    // Now set up the left and right hand sides of the equation for the constraint.
    int constraint = 0;
    for (int j = 0; j < N->num_points; ++j)
    {
      bool constrained = false;
      for (int cc = 0; cc < num_constraints; ++cc)
      {
        if (j == constraint_indices[cc]) 
        {
          constrained = true;
          break;
        }
      }
      if (constrained) 
      {
        amat[num_constraints*constraint+i] = a[i]*N_vals[j] + b[i]*N_grads[j].x + c[i]*N_grads[j].y + d[i]*N_grads[j].z;
        A[num_constraints*j+i] = 0.0;
        constraint++;
      }
      else
      {
        A[num_constraints*j+i] = -a[i]*N_vals[j] - b[i]*N_grads[j].x - c[i]*N_grads[j].y - d[i]*N_grads[j].z;
      }
    }
  }

  // Compute A by solving the linear system.
  int pivot[num_constraints], info;
  dgetrf(&num_constraints, &num_constraints, amat, &num_constraints, pivot, &info);
  ASSERT(info == 0);
  char no_trans = 'N';
  dgetrs(&no_trans, &num_constraints, &N->num_points, amat, &num_constraints, pivot, A, &num_constraints, &info);
  ASSERT(info == 0);

  // Compute B = amatinv * e.
  int one = 1;
  memcpy(B, e, sizeof(double)*num_constraints);
  dgetrs(&no_trans, &num_constraints, &one, amat, &num_constraints, pivot, B, &num_constraints, &info);
  ASSERT(info == 0);
}

typedef struct 
{
  int A;
  double B;
} simple_weighting_func_params_t;

static void simple_weighting_func(void* context, point_t* x, point_t* x0, double* W, vector_t* gradient)
{
  simple_weighting_func_params_t* params = (simple_weighting_func_params_t*)context;
  double D = point_distance(x, x0);
  *W = 1.0 / (pow(D, params->A) + pow(params->B, params->A));
  if (D == 0.0)
  {
    gradient->x = gradient->y = gradient->z = 0.0;
  }
  else
  {
    double dDdx = x->x / D, dDdy = x->y / D, dDdz = x->z / D;
    double deriv_term = -(*W)*(*W) * params->A * pow(D, params->A-1);
    gradient->x = deriv_term * dDdx;
    gradient->y = deriv_term * dDdy;
    gradient->z = deriv_term * dDdz;
  }
}

void poly_ls_shape_set_simple_weighting_func(poly_ls_shape_t* N, int A, double B)
{
  ASSERT(A > 0);
  ASSERT(B > 0);
  N->weighting_func = &simple_weighting_func;
  simple_weighting_func_params_t* params = malloc(sizeof(simple_weighting_func_params_t));
  params->A = A;
  params->B = B;
  N->w_context = params;
  N->w_dtor = &free;
}

#ifdef __cplusplus
}
#endif

