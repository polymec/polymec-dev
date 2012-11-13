#include <string.h>
#include <stdlib.h>
#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"
#include "core/unordered_map.h"
#include "core/least_squares.h"
#include "core/constant_st_func.h"
#include "geometry/create_cubic_lattice_mesh.h"
#include "geometry/create_cvt.h"
#include "geometry/plane.h"
#include "geometry/cylinder.h"
#include "geometry/sphere.h"
#include "geometry/intersection.h"
#include "poisson/poisson_model.h"
#include "poisson/laplacian_op.h"

#ifdef __cplusplus
extern "C" {
#endif

// Boundary condition structure.
// This represents a generic boundary condition: 
// alpha * phi + beta * dphi/dn = F
typedef struct
{
  double alpha, beta;
  st_func_t* F;
} poisson_bc_t;

static poisson_bc_t* create_bc(double alpha, double beta, st_func_t* F)
{
  ASSERT(F != NULL);
  poisson_bc_t* bc = malloc(sizeof(poisson_bc_t));
  bc->alpha = alpha;
  bc->beta = beta;
  bc->F = F;
  return bc;
}

static void free_bc(void* bc)
{
  poisson_bc_t* pbc = (poisson_bc_t*)bc;
  pbc->F = NULL;
}

// Metadata for boundary cells (those cells that touch boundary faces).
typedef struct
{
  int num_neighbor_cells;
  int* neighbor_cells;

  int num_boundary_faces;
  int* boundary_faces;
  poisson_bc_t** bc_for_face;
} boundary_cell_t;

// Poisson model context structure.
typedef struct 
{
  mesh_t* mesh;             // Mesh.
  st_func_t* rhs;           // Right-hand side function.
  double* phi;              // Solution array.
  lin_op_t* L;              // Laplacian operator.

  str_ptr_unordered_map_t* bcs; // Boundary conditions.
  poly_ls_shape_t* shape;       // Least-squares shape functions.

  KSP solver;               // Linear solver.
  Mat A;                    // Laplacian matrix.
  Vec x;                    // Solution vector.
  Vec b;                    // RHS vector. 

  int* boundary_cells;      // List of cells adjoining boundary faces.
  int num_boundary_cells;   // Number of boundary cells.
  int_ptr_unordered_map_t*  boundary_cell_info;

  bool initialized;         // Initialized flag.
  MPI_Comm comm;            // MPI communicator.

} poisson_t;

// A proper constructor.
static model_t* create_poisson(mesh_t* mesh, st_func_t* rhs, str_ptr_unordered_map_t* bcs)
{
  model_t* poisson = poisson_model_new(NULL);
  poisson_t* p = (poisson_t*)model_context(poisson);
  p->mesh = mesh;
  p->rhs = rhs;
  p->bcs = bcs;
  p->L = laplacian_op_new(p->mesh);
  return poisson;
}

//------------------------------------------------------------------------
//                            Benchmarks
//------------------------------------------------------------------------

// Mesh generation for benchmark problems.
static mesh_t* create_cube_mesh(int dim, int N)
{
  // Get the characteristic resolution.
  int N3[3] = {N, 1, 1};
  for (int d = 0; d < dim; ++d)
    N3[d] = N;

  return create_cubic_lattice_mesh(N3[0], N3[1], N3[2], 1);
}

static void run_analytic_problem(mesh_t* mesh, st_func_t* rhs, str_ptr_unordered_map_t* bcs, double t1, double t2, sp_func_t* solution, double* Lpnorms)
{
  // Create the model.
  model_t* model = create_poisson(mesh, rhs, bcs);

  // Run the thing.
  model_run(model, t1, t2);

  // Calculate the Lp norm of the error and write it to Lpnorms.

  // Clean up.
  model_free(model);
}

static void poisson_run_laplace_sov(int variant)
{
  // Dimension.
  int dim;
  // FIXME

  // RHS function.
  st_func_t* rhs;
  // FIXME

  // Boundary conditions.
  str_ptr_unordered_map_t* bcs = str_ptr_unordered_map_new();
  if ((variant % 3) == 0)
    str_ptr_unordered_map_insert(bcs, "boundary", create_bc(1.0, 0.0, NULL));
  else if ((variant % 3) == 1)
    str_ptr_unordered_map_insert(bcs, "boundary", create_bc(0.0, 1.0, NULL));
  else // if ((variant % 3) == 2)
    str_ptr_unordered_map_insert(bcs, "boundary", create_bc(1.0, 1.0, NULL));

  // Analytic solution.
  sp_func_t* sol;
  // FIXME = create_paraboloid_solution(variant);

  // Start/end times.
  double t1, t2;
  // FIXME

  // Base resolution.
  int N0;
  // FIXME

  // Number of refinements.
  int num_refinements;
  // FIXME

  // Do a convergence study.
  double Lpnorms[num_refinements][3];
  for (int iter = 0; iter < num_refinements; ++iter)
  {
    int N = pow(N0, iter+1);
    mesh_t* mesh = create_cube_mesh(dim, N);
    run_analytic_problem(mesh, rhs, bcs, t1, t2, sol, Lpnorms[iter]);
  }
}

static void poisson_run_paraboloid(int variant)
{
  // Dimension.
  int dim;
  // FIXME

  // RHS function.
  st_func_t* rhs;
  // FIXME = create_paraboloid_rhs(variant);

  // Boundary conditions.
  str_ptr_unordered_map_t* bcs = str_ptr_unordered_map_new();
  if ((variant % 3) == 0)
    str_ptr_unordered_map_insert(bcs, "boundary", create_bc(1.0, 0.0, NULL));
  else if ((variant % 3) == 1)
    str_ptr_unordered_map_insert(bcs, "boundary", create_bc(0.0, 1.0, NULL));
  else // if ((variant % 3) == 2)
    str_ptr_unordered_map_insert(bcs, "boundary", create_bc(1.0, 1.0, NULL));

  // Analytic solution.
  sp_func_t* sol;
  // FIXME = create_paraboloid_solution(variant);

  // Start/end times.
  double t1, t2;
  // FIXME

  // Base resolution.
  int N0;
  // FIXME

  // Number of refinements.
  int num_refinements;
  // FIXME

  // Do a convergence study.
  double Lpnorms[num_refinements][3];
  for (int iter = 0; iter < num_refinements; ++iter)
  {
    int N = pow(N0, iter+1);
    mesh_t* mesh = create_cube_mesh(dim, N);
    run_analytic_problem(mesh, rhs, bcs, t1, t2, sol, Lpnorms[iter]);
  }
}

//------------------------------------------------------------------------
//                        Model implementation
//------------------------------------------------------------------------

// Apply boundary conditions to a set of boundary cells.
static void apply_bcs(int* boundary_cells,
                      int num_boundary_cells,
                      int_ptr_unordered_map_t* boundary_cell_info,
                      mesh_t* mesh,
                      poly_ls_shape_t* shape,
                      double t,
                      Mat A,
                      Vec b)
{
  // Zero the rows corresponding to boundary cells.
  {
    MatZeroRows(A, num_boundary_cells, boundary_cells, 0.0, NULL, NULL);
    VecAssemblyBegin(b);
    double zeros[num_boundary_cells];
    VecSetValues(b, num_boundary_cells, boundary_cells, zeros, INSERT_VALUES);
    VecAssemblyEnd(b);
  }

  // Go over the boundary cells, enforcing boundary conditions on each 
  // face.
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  for (int c = 0; c < num_boundary_cells; ++c)
  {
    int bcell = boundary_cells[c];
    cell_t* cell = &mesh->cells[bcell];
    boundary_cell_t* cell_info = (boundary_cell_t*)int_ptr_unordered_map_get(boundary_cell_info, bcell);

    // Construct a polynomial least-squares fit for the cell and its 
    // neighbors (about its center), plus any boundary faces. That is,
    // construct a fit that satisfies the boundary conditions!

    // Number of points = number of neighbors + self + boundary faces
    int nb = cell_info->num_boundary_faces;
    int num_points = cell_info->num_neighbor_cells + 1 + nb;
    point_t points[num_points];
    points[0].x = mesh->cells[bcell].center.x;
    points[0].y = mesh->cells[bcell].center.y;
    points[0].z = mesh->cells[bcell].center.z;
    for (int n = 0; n < cell_info->num_neighbor_cells; ++n)
    {
      int neighbor = cell_info->neighbor_cells[n];
      points[n+1].x = mesh->cells[neighbor].center.x;
      points[n+1].y = mesh->cells[neighbor].center.y;
      points[n+1].z = mesh->cells[neighbor].center.z;
    }
    int boundary_point_indices[nb];
    for (int n = 0; n < nb; ++n)
    {
      int bface = cell_info->boundary_faces[n];
      face_t* face = &mesh->faces[bface];
      int offset = 1 + cell_info->num_neighbor_cells;
      boundary_point_indices[n] = n+offset;
      points[n+offset].x = face->center.x;
      points[n+offset].y = face->center.y;
      points[n+offset].z = face->center.z;
    }
    poly_ls_shape_set_domain(shape, &points[0], points, num_points);

    // Now traverse the boundary faces of this cell and enforce the 
    // appropriate boundary condition on each. To do this, we compute 
    // the elements of an affine linear transformation that allows us 
    // to compute boundary values of the solution in terms of the 
    // interior values.
    vector_t face_normals[nb];
    double aff_matrix[nb*num_points], aff_vector[nb];
    {
      double a[nb], b[nb], c[nb], d[nb], e[nb];
      for (int f = 0; f < nb; ++f)
      {
        // Retrieve the boundary condition for the face.
        poisson_bc_t* bc = cell_info->bc_for_face[f];

        // Compute the face normal and the corresponding coefficients,
        // enforcing alpha * phi + beta * dphi/dn = F on the face.
        face_t* face = &mesh->faces[cell_info->boundary_faces[f]];
        vector_t* n = &face_normals[f];
        n->x = face->center.x - cell->center.x;
        n->y = face->center.y - cell->center.y;
        n->z = face->center.z - cell->center.z;
        a[f] = bc->alpha;
        b[f] = bc->beta*n->x, c[f] = bc->beta*n->y, d[f] = bc->beta*n->z;

        // Compute F at the face center.
        st_func_eval(bc->F, &face->center, t, &e[f]);
      }

      // Now that we've gathered information about all the boundary
      // conditions, compute the affine transformation that maps the 
      // (unconstrained) values of the solution to the constrained values. 
      poly_ls_shape_compute_constraint_transform(shape, 
          boundary_point_indices, nb, a, b, c, d, e, aff_matrix, aff_vector);
    }

    // Compute the flux through each boundary face and alter the 
    // linear system accordingly.
    int ij[num_points];
    double N[num_points], Aij[num_points], bi;
    vector_t grad_N[num_points]; 
    for (int f = 0; f < nb; ++f)
    {
      // Compute the shape function values and gradients at the face center.
      face_t* face = &mesh->faces[cell_info->boundary_faces[f]];
      poly_ls_shape_compute_gradients(shape, &face->center, N, grad_N);

      // Add the dphi/dn terms for face f to the matrix.
      vector_t* n = &face_normals[f];
      for (int j = 0; j < num_points; ++j)
      {
        ij[j] = cell_info->neighbor_cells[j];
        Aij[j] = aff_matrix[nb*j+f]*vector_dot(n, &grad_N[j]);
        bi = -aff_vector[f];
      }
      MatSetValuesLocal(A, 1, &bcell, num_points, ij, Aij, ADD_VALUES);
      VecSetValues(b, 1, &bcell, &bi, ADD_VALUES);
    }
  }
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyEnd(b);
}

static void poisson_run_benchmark(const char* benchmark)
{
  char* variant_str = NULL;
  if ((variant_str = strstr(benchmark, "laplace_sov")) != NULL)
  {
    int variant = atoi(variant_str + strlen("laplace_sov"));
    poisson_run_laplace_sov(variant);
  }
  else if ((variant_str = strstr(benchmark, "paraboloid")) != NULL)
  {
    int variant = atoi(variant_str + strlen("paraboloid"));
    poisson_run_paraboloid(variant);
  }
  else
  {
    char err[1024];
    snprintf(err, 1024, "poisson: unknown benchmark: '%s'", benchmark);
    arbi_error(err);
  }
}

static void poisson_advance(void* context, double t, double dt)
{
  poisson_t* p = (poisson_t*)context;

  // Compute the RHS vector.
  double values[p->mesh->num_cells];
  int indices[p->mesh->num_cells];
  VecAssemblyBegin(p->b);
  for (int c = 0; c < p->mesh->num_cells; ++c)
  {
    indices[c] = c;
    point_t xc = {.x = 0.0, .y = 0.0, .z = 0.0};
    st_func_eval(p->rhs, &xc, t+dt, &values[c]);
  }
  VecSetValues(p->b, p->mesh->num_cells, indices, values, INSERT_VALUES);
  VecAssemblyEnd(p->b);

  // Make sure that boundary conditions are satisfied.
  apply_bcs(p->boundary_cells, p->num_boundary_cells, p->boundary_cell_info, 
            p->mesh, p->shape, t+dt, p->A, p->b);

  // Set up the linear solver.
  KSPSetOperators(p->solver, p->A, p->A, SAME_NONZERO_PATTERN);

  // Solve the linear system.
  KSPSolve(p->solver, p->b, p->x);

  // Copy the values from x to our solution array.
  double* x;
  VecGetArray(p->x, &x);
  memcpy(p->phi, x, sizeof(double)*p->mesh->num_cells);
  VecRestoreArray(p->x, &x);
}

static void poisson_init(void* context, double t)
{
  poisson_t* p = (poisson_t*)context;

  if (p->initialized)
  {
    KSPDestroy(&p->solver);
    MatDestroy(&p->A);
    VecDestroy(&p->x);
    VecDestroy(&p->b);
    free(p->phi);
    p->initialized = false;
  }

  // Initialize the linear solver and friends.
  MatCreate(p->comm, &p->A);
  MatSetOption(p->A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
//  MatSetOption(p->A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
//  MatSetOption(p->A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  MatSetType(p->A, MATSEQAIJ);
  MatSetSizes(p->A, p->mesh->num_cells, p->mesh->num_cells, PETSC_DETERMINE, PETSC_DETERMINE);
  VecCreate(p->comm, &p->x);
  VecSetType(p->x, VECSEQ);
  VecCreate(p->comm, &p->b);
  VecSetType(p->b, VECSEQ);
  KSPCreate(p->comm, &p->solver);

  // Initialize the solution vector.
  p->phi = malloc(sizeof(double)*p->mesh->num_cells);

  // Pre-allocate matrix entries.
  PetscInt nz = 0, nnz[p->mesh->num_cells];
  for (int i = 0; i < p->mesh->num_cells; ++i)
    nnz[i] = lin_op_stencil_size(p->L, i);
  MatSeqAIJSetPreallocation(p->A, nz, nnz);

  // Set matrix entries.
  MatAssemblyBegin(p->A, MAT_FINAL_ASSEMBLY);
  for (int i = 0; i < p->mesh->num_cells; ++i)
  {
    int indices[nnz[i]];
    double values[nnz[i]];
    lin_op_compute_stencil(p->L, i, indices, values);
    for (int j = 0; j < nnz[i]; ++j)
      indices[j] += i;
    MatSetValuesLocal(p->A, 1, &i, nnz[i], indices, values, INSERT_VALUES);
  }
  MatAssemblyEnd(p->A, MAT_FINAL_ASSEMBLY);

  // Gather information about boundary cells.
  // FIXME
//  boundary_cell_t* cell_info = (boundary_cell_t*)int_ptr_unordered_map_get(boundary_cell_info, bcell);

  // We simply solve the problem for t = 0.
  poisson_advance((void*)p, t, 0.0);

  // We are now initialized.
  p->initialized = true;
}

static void poisson_dump(void* context, io_interface_t* io, double t, int step)
{
  ASSERT(context != NULL);
  poisson_t* p = (poisson_t*)context;

  io_dataset_t* dataset = io_dataset_new("default", 1, 0);
  io_dataset_write_mesh(dataset, p->mesh);
  io_dataset_write_field(dataset, "phi", p->phi, 1, MESH_CELL);
  io_append_dataset(io, dataset);
}

static void poisson_dtor(void* ctx)
{
  poisson_t* p = (poisson_t*)ctx;
  mesh_free(p->mesh);
  int pos = 0;
  char* key;
  void* val;
  while (str_ptr_unordered_map_next(p->bcs, &pos, &key, &val))
  {
    free(key);
    free_bc(val);
  }
  free(p->boundary_cells);
  int_ptr_unordered_map_free(p->boundary_cell_info);
  free(p);
}

model_t* poisson_model_new(options_t* options)
{
  model_vtable vtable = { .run_benchmark = &poisson_run_benchmark,
                          .init = &poisson_init,
                          .advance = &poisson_advance,
                          .dump = &poisson_dump,
                          .dtor = &poisson_dtor};
  poisson_t* context = malloc(sizeof(poisson_t));
  context->phi = NULL;
  context->initialized = false;
  context->comm = MPI_COMM_WORLD;
  context->boundary_cells = NULL;
  context->num_boundary_cells = 0;
  context->boundary_cell_info = int_ptr_unordered_map_new();
  model_t* model = model_new("poisson", context, vtable);
  static const char* benchmarks[] = {"laplace_sov1", "laplace_sov2", "laplace_sov3", 
                                     "paraboloid1", "paraboloid2", "paraboloid3",
                                     NULL};
  model_register_benchmarks(model, benchmarks);
  if (options != NULL)
  {
  }
  return model;
}


#ifdef __cplusplus
}
#endif

