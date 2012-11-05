#include <string.h>
#include <stdlib.h>
#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"
#include "core/unordered_map.h"
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

// LAPACK prototypes.

// LU factorization.
void dgetrf(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
// Matrix inverse.
void dgetri(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

// Boundary condition structure.
typedef enum { DIRICHLET, NEUMANN, ROBIN } poisson_bc_type;
typedef struct
{
  poisson_bc_type type;
  st_func_t* func;
} poisson_bc_t;

static poisson_bc_t* create_dirichlet_bc(st_func_t* func)
{
  ASSERT(func != NULL);
  poisson_bc_t* bc = malloc(sizeof(poisson_bc_t));
  bc->type = DIRICHLET;
  bc->func = func;
  return bc;
}

static poisson_bc_t* create_neumann_bc(st_func_t* func)
{
  ASSERT(func != NULL);
  poisson_bc_t* bc = malloc(sizeof(poisson_bc_t));
  bc->type = NEUMANN;
  bc->func = func;
  return bc;
}

static poisson_bc_t* create_robin_bc(st_func_t* func)
{
  ASSERT(func != NULL);
  poisson_bc_t* bc = malloc(sizeof(poisson_bc_t));
  bc->type = ROBIN;
  bc->func = func;
  return bc;
}

static void free_bc(void* bc)
{
  poisson_bc_t* pbc = (poisson_bc_t*)bc;
  pbc->func = NULL;
}

// Context structure.
typedef struct 
{
  mesh_t* mesh;             // Mesh.
  st_func_t* rhs;           // Right-hand side function.
  str_ptr_unordered_map_t* bcs;  // Boundary conditions.
  double* phi;              // Solution array.
  lin_op_t* L;              // Laplacian operator.

  KSP solver;               // Linear solver.
  Mat A;                    // Laplacian matrix.
  Vec x;                    // Solution vector.
  Vec b;                    // RHS vector. 

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
  p->phi = NULL;
  p->L = laplacian_op_new(p->mesh);
  p->comm = MPI_COMM_WORLD;
  p->initialized = false;
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
    str_ptr_unordered_map_insert(bcs, "boundary", create_dirichlet_bc(NULL));
  else if ((variant % 3) == 1)
    str_ptr_unordered_map_insert(bcs, "boundary", create_neumann_bc(NULL));
  else // if ((variant % 3) == 2)
    str_ptr_unordered_map_insert(bcs, "boundary", create_robin_bc(NULL));

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
    str_ptr_unordered_map_insert(bcs, "boundary", create_dirichlet_bc(NULL));
  else if ((variant % 3) == 1)
    str_ptr_unordered_map_insert(bcs, "boundary", create_neumann_bc(NULL));
  else // if ((variant % 3) == 2)
    str_ptr_unordered_map_insert(bcs, "boundary", create_robin_bc(NULL));

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

static void apply_dirichlet_bc(mesh_t* mesh, int* boundary_faces, int num_boundary_faces, 
                               st_func_t* func, double t, Mat A, Vec b)
{
  // Make a list of cells adjacent to boundary faces.
  int_unordered_set_t* bcells = int_unordered_set_new();
  // FIXME

  // Traverse the boundary cells and generate replacement rows in the 
  // linear system that enforce the Dirichlet boundary condition.
  int pos = 0, bcell;
  while (int_unordered_set_next(bcells, &pos, &bcell))
  {
    // Generate a list of neighboring cells for this cell, plus a list of 
    // boundary faces attached to it.
    // FIXME

    // Generate the moment matrix corresponding to a least-squares linear 
    // polynomial fit of the solution about this cell.
    // FIXME

    // Factor the moment matrix and compute its inverse.
    // FIXME

    // Compute the vector sum of the face-normal-area products for boundary 
    // faces attached to this cell. This will be subtracted from the RHS.
    // FIXME

    // Compute the vector sum of the face-normal-area products for interior
    // faces attached to this cell.
    // FIXME

    // Generate the row in the matrix A that corresponds to this cell.
    // FIXME

    // Subtract off the boundary contribution from the RHS.
    // FIXME
  }

  // Clean up.
  int_unordered_set_free(bcells);
}

static void apply_neumann_bc(mesh_t* mesh, int* boundary_faces, int num_boundary_faces,
                             st_func_t* func, double t, Mat A, Vec b)
{
  int indices[num_boundary_faces];
  double values[num_boundary_faces];
  for (int f = 0; f < num_boundary_faces; ++f)
  {
    int bface = boundary_faces[f];
    face_t* face = &mesh->faces[bface];
    double bvalue;
    st_func_eval(func, &face->center, t, &bvalue);
    bvalue = -bvalue * face->area;
    indices[f] = bface;
    values[f] = bvalue;
  }
  VecAssemblyBegin(b);
  VecSetValues(b, mesh->num_cells, indices, values, ADD_VALUES);
  VecAssemblyEnd(b);
}

static void apply_robin_bc(mesh_t* mesh, int* boundary_faces, int num_boundary_faces,
                           st_func_t* func, double t, Mat A, Vec b)
{
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
  int pos = 0;
  char* key;
  void* val;
  while (str_ptr_unordered_map_next(p->bcs, &pos, &key, &val))
  {
    // We enforce boundary conditions with ghost cells. For each face on the 
    // boundary, cell1 is the interior cell and cell2 is the ghost cell.
    int num_faces;
    int* boundary_faces = mesh_tag(p->mesh->face_tags, (const char*)key, &num_faces);
    poisson_bc_t* bc = (poisson_bc_t*)val;

    if (bc->type == DIRICHLET)
      apply_dirichlet_bc(p->mesh, boundary_faces, num_faces, bc->func, t+dt, p->A, p->b);
    else if (bc->type == NEUMANN)
      apply_neumann_bc(p->mesh, boundary_faces, num_faces, bc->func, t+dt, p->A, p->b);
    else 
    {
      ASSERT(bc->type == ROBIN);
      apply_robin_bc(p->mesh, boundary_faces, num_faces, bc->func, t+dt, p->A, p->b);
    }
  }

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

