#include <stdlib.h>
#include "geometry/create_cvt.h"
#include "geometry/voronoi.h"
#include "tao.h"

#ifdef __cplusplus
extern "C" {
#endif

void cvt_energy_function(mesh_t* mesh, sp_func_t* density, sp_func_t* boundary,
                         double* F, double* gradF)
{
}

// Context for TAO's optimization.
typedef struct 
{
  point_t* points;
  int num_points;
  mesh_t* mesh;
  cvt_obj_func_t F;
  sp_func_t* rho;
  sp_func_t* B;
} cvt_opt_context_t;

// This computes the function F and its gradient in the parlance of TAO.
static PetscErrorCode compute_F_and_gradF(TaoSolver solver, 
                                          Vec X, 
                                          PetscReal* F,
                                          Vec gradF,
                                          void* context)
{
  cvt_opt_context_t* data = (cvt_opt_context_t*)context;

  // Copy local data from X to the points.
  // FIXME

  // Figure out the ghost points if we're in parallel.
  int Nghost = 0;
  point_t* ghost_points = NULL;
  // FIXME

  // Re-tessellate.
  if (data->mesh != NULL)
    mesh_free(data->mesh);
  data->mesh = voronoi_tessellation(data->points, data->num_points, 
                                    ghost_points, Nghost);

  // Evaluate the objective function.
  double Fval;
  double* gradFval; // FIXME
  data->F(data->mesh, data->rho, data->B, &Fval, gradFval);

  // Clean up.
  if (ghost_points != NULL)
    free(ghost_points);

  return -1;
}

mesh_t* create_bounded_cvt(int N, 
                           sp_func_t* rho, 
                           cvt_obj_func_t F, 
                           sp_func_t* B)
{
  // Generate N points at random within the domain.
  unsigned int seed = 10;
  srandom(seed);
  point_t* points = malloc(N*sizeof(point_t));
  for (int i = 0; i < N; ++i)
  {
    // Assign random components.
    points[i].x = random(), points[i].y = random(), points[i].z = random();

    // Try again if we're out of bounds.
    double Bi;
    sp_func_eval(B, &(points[i]), &Bi);
    while (Bi > 0.0)
    {
      points[i].x = random(), points[i].y = random(), points[i].z = random();
      sp_func_eval(B, &(points[i]), &Bi);
    } 
  }

  // Allocate storage for the solution and populate it with initial data.
  Vec X; 
  VecCreateSeq(PETSC_COMM_SELF, 3*N, &X);
  // FIXME

  // Set up a context for Tao.
  cvt_opt_context_t data = {.points = points, .num_points = N, .mesh = NULL, 
                            .F = F, .rho = rho, .B = B};

  // Set up a Tao solver.
  TaoSolver optimizer;
  TaoCreate(PETSC_COMM_SELF, &optimizer);
  TaoSetInitialVector(optimizer, X);
  TaoSetObjectiveAndGradientRoutine(optimizer, &compute_F_and_gradF, &data);

  // Do some nonlinear optimization!
  int status = TaoSolve(optimizer);
  if (status != 0)
    arbi_error("Optimization failed!");

  // Clean up.
  TaoDestroy(&optimizer);
  VecDestroy(&X);
  free(data.points);

  // Pluck the mesh from our context.
  return data.mesh;
}

#ifdef __cplusplus
}
#endif

