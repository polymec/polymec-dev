#include "advect/advect_diffusion_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

// Advective diffusion solver type
typedef struct
{
  lin_op_t* diff_op;
  mesh_t* mesh;
} ad_solver_t;

static void ad_create_matrix(void* context, Mat* mat)
{
}

static void ad_create_vector(void* context, Vec* vec)
{
}

static void ad_create_ksp(void* context, KSP* ksp)
{
}

static void ad_compute_diffusion_matrix(void* context, Mat D, double t)
{
}

static void ad_compute_source_vector(void* context, Vec S, double t)
{
}

static void ad_apply_bcs(void* context, Mat A, Vec b, double t)
{
}

static void ad_dtor(void* context)
{
}

diffusion_solver_t* advect_diffusion_solver_new(lin_op_t* diffusion_op,
                                                mesh_t* mesh)
{
  diffusion_solver_vtable vtable = 
  { .create_matrix            = ad_create_matrix,
    .create_vector            = ad_create_vector,
    .create_ksp               = ad_create_ksp,
    .compute_diffusion_matrix = ad_compute_diffusion_matrix,
    .compute_source_vector    = ad_compute_source_vector,
    .apply_bcs                = ad_apply_bcs,
    .dtor                     = ad_dtor
  };
  ad_solver_t* a = malloc(sizeof(ad_solver_t));
  a->diff_op = diffusion_op;
  a->mesh = mesh;
  diffusion_solver_t* solver = diffusion_solver_new("advective diffusion solver", a, vtable);
  return solver;
}

#ifdef __cplusplus
}
#endif

