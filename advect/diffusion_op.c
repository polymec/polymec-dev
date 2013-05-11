#include "advect/diffusion_op.h"

typedef struct
{
  st_func_t* diffusivity;
  double t; // Simulation time.
} diffusion_op_t;

static void diffusion_op_dtor(void* context)
{
  diffusion_op_t* diff = (diffusion_op_t*)context;
  diff->diffusivity = NULL;
  free(diff);
}

static int diffusion_op_stencil_size(void* context, mesh_t* mesh, int index)
{
  int size = 1;
  cell_t* cell = &mesh->cells[index];
  for (int f = 0; f < cell->num_faces; ++f)
  {
    if (face_opp_cell(cell->faces[f], cell) != NULL)
      ++size;
  }
  return size;
}

static void diffusion_op_compute_stencil(void* context, mesh_t* mesh, int index, int* offsets, double* weights)
{
  diffusion_op_t* diff = (diffusion_op_t*)context;

  // "Diagonal" term
  offsets[0] = 0;
  weights[0] = 0.0;
  cell_t* cell = &mesh->cells[index];
  int i = 1;
  for (int f = 0; f < cell->num_faces; ++f)
  {
    face_t* face = cell->faces[f];
    double A = face->area;
    cell_t* opp_cell = face_opp_cell(face, cell);
    if (opp_cell != NULL)
    {
      offsets[i] = (opp_cell - &mesh->cells[0]) - index;
      double L = point_distance(&opp_cell->center, &cell->center);
      double D;
      st_func_eval(diff->diffusivity, &face->center, diff->t, &D);
      weights[i] = D*A/L;
      weights[0] -= D*A/L;
      ++i;
    }
  }
}

static void diffusion_op_apply(void* context, mesh_t* mesh, double* field, double* Lfield)
{
  diffusion_op_t* diff = (diffusion_op_t*)context;
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    cell_t* cell = &mesh->cells[c];
    Lfield[c] = 0.0;
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
      double A = face->area;
      cell_t* opp_cell = face_opp_cell(face, cell);
      if (opp_cell != NULL)
      {
        int cc = opp_cell - &mesh->cells[0];
        double L = point_distance(&opp_cell->center, &cell->center);
        double D;
        st_func_eval(diff->diffusivity, &face->center, diff->t, &D);
        Lfield[c] += D * A/L * (field[cc] - field[c]);
      }
    }
  }
}

lin_op_t* diffusion_op_new(mesh_t* mesh, st_func_t* diffusivity)
{
  ASSERT(diffusivity != NULL);
  lin_op_vtable vtable = {.stencil_size = &diffusion_op_stencil_size,
                          .compute_stencil = &diffusion_op_compute_stencil,
                          .apply = &diffusion_op_apply,
                          .dtor = &diffusion_op_dtor};
  diffusion_op_t* diff = malloc(sizeof(diffusion_op_t));
  diff->t = 0.0;
  diff->diffusivity = diffusivity;
  return lin_op_new("diffusion", (void*)diff, vtable, mesh);
}

void diffusion_op_set_time(lin_op_t* op, double t)
{
  diffusion_op_t* diff = lin_op_context(op);
  diff->t = t;
}

