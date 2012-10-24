#include "poisson/laplacian_op.h"

#ifdef __cplusplus
extern "C" {
#endif

static int laplacian_op_stencil_size(void* context, mesh_t* mesh, int index)
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

static void laplacian_op_compute_offsets(void* context, mesh_t* mesh, int index, int* offsets)
{
  int i = 0;
  offsets[i++] = 0; // "diagonal" term
  cell_t* cell = &mesh->cells[index];
  for (int f = 0; f < cell->num_faces; ++f)
  {
    cell_t* opp_cell = face_opp_cell(cell->faces[f], cell);
    if (opp_cell != NULL)
      offsets[i++] = (opp_cell - &mesh->cells[0]) - index;
  }
}

static void laplacian_op_compute_weights(void* context, mesh_t* mesh, int index, int* offsets, double* weights)
{
  int i = 0;
  weights[i++] = 0.0; // "diagonal" term
  cell_t* cell = &mesh->cells[index];
  for (int f = 0; f < cell->num_faces; ++f)
  {
    cell_t* opp_cell = face_opp_cell(cell->faces[f], cell);
    if (opp_cell != NULL)
    {
      weights[i++] = 1.0;
      weights[0] -= 1.0;
    }
  }
}

lin_op_t* laplacian_op_new(mesh_t* mesh)
{
  lin_op_vtable vtable = {.stencil_size = &laplacian_op_stencil_size,
                          .compute_offsets = &laplacian_op_compute_offsets,
                          .compute_weights = &laplacian_op_compute_weights};
  return lin_op_new("laplacian", NULL, vtable, mesh);
}

#ifdef __cplusplus
}
#endif

