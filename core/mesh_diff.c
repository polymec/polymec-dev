#include "core/mesh_diff.h"

#ifdef __cplusplus
extern "C" {
#endif

struct mesh_diff_t
{
  int num_deltas, capacity;
  mesh_delta_t** deltas;
};

static void mesh_diff_reserve(mesh_diff_t* diff, int capacity)
{
  if (diff->capacity < capacity)
  {
    int new_size = MAX(diff->capacity, 1);
    while (new_size < capacity)
      new_size *= 2;
    diff->deltas = realloc(diff->deltas, sizeof(mesh_delta_t*)*new_size);
    diff->capacity = new_size;
  }
}

mesh_diff_t* mesh_diff_new()
{
  mesh_diff_t* diff = malloc(sizeof(mesh_diff_t));
  diff->num_deltas = 0;
  diff->capacity = 0;
  diff->deltas = NULL;
  mesh_diff_reserve(diff, 32);
  return diff;
}

void mesh_diff_free(mesh_diff_t* diff)
{
  for (int d = 0; d < diff->num_deltas; ++d)
    diff->deltas[d] = NULL;
  free(diff->deltas);
}

void mesh_diff_append(mesh_diff_t* diff, mesh_delta_t* delta)
{
  ASSERT(delta != NULL);
  mesh_diff_reserve(diff, diff->num_deltas+1);
  diff->deltas[diff->num_deltas] = delta;
  diff->num_deltas += 1;
}

void mesh_diff_apply(mesh_diff_t* diff, mesh_t* mesh)
{
  for (int d = 0; d < diff->num_deltas; ++d)
    mesh_delta_apply(diff->deltas[d], mesh);
}

bool mesh_diff_rollback(mesh_diff_t* diff, mesh_t* mesh)
{
  mesh_diff_t* inv = mesh_diff_inverse(diff);
  mesh_diff_apply(inv, mesh);
  return true; // FIXME
}

mesh_diff_t* mesh_diff_inverse(mesh_diff_t* diff)
{
  mesh_diff_t* inv = mesh_diff_new();
  for (int d = diff->num_deltas-1; d >= 0; --d)
    mesh_diff_append(inv, diff->deltas[d]);
  return inv;
}

void mesh_diff_fprintf(mesh_diff_t* diff, FILE* file)
{
  fprintf(file, "Mesh diff (%d deltas):\n", diff->num_deltas);
  for (int d = 0; d < diff->num_deltas; ++d)
  {
    fprintf(file, "  ");
    mesh_delta_fprintf(diff->deltas[d], file);
    fprintf(file, "\n");
  }
}

#ifdef __cplusplus
}
#endif

