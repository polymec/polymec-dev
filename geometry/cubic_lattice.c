#include "core/polymec.h"
#include "geometry/cubic_lattice.h"
#include <gc/gc.h>

#ifdef __cplusplus
extern "C" {
#endif

cubic_lattice_t* cubic_lattice_new(int nx, int ny, int nz)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  cubic_lattice_t* l = GC_MALLOC(sizeof(cubic_lattice_t));
  l->nx = nx;
  l->ny = ny;
  l->nz = nz;
  return l;
}

int_int_unordered_map_t* cubic_lattice_generate_x_periodic_map(void* context, mesh_t* mesh, char* tag1, char* tag2)
{
  cubic_lattice_t* lattice = mesh_property(mesh, "lattice");
  ASSERT(lattice != NULL);
  int_int_unordered_map_t* map = int_int_unordered_map_new();
 
  // Go over the x-faces and map everything.
  for (int j = 0; j < lattice->ny; ++j)
  {
    for (int k = 0; k < lattice->nz; ++k)
    {
      int low = cubic_lattice_x_face(lattice, 0, j, k);
      int high = cubic_lattice_x_face(lattice, lattice->nx, j, k);
      int_int_unordered_map_insert(map, low, high);
      int_int_unordered_map_insert(map, high, low);
    }
  }

  return map;
}

int_int_unordered_map_t* cubic_lattice_generate_y_periodic_map(void* context, mesh_t* mesh, char* tag1, char* tag2)
{
  cubic_lattice_t* lattice = mesh_property(mesh, "lattice");
  ASSERT(lattice != NULL);
  int_int_unordered_map_t* map = int_int_unordered_map_new();
 
  // Go over the y-faces and map everything.
  for (int i = 0; i < lattice->nx; ++i)
  {
    for (int k = 0; k < lattice->nz; ++k)
    {
      int low = cubic_lattice_y_face(lattice, i, 0, k);
      int high = cubic_lattice_y_face(lattice, i, lattice->ny, k);
      int_int_unordered_map_insert(map, low, high);
      int_int_unordered_map_insert(map, high, low);
    }
  }

  return map;
}

int_int_unordered_map_t* cubic_lattice_generate_z_periodic_map(void* context, mesh_t* mesh, char* tag1, char* tag2)
{
  cubic_lattice_t* lattice = mesh_property(mesh, "lattice");
  ASSERT(lattice != NULL);
  int_int_unordered_map_t* map = int_int_unordered_map_new();
 
  // Go over the z-faces and map everything.
  for (int i = 0; i < lattice->nx; ++i)
  {
    for (int j = 0; j < lattice->ny; ++j)
    {
      int low = cubic_lattice_z_face(lattice, i, j, 0);
      int high = cubic_lattice_z_face(lattice, i, j, lattice->nz);
      int_int_unordered_map_insert(map, low, high);
      int_int_unordered_map_insert(map, high, low);
    }
  }

  return map;
}

#ifdef __cplusplus
}
#endif

