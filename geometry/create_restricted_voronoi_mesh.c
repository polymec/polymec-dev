#include "geometry/create_restricted_voronoi_mesh.h"
#include "geometry/create_unbounded_voronoi_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

mesh_t* create_restricted_voronoi_mesh(point_set_t* generators,
                                       point_set_t* ghost_generators,
                                       surface_mesh_t* surface_mesh)
{
  // First, construct an unbounded mesh.
  mesh_t* mesh = NULL;
  {
    point_set_pos_t pos = point_set_start(generators);
    int index;
    int num_gens = point_set_size(generators);
    int num_ghost_gens = point_set_size(ghost_generators);
    point_t gens[num_gens], ghost_gens[num_ghost_gens];
    double coords[3];
    while (point_set_next(generators, &pos, &index, coords))
    {
      gens[index].x = coords[0]; 
      gens[index].y = coords[1]; 
      gens[index].z = coords[2];
    }
    pos = point_set_start(ghost_generators);
    while (point_set_next(ghost_generators, &pos, &index, coords))
    {
      ghost_gens[num_gens + index].x = coords[0];
      ghost_gens[num_gens + index].y = coords[1]; 
      ghost_gens[num_gens + index].z = coords[2];
    }
    mesh = create_unbounded_voronoi_mesh(gens, num_gens, ghost_gens, num_ghost_gens);
  }

  // FIXME

  return mesh;
}

#ifdef __cplusplus
}
#endif

