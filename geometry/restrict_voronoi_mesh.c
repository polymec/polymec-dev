#include "core/point_set.h"
#include "geometry/restrict_voronoi_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

void restrict_voronoi_mesh(mesh_t* mesh, surface_mesh_t* surface_mesh)
{
  // Look for a point set containing the generators within the mesh.
  point_set_t* generators = mesh_property(mesh, "generators");
  ASSERT(generators != NULL);
  // FIXME
}

#ifdef __cplusplus
}
#endif

