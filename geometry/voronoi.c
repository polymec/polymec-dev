#include "geometry/voronoi.h"
#include "mpi.h"
#include "ctetgen.h"

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------
mesh_t* voronoi_tessellation(point_t* points, int num_points, 
                             point_t* ghost_points, int num_ghost_points)
{
  ASSERT(points != NULL);
  ASSERT(num_points >= 2);
  ASSERT(num_ghost_points >= 0);

  // Not yet implemented!
  return NULL;
}
//------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

