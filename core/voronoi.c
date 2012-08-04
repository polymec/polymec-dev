#include "core/voronoi.h"
#include "mpi.h"
#include "ctetgen.h"

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------
mesh_t* voronoi_tessellation(point_t* points, int num_points, bbox_t* bounding_box)
{
  ASSERT(points != NULL);
  ASSERT(num_points >= 2);
  ASSERT(bounding_box != NULL);

  // Not yet implemented!
  return NULL;
}
//------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

