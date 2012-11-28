#ifndef POLYMEC_REPARTITION_H
#define POLYMEC_REPARTITION_H

#include "point.h"
#include "exchanger.h"
#include "mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function repartitions the given set of points with the given 
// weights, alloting them to parallel domains to balance their load.
// If weights is NULL, the points are assigned equal weights.
// It creates and returns an exchanger object that can be used to migrate 
// data from the old partition to the new.
exchanger_t* repartition_points(MPI_Comm comm, point_t* points, double* weights, int num_points);

// This function repartitions the given mesh, alloting the cells to parallel 
// domains to balance the load.
// It creates and returns an exchanger object that can be used to migrate 
// data from the old partition to the new.
exchanger_t* repartition_mesh(MPI_Comm comm, mesh_t* mesh);

#ifdef __cplusplus
}
#endif

#endif

