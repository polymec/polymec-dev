#ifndef POLYMEC_MERGE_MESH_NODES_H
#define POLYMEC_MERGE_MESH_NODES_H

#include "core/mesh.h"

// This merges any nodes within the given mesh that are closer together 
// than the given tolerance. This is not a reversible operation.
void merge_mesh_nodes(mesh_t* mesh, double tolerance);

#endif

