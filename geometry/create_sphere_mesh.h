#ifndef ARBI_CREATE_SPHERE_MESH_H
#define ARBI_CREATE_SPHERE_MESH_H

#include "core/mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

mesh_t* create_sphere_mesh(int n_center, int n_radial, 
                           double box_length, double radius);

#ifdef __cplusplus
}
#endif

#endif

