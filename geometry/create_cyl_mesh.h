#ifndef POLYMEC_CREATE_CYL_MESH_H
#define POLYMEC_CREATE_CYL_MESH_H

#include "core/mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

mesh_t* create_cyl_mesh(int n_center, int n_radius, int n_axial, 
                        double box_len, double radius, double axial_len);

#ifdef __cplusplus
}
#endif

#endif

