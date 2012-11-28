#ifndef POLYMEC_CREATE_BOX_MESH_H
#define POLYMEC_CREATE_BOX_MESH_H

#include "core/mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

mesh_t* create_box_mesh(int N[], double low[], double high[]);

#ifdef __cplusplus
}
#endif

#endif

