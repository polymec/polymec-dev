#ifndef POLYMEC_PERIODIC_BC_H
#define POLYMEC_PERIODIC_BC_H

#include "core/mesh.h"
#include "core/unordered_map.h"

#ifdef __cplusplus
extern "C" {
#endif

// This is an opaque type that represents a periodic boundary condition.
// Objects of this type are garbage-collected.
typedef struct periodic_bc_t periodic_bc_t;

// Creates a periodic boundary condition that identifies faces belonging 
// to tag1 with corresponding faces belonging to tag2. It is not necessary 
// that tag1[i] corresponds to tag2[i]--the geometric information about the 
// faces in the mesh will determine the correspondence.
periodic_bc_t* periodic_bc_new(const char* tag1, const char* tag2);

// Creates a periodic boundary condition with a function (and a context) 
// for generating periodic mappings between faces.
typedef int_int_unordered_map_t* (*periodic_bc_mapping_func)(void*, mesh_t*, char*, char*);
periodic_bc_t* periodic_bc_new_with_map_func(const char* tag1, const char* tag2, periodic_bc_mapping_func mapping_func, void* mapping_context);

// Returns true if the given periodic_bc object is actually a periodic_bc object,
// false if it has been improperly cast to one.
bool periodic_bc_is_valid(periodic_bc_t* bc);

// Returns true if the given pointer can be cast to a periodic_bc object, 
// false if not.
bool pointer_is_periodic_bc(void* ptr);

// Retrieves the two tags from the periodic boundary condition. tag1 and 
// tag2 will be pointed to internally-stored strings (which must not be deleted).
void periodic_bc_get_tags(periodic_bc_t* bc, char** tag1, char** tag2);

// Generates a face-to-face mapping for a periodic boundary condition.
int_int_unordered_map_t* periodic_bc_generate_map(periodic_bc_t* bc, mesh_t* mesh);

#ifdef __cplusplus
}
#endif

#endif
