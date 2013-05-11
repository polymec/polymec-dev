#ifndef POLYMEC_MESH_DIFF_H
#define POLYMEC_MESH_DIFF_H

#include "core/mesh_delta.h"

// A mesh_diff is a sequence of mesh_delta objects that form an atomic 
// transaction implementing a mesh update.
typedef struct mesh_diff_t mesh_diff_t;

// Creates a new mesh_diff containing no deltas.
mesh_diff_t* mesh_diff_new();

// Destroys a mesh_diff and all of its deltas.
void mesh_diff_free(mesh_diff_t* diff);

// Appends the given delta to the mesh_diff.
void mesh_diff_append(mesh_diff_t* diff, mesh_delta_t* delta);

// Applies the mesh_diff to the mesh, enacting all of the deltas.
void mesh_diff_apply(mesh_diff_t* diff, mesh_t* mesh);

// Rolls back this transaction on the mesh, assuming it was the last one 
// applied. Returns true if the rollback succeeded, false if not.
bool mesh_diff_rollback(mesh_diff_t* diff, mesh_t* mesh);

// Constructs a mesh_diff representing the inverse of this diff. This is 
// used in the rollback mechanism.
mesh_diff_t* mesh_diff_inverse(mesh_diff_t* diff);

// Writes a text representation of the mesh_diff to the given file.
void mesh_diff_fprintf(mesh_diff_t* diff, FILE* file);

#endif

