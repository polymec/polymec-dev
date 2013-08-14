#ifndef POLYMEC_CREATE_CVT_H
#define POLYMEC_CREATE_CVT_H

#include "core/mesh.h"
#include "core/point.h"

// This type defines a scheme for moving points around in a Centroidal 
// Voronoi Tessellation iteration.
typedef struct cvt_iterator_t cvt_iterator_t;

// To define a specific CVT iterator, the following methods must be 
// implemented:

// A function for initializing data within the context object of an iterator.
typedef void (*cvt_iterator_init_func)(void* context, 
                                       point_t* stationary_generators, 
                                       int num_stationary_generators,
                                       point_t* mobile_generators, 
                                       int num_mobile_generators);

// A function for moving a set of mobile points/generators.
typedef void (*cvt_iterator_move_points_func)(void* context, 
                                              point_t* mobile_generators, 
                                              int num_mobile_generators);

// A function that returns true if the CVT iteration is finished based 
// on a given tessellation.
typedef bool (*cvt_iterator_is_finished_func)(void* context, 
                                              mesh_t* mesh,
                                              int iteration);

// A function for destroying the context within a CVT iterator.
typedef void (*cvt_dtor)(void*);

typedef struct
{
  cvt_iterator_init_func        init;
  cvt_iterator_move_points_func move_points;
  cvt_iterator_is_finished_func is_finished;
  cvt_dtor                      dtor;
} cvt_iterator_vtable;

// Creates an instance a model with the given name and virtual table.
cvt_iterator_t* cvt_iterator_new(const char* name, void* context, cvt_iterator_vtable vtable);

// Destroys the CVT iterator.
void cvt_iterator_free(cvt_iterator_t* cvt_iter);

// Returns the name of the CVT iterator.
char* cvt_iterator_name(cvt_iterator_t* cvt_iter);

// Returns the context object for the CVT iterator.
void* cvt_iterator_context(cvt_iterator_t* cvt_iter);

// This function creates a Centroidal Voronoi Tessellation (CVT) given a set 
// of stationary and mobile points, plus an iteration method for achieving 
// the tessellation. The CVT iterator is consumed in the process.
mesh_t* create_cvt(point_t* stationary_generators, int num_stationary_generators, 
                   point_t* mobile_generators, int num_mobile_generators,
                   cvt_iterator_t* cvt_iter);

#endif

