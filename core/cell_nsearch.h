#ifndef ARBI_CELL_NSEARCH_H
#define ARBI_CELL_NSEARCH_H

#include "arbi.h"

#ifdef __cplusplus
extern "C" {
#endif

// A neighbor search algorithm for finding the "neighbor" cells of a given 
// cell. Objects of this type are garbage-collected.
typedef struct cell_nsearch_t cell_nsearch_t;

// Creates a cell_nsearch object that searches for the given minimum number 
// of neighbors for a cell.
cell_nsearch_t* cell_nsearch_new(int min_neighbors_sought);

#ifdef __cplusplus
}
#endif

#endif

