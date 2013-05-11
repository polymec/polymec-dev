#include <gc/gc.h>
#include <stdlib.h>
#include "core/cell_nsearch.h"

// A neighbor search algorithm for finding the "neighbor" cells of a given 
// cell. Objects of this type are garbage-collected.
struct cell_nsearch_t 
{
  int min_neighbors;
};

static void cell_nsearch_free(void* ctx, void* dummy)
{
  cell_nsearch_t* ns = (cell_nsearch_t*)ctx;
  free(ns);
}

cell_nsearch_t* cell_nsearch_new(int min_neighbors_sought)
{
  cell_nsearch_t* ns = GC_MALLOC(sizeof(cell_nsearch_t));
  ns->min_neighbors = min_neighbors_sought;
  GC_register_finalizer(ns, &cell_nsearch_free, ns, NULL, NULL);
  return ns;
}

