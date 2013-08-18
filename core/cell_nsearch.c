// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

