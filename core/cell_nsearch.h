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

#ifndef POLYMEC_CELL_NSEARCH_H
#define POLYMEC_CELL_NSEARCH_H

#include "polymec.h"

// A neighbor search algorithm for finding the "neighbor" cells of a given 
// cell. Objects of this type are garbage-collected.
typedef struct cell_nsearch_t cell_nsearch_t;

// Creates a cell_nsearch object that searches for the given minimum number 
// of neighbors for a cell.
cell_nsearch_t* cell_nsearch_new(int min_neighbors_sought);

#endif

