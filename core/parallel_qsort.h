// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PARALLEL_QSORT_H
#define POLYMEC_PARALLEL_QSORT_H

#include "core/rng.h"

// This function implements a parallel quicksort using regular or 
// random sampling. Data is assumed to exist for each process belonging to 
// the communicator comm, in an array pointed to by base, each having a 
// (locally-stored) number of elements nel, each of a size width in bytes. 
// The comparator compar is used to perform comparisons between local elements 
// in each array. If rng is NULL, sampling for pivots is performed at regular 
// intervals. Otherwise, rng is used to randomly select pivots on each process.
// On output, base points to the same locally-stored array as before, but
// its data has been sorted on each process p such that processes preceding p 
// contain sorted data preceding that on p, and processes following p contain 
// sorted data following that on p.
void parallel_qsort(MPI_Comm comm, 
                    void* base, 
                    size_t nel, 
                    size_t width,
                    int (*compar)(const void* left, const void* right),
                    rng_t* rng);

#endif

