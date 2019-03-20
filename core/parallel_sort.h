// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PARALLEL_SORT_H
#define POLYMEC_PARALLEL_SORT_H

#include "core/polymec.h"

/// \addtogroup core core
///@{

/// This function sorts an array whose data is distributed across the given
/// MPI communicator, using an odd-even transposition sort.
/// array pointed to by base, each having a (locally-stored) number of elements
/// The comparator comp is used to perform comparisons between local elements
/// in each array. On output, base points to the same locally-stored array as
/// before, but its data has been sorted on each process p such that processes
/// preceding p contain sorted data preceding that on p, and processes
/// following p contain sorted data following that on p.
/// \collective Collective on comm.
void parallel_sort(MPI_Comm comm,
                   void* base,
                   size_t nel,
                   size_t width,
                   int (*comp)(const void* left, const void* right));

///@}

#endif

