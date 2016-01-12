// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_DENSE_LOCAL_MATRIX_H
#define POLYMEC_DENSE_LOCAL_MATRIX_H

#include "core/local_matrix.h"

// This returns a (square) dense local matrix of the given size N.
local_matrix_t* dense_local_matrix_new(int N);

#endif
