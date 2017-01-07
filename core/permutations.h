// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PERMUTE_H
#define POLYMEC_PERMUTE_H

#include "core/rng.h"

// Generates a random permutation of N elements using Knuth shuffles and the 
// given random number generator.
void random_permutation(int N, rng_t* rng, int* permutation);

#endif
