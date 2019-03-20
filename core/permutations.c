// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/permutations.h"

void random_permutation(int N, rng_t* rng, int* permutation)
{
  for (int i = 0; i < N; ++i)
  {
    int j = rng_uniform_int(rng, N+1);
    permutation[j] = i;
  }
}

