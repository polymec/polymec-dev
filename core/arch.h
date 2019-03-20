// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ARCH_H
#define POLYMEC_ARCH_H

#include <stdio.h>

// This file contains functions that are unavailable on certain target
// architectures.

// This is a port of the fmemopen function, which is not part of the C standard.
FILE* fmemopen(void *buf, size_t size, const char *mode);

// This is a port of open_memstream, which is not part of the C standard.
FILE* open_memstream(char **buf, size_t *len);

#endif
