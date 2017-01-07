// Copyright (c) 2012-2017, Jeffrey N. Johnson
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

#include <pthread.h>

// This is a port of the fmemopen function, which is not part of the C standard.
FILE* fmemopen(void *buf, size_t size, const char *mode);

// This is a port of open_memstream, which is not part of the C standard.
FILE* open_memstream(char **buf, size_t *len);

#ifdef APPLE
//------------------------------------------------------------------------
// Below is an implementation of the non-standard PThreads barrier type.
//------------------------------------------------------------------------

// Descriptor for PThreads barrier attributes -- not used in this implementation.
typedef int pthread_barrierattr_t;

// Barrier type.
typedef struct pthread_barrier_t pthread_barrier_t;

// Initializes a barrier that blocks the given number of threads.
int pthread_barrier_init(pthread_barrier_t *barrier, 
                         const pthread_barrierattr_t *attr, 
                         unsigned int count);

// Destroys a barrier created with pthread_barrier_init().
int pthread_barrier_destroy(pthread_barrier_t *barrier);

// Causes the barrier to block until all of its threads finish their tasks.
int pthread_barrier_wait(pthread_barrier_t *barrier);

#endif
#endif
