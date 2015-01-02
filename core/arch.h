// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_ARCH_H
#define POLYMEC_ARCH_H

#include <stdio.h>

// This file contains functions that are unavailable on certain target 
// architectures.

#ifdef APPLE

#include <pthread.h>

// This is a port of the fmemopen function (available on Linux).
FILE* fmemopen(void *buf, size_t size, const char *mode);

// This is a port of open_memstream (available on Linux).
FILE* open_memstream(char **buf, size_t *len);

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
