#ifndef POLYMEC_ARCH_H
#define POLYMEC_ARCH_H

#include <stdio.h>

// This file contains functions that are unavailable on certain target 
// architectures.

#ifdef APPLE

// This is a port of the fmemopen function (available on Linux).
FILE* fmemopen(void *buf, size_t size, const char *mode);

// This is a port of open_memstream (available on Linux).
FILE* open_memstream(char **buf, size_t *len);

#endif 
#endif
