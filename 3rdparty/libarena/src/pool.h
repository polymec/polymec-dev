/* ==========================================================================
 * libarena/pool.h.in - Custom Memory Allocator Interface
 * --------------------------------------------------------------------------
 * Copyright (c) 2006  William Ahern
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to permit
 * persons to whom the Software is furnished to do so, subject to the
 * following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
 * NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 * ==========================================================================
 */
#ifndef ARENA_POOL_H
#define ARENA_POOL_H

#include <stdarg.h>	/* Be helpful with va_list */


/*
 * Don't require arena/proto.h
 */
struct arena_prototype;


typedef struct pool POOL;

struct pool_bucket_options {
	size_t bucketlen;
	size_t nperblock;
};

extern const struct pool_options {
	size_t alignment;

	const struct pool_bucket_options *buckets;
} pool_defaults;


POOL *pool_open(const struct pool_options *, const struct arena_prototype *);

void pool_close(POOL *);

const struct arena_prototype *pool_export(POOL *);

void *pool_get(POOL *, size_t, size_t);

void pool_put(POOL *, void *);

void *pool_realloc(POOL *, void *, size_t, size_t);

struct pool *pool_import(const struct arena_prototype *);

const char *pool_strerror(POOL *);

void pool_clearerr(POOL *);

char *pool_strdup(POOL *, const char *);

char *pool_strndup(POOL *, const char *, size_t);

void *pool_memdup(POOL *, const void *, size_t);

#ifndef __GNUC__
#ifndef __attribute__
#define __attribute__(x)
#endif
#endif

int pool_vasprintf(POOL *, char **, const char *, va_list)
	__attribute__((__format__ (printf, 3, 0)));

int pool_asprintf(POOL *, char **, const char *, ...)
	__attribute__((__format__ (printf, 3, 4)));

char *pool_vsprintf(POOL *, const char *, va_list)
	__attribute__((__format__ (printf, 2, 0)));

char *pool_sprintf(POOL *, const char *, ...)
	__attribute__((__format__ (printf, 2, 3)));

#endif /* ARENA_POOL_H */
