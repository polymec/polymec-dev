/* ==========================================================================
 * arena/arena.h.in - Custom Memory Allocator Interface
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
#ifndef ARENA_ARENA_H
#define ARENA_ARENA_H


#include <stdarg.h>	/* Be helpful with va_list */


typedef struct arena ARENA;


extern const struct arena_options {
	size_t alignment;
	size_t blocklen;
} arena_defaults;


ARENA *arena_open(const struct arena_options *, const struct arena_prototype *);

void arena_close(ARENA *);

const struct arena_prototype *arena_export(ARENA *);

void *arena_malloc(ARENA *, size_t, size_t);

void *arena_realloc(ARENA *, void *, size_t, size_t);

void arena_free(ARENA *, void *);

void arena_mark(ARENA *, void **);

void arena_reset(ARENA *, void *);

ARENA *arena_import(const struct arena_prototype *);

const char *arena_strerror(ARENA *);

void arena_clearerr(ARENA *);

char *arena_strdup(ARENA *, const char *);

char *arena_strndup(ARENA *, const char *, size_t);

void *arena_memdup(ARENA *, const void *, size_t);

#ifndef __GNUC__
#ifndef __attribute__
#define __attribute__(x)
#endif
#endif

int arena_vasprintf(ARENA *, char **, const char *, va_list)
	__attribute__((__format__ (printf, 3, 0)));

int arena_asprintf(ARENA *, char **, const char *, ...)
	__attribute__((__format__ (printf, 3, 4)));

char *arena_vsprintf(ARENA *, const char *, va_list)
	__attribute__((__format__ (printf, 2, 0)));

char *arena_sprintf(ARENA *, const char *, ...)
	__attribute__((__format__ (printf, 2, 3)));

#endif /* ARENA_ARENA_H */

