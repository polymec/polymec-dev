/* ==========================================================================
 * libarena/util.h - Custom Memory Allocator Interface
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
#ifndef ARENA_UTIL_H
#define ARENA_UTIL_H


/*
 * Shortcuts to the arena_util API.
 */
#ifndef ARENA_UTIL_NOAP

#define ap_strldup(a, b, c)	arena_util_strldup(a, b, c)

#define ap_strndup(a, b, c)	arena_util_strndup(a, b, c)

#define ap_strdup(a, b)		arena_util_strdup(a, b)

#define ap_memdup(a, b, c)	arena_util_memdup(a, b, c)

#define ap_vasprintf(a, b, c, d) \
				arena_util_vasprintf(a, b, c, d)

#define ap_asprintf(a, b, c, ...) \
				arena_util_asprintf(a, b, c, __VA_ARGS__)

#define ap_vsprintf(a, b, c)	arena_util_vsprintf(a, b, c)

#define ap_sprintf(a, b, ...)	arena_util_sprintf(a, b, __VA_ARGS__)

#define ap_split(a, b, c, d, e)	arena_util_split(a, b, c, d, e)

#endif /* ARENA_UTIL_NOAP */


char *arena_util_strldup(const struct arena_prototype *a, const char *src, size_t srclen);

char *arena_util_strndup(const struct arena_prototype *a, const char *src, size_t maxlen);

char *arena_util_strdup(const struct arena_prototype *a, const char *src);

void *arena_util_memdup(const struct arena_prototype *a, const void *src, size_t srclen);

#ifndef __GNUC__
#ifndef __attribute__
#define __attribute__(x)
#endif
#endif

int arena_util_vasprintf(const struct arena_prototype *a, char **dstp, const char *fmt, va_list ap)
	__attribute__((__format__ (printf, 3, 0)));

int arena_util_asprintf(const struct arena_prototype *a, char **dstp, const char *fmt, ...)
	__attribute__((__format__ (printf, 3, 4)));

char *arena_util_vsprintf(const struct arena_prototype *a, const char *fmt, va_list ap)
	__attribute__((__format__ (printf, 2, 0)));

char *arena_util_sprintf(const struct arena_prototype *a, const char *fmt, ...)
	__attribute__((__format__ (printf, 2, 3)));

#if HAVE_STRSEP
int arena_util_split(const struct arena_prototype *a, int *argcp, char ***argvp, char **strp, const char *sep);
#endif


#endif /* ARENA_UTIL_H */
