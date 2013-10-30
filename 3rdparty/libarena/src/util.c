/* ==========================================================================
 * libarena/src/util.c - Custom Memory Allocator Interface
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
#include <stdio.h>	/* vsnprintf(3) */
#include <stdarg.h>	/* va_list va_copy va_start va_end */

#include <string.h>	/* memchr(3) memcpy(3) strsep(3) strerror(3) */

#include <errno.h>	/* EFAULT errno */

#include "proto.h"	/* struct arena_prototype */
#include "util.h"


#ifndef va_copy
#ifdef __va_copy
#define va_copy(a, b)	__va_copy((a), (b)
#else
#define va_copy(a, b)	((a) = (b))
#endif
#endif


/*
 * Intermediate stack buffers
 */
#ifndef ARENA_UTIL_BUFSZ
#define ARENA_UTIL_BUFSZ	256
#endif


char *arena_util_strldup(const struct arena_prototype *a, const char *src, size_t srclen) {
	char *dst;

	if (!(dst = a->malloc(a,srclen + 1,1)))
		return 0;

	dst[srclen]	= '\0';

	return memcpy(dst,src,srclen);
} /* arena_util_strldup() */


char *arena_util_strndup(const struct arena_prototype *a, const char *src, size_t maxlen) {
	char *nul;
	size_t srclen;

	if ((nul = memchr(src,'\0',maxlen)))
		srclen	= nul - src;
	else
		srclen	= maxlen;

	return arena_util_strldup(a,src,srclen);
} /* arena_util_strndup() */


char *arena_util_strdup(const struct arena_prototype *a, const char *src) {
	size_t srclen	= strlen(src);
	char *dst;

	if (!(dst = a->malloc(a,srclen + 1,1)))
		return 0;

	(void)memcpy(dst,src,srclen);

	dst[srclen]	= '\0';

	return dst;
} /* arena_util_strndup() */


void *arena_util_memdup(const struct arena_prototype *a, const void *src, size_t srclen) {
	void *dst;

	if ((dst = a->malloc(a,srclen,0)))
		return memcpy(dst,src,srclen);
	else
		return 0;
} /* arena_util_memdup() */


int arena_util_vasprintf(const struct arena_prototype *a, char **dstp, const char *fmt, va_list ap) {
	char buf[ARENA_UTIL_BUFSZ];
	int len;
	va_list tap;

	*dstp	= 0;

	va_copy(tap, ap);

	len	= vsnprintf(buf,sizeof buf,fmt,tap);

	if (len < 0)
		return -1;

	if (!(*dstp = a->malloc(a,len + 1,1)))
		return -1;

	if ((size_t)len < sizeof buf)
		return memcpy(*dstp,buf,len + 1), len;

	va_copy(tap,ap);

	if (len	!= vsnprintf(*dstp,len + 1,fmt,tap)) {
		int errsav	= errno;

		a->free(a,*dstp);

		*dstp	= 0;

		errno	= errsav;

		return -1;
	}

	return len;
} /* arena_util_vasprintf() */


int arena_util_asprintf(const struct arena_prototype *a, char **dstp, const char *fmt, ...) {
	va_list ap;
	int len;

	va_start(ap,fmt);

	len	= arena_util_vasprintf(a,dstp,fmt,ap);

	va_end(ap);

	return len;
} /* arena_util_asprintf() */


char *arena_util_vsprintf(const struct arena_prototype *a, const char *fmt, va_list ap) {
	char *dst;

	if (-1 != arena_util_vasprintf(a,&dst,fmt,ap))
		return dst;
	else
		return 0;
} /* arena_util_vsprintf() */


char *arena_util_sprintf(const struct arena_prototype *a, const char *fmt, ...) {
	va_list ap;
	char *dst;

	va_start(ap,fmt);

	if (-1 == arena_util_vasprintf(a,&dst,fmt,ap))
		dst	= 0;

	va_end(ap);

	return dst;
} /* arena_util_sprintf() */


#if HAVE_STRSEP
int arena_util_split(const struct arena_prototype *a, int *argcp, char ***argvp, char **strp, const char *sep) {
	#define argc	(*argcp)
	#define argv	(*argvp)

	char *str	= *strp;
	char *tmpv[ARENA_UTIL_BUFSZ / sizeof (char *)];
	int nargv;
	char **ap;
	char *nxt;
	char *sub;

	argc	= 0;
	argv	= tmpv;
	nargv	= sizeof tmpv / sizeof *tmpv;

	if (!sep)
		sep	= " \t\r\n";

split:
	if (!(nxt = *strp = arena_util_strdup(a,str)))
		goto fail;

	for (argc = 0, ap = argv;;) {
		if (0 == (sub = strsep(&nxt,sep)))
			break;

		if (argc < nargv)
			*ap	= sub;

		if (*sub != '\0')
			ap++, argc++;
	}

	if (argc >= nargv) {
		if (argv == tmpv) {
			errno	= EFAULT;

			goto fail;
		}

		a->free(a,*strp);
		*strp	= 0;

		if (!(argv = a->malloc(a,sizeof *argv * (argc + 1),0)))
			goto fail;

		nargv	= argc;

		goto split;
	} else if (argv == tmpv) {
		if (!(argv = arena_util_memdup(a,tmpv,sizeof *tmpv * (argc + 1))))
			goto fail;
	}

	argv[argc]	= 0;

	return argc;
fail:
	if (argv != tmpv)
		a->free(a,argv), argv = 0;

	a->free(a,*strp), *strp = 0;

	return -1;

	#undef argc
	#undef argv
} /* arena_util_split() */
#endif


#if 0 && HAVE_STRSEP
int main(int argc, char **argv) {
	const struct arena_prototype *a	= ARENA_STDLIB;
	char *argbuf	= argv[1];

	if (-1 == arena_util_split(a,&argc,&argv,&argbuf,0))
		return EXIT_FAILURE;

	for (int i = 0; i < argc; i++) {
		printf("[%.3d] %s\n",i,argv[i]);
	}

	return 0;
}
#endif
