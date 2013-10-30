/* ==========================================================================
 * libarena/src/proto.c - Custom Memory Allocator Interface
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
#include <stdlib.h>	/* size_t malloc(3) realloc(3) free(3) */
#include <stddef.h>	/* offsetof() */

#include <string.h>	/* strerror(3) memmove(3) */

#include <errno.h>	/* ENOMEM errno */

#include "align.h"
#include "proto.h"	/* struct arena_prototype */
#include "rbits.h"

#ifndef MIN
#define MIN(a, b)	(((a) < (b))? (a) : (b))
#endif

#ifndef MAX
#define MAX(a, b)	(((a) > (b))? (a) : (b))
#endif


/*
 * ARENA_SYSTEM_ALIGNMENT is the lowest common alignment denominator for all
 * datatypes. However, the system malloc implementation may align even
 * stricter. For instance, OpenBSD's malloc aligns on 16 byte boundaries on
 * all implementations, even though most only strictly require an 8 byte
 * boundary.
 */
#ifndef ARENA_SYSTEM_MALLOC_ALIGNS
#define ARENA_SYSTEM_MALLOC_ALIGNS	ARENA_SYSTEM_ALIGNMENT
#endif


const char *ARENA_NOERROR	= "Success";


static int sys_errno;


static void *sys_malloc(const struct arena_prototype *sys, size_t size, size_t align) {
	unsigned char *p;
	size_t ptroff;

	if (align == 0)
		align	= ARENA_SYSTEM_ALIGNMENT;

	ptroff	= rbits_ptroffset((unsigned char *)ARENA_SYSTEM_MALLOC_ALIGNS,size,align);

	if (!(p = malloc(ptroff + size)))
		sys_errno = errno;

	(void)rbits_put(p,ptroff,size,0);

	return p + ptroff;
} /* sys_malloc() */


static void *sys_realloc(const struct arena_prototype *sys, void *q, size_t dstlen, size_t align) {
	unsigned char *p;
	size_t srcoff, dstoff, srclen;

	if (align == 0)
		align	= ARENA_SYSTEM_ALIGNMENT;

	if (q) {
		srclen	= rbits_get((unsigned char *)q - 1,&p);
		srcoff	= (unsigned char *)q - p;
		p	= (unsigned char *)q - srcoff;
	} else {
		srclen	= 0;
		srcoff	= 0;
		p	= q;
	}

	if (dstlen == 0) {
		free(p);

		return 0;
	}

	dstoff	= MAX(srcoff,rbits_ptroffset((unsigned char *)ARENA_SYSTEM_MALLOC_ALIGNS,dstlen,align));

	if (!(p = realloc(p,dstoff + dstlen))) {
		sys_errno	= errno;

		return 0;
	}

	/* Do we need to shift the contents to write out the new size? */
	if (dstoff > srcoff) {
		(void)memmove(p + dstoff,p + srcoff,MIN(srclen,dstlen));
	}

	(void)rbits_put(p,dstoff,dstlen,0);

	return p + dstoff;
} /* sys_realloc() */


static void sys_free(const struct arena_prototype *sys, void *q) {
	unsigned char *p;

	if (q == 0)
		return /* void */;

	(void)rbits_get((unsigned char *)q - 1,&p);

	free(p);

	return /* void */;
} /* sys_free() */


static const char sys_name[] = "stdlib + alignment";

static const char *sys_instanceof(const struct arena_prototype *sys) {
	return &sys_name[0];
} /* sys_instanceof() */


static const char *sys_strerror(const struct arena_prototype *sys) {
	if (sys_errno)
		return strerror(sys_errno);

	return ARENA_NOERROR;
} /* sys_strerror() */


static void sys_clearerr(const struct arena_prototype *sys) {
	sys_errno	= 0;
} /* sys_clearerr() */


static const struct arena_prototype arena_stdlib_aligned = {
	/* .malloc	= */ sys_malloc,
	/* .realloc	= */ sys_realloc,
	/* .free	= */ sys_free,
	/* .instanceof	= */ sys_instanceof,
	/* .strerror	= */ sys_strerror,
	/* .clearerr	= */ sys_clearerr,
};


const struct arena_prototype *ARENA_STDLIB_ALIGNED	= &arena_stdlib_aligned;



static int null_errno;


static void *null_malloc(const struct arena_prototype *null, size_t size, size_t align) {
	void *p;

	if (!(p = malloc(size)))
		return null_errno = errno, (void *)0;

#if 0 /* If we can't do it for null_realloc, may as well not do it here. */
	if (0 != ((uintptr_t)p % align))
		return free(p), null_errno = errno, (void *)0;
#endif

	return p;
} /* null_malloc() */


static void *null_realloc(const struct arena_prototype *null, void *p, size_t size, size_t align) {
	if (!(p = realloc(p, size)))
		return null_errno = errno, (void *)0;

	return p;
} /* null_realloc() */


static void null_free(const struct arena_prototype *null, void *p) {
	free(p);

	return /* void */;
} /* null_free() */


static const char null_name[] = "null";

static const char *null_instanceof(const struct arena_prototype *null) {
	return &null_name[0];
} /* null_instanceof() */


static const char *null_strerror(const struct arena_prototype *null) {
	if (null_errno)
		return strerror(null_errno);

	return ARENA_NOERROR;
} /* null_strerror() */


static void null_clearerr(const struct arena_prototype *null) {
	null_errno	= 0;
} /* null_clearerr() */


static const struct arena_prototype arena_stdlib = {
	/* .malloc	= */ null_malloc,
	/* .realloc	= */ null_realloc,
	/* .free	= */ null_free,
	/* .instanceof	= */ null_instanceof,
	/* .strerror	= */ null_strerror,
	/* .clearerr	= */ null_clearerr,
};

const struct arena_prototype *ARENA_STDLIB	= &arena_stdlib;



int arena_debug(void) {
	extern char *getenv(const char *);
	static int debug	= 0;

	if (!debug)
		debug	= (getenv("ARENA_DEBUG"))? 1 : -1;

	return (debug > 0)? 1 : 0;
} /* arena_debug() */
