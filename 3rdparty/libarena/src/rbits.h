/* ==========================================================================
 * libarena/rbits.h - Custom Memory Allocator Interface
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
 * --------------------------------------------------------------------------
 * Reverse Variable-length Bit-strings
 *
 * Marshal and unmarshal an rbitsint_t integer into a reverse bit-string;
 * that is, written out right-to-left (in cardinality of memory space).
 * ==========================================================================
 */
#ifndef ARENA_RBITS_H
#define ARENA_RBITS_H

#ifdef _WIN32
#define inline
#endif

#include <stddef.h>	/* size_t */

#ifdef _WIN32
typedef unsigned long uintptr_t;
#else
#include <inttypes.h>	/* uintptr_t */
#endif

#include <limits.h>	/* CHAR_BIT */


#ifndef rbitsint_t
#define rbitsint_t	size_t
#endif


/* Maximum space needed to store an rbitsint_t using CHAR_BIT - 1 bits */
#define RBITS_MAXLEN	(sizeof (rbitsint_t) + (((sizeof (rbitsint_t) * CHAR_BIT) - (sizeof (rbitsint_t) * (CHAR_BIT - 1))) / CHAR_BIT) + 1)


/*
 * Store the bit value representation of an rbitsint_t integer across buf,
 * starting from the end, preserving the highest bit of each unsigned char
 * in buf for use as a delimiter. Returns the last part of buf (lowest
 * unsigned char *) that was written to. If compact is set this will be the
 * part that holds the highest order bit of size (equal or greater than
 * buf), otherwise buf.
 */
static inline unsigned char *rbits_put(unsigned char *buf, size_t buflen, rbitsint_t i, int compact) {
	unsigned char *c	= &buf[buflen];
	unsigned char *last	= c;

	/* Iterate backwards, storing the size in all but the highest bit of
	 * each byte. The highest bit serves as a marker telling us when to
	 * stop.
	 */
	do {
		c--;

		/* Assign all but the highest bit, which is preserved. */
		if ((*c = i & ~(1U << (CHAR_BIT - 1))))
			last	= c;

		i >>= CHAR_BIT - 1;
	} while (c > buf);

	if (!compact)
		last	= buf;

	/* Tag our terminal byte */
	*last	|= 1U << (CHAR_BIT - 1);

	return last;
} /* rbits_put() */


/*
 * Return the buffer size required to hold an rbitsint_t value bit-string
 * representation.
 */
static inline size_t rbits_len(rbitsint_t i) {
	unsigned char buf[RBITS_MAXLEN];
	unsigned char *pos	= rbits_put(buf,sizeof buf,i,1);

	return &buf[sizeof buf] - pos;
} /* rbits_len() */


/*
 * Return the offset from p required to 1) store the requested size and 2)
 * align to the desired alignment.
 */
static inline size_t rbits_ptroffset(unsigned char *p, size_t size, size_t align) {
	unsigned char lenbuf[RBITS_MAXLEN];
	unsigned char *lenend	= &lenbuf[sizeof lenbuf - 1];
	unsigned char *lenpos;
	uintptr_t ptrpos	= (uintptr_t)p;

	lenpos	= rbits_put(lenbuf,sizeof lenbuf,size,1);

	ptrpos	+= (lenend - lenpos) + 1;
	ptrpos	+= ARENA_BOUNDARY_OFFSETOF(ptrpos,align); /* Needs power of 2. */

	return ptrpos - (uintptr_t)p;
} /* rbits_ptroffset() */


/*
 * Beginning from *p, work backwards reconstructing the value of an
 * rbitsint_t integer. Stop when the highest order bit of *p is set, which
 * should have been previously preserved as a marker. Return the
 * reconstructed value, setting *end to the last position used of p.
 */
static inline rbitsint_t rbits_get(unsigned char *p, unsigned char **end) {
	rbitsint_t i	= 0;
	int n		= 0;

	do {
		i	|= (*p & ~(1 << (CHAR_BIT - 1))) << (n++ * (CHAR_BIT - 1));
	} while (!(*(p--) & (1 << (CHAR_BIT - 1))));

	*end	= p + 1;

	return i;
} /* rbits_get() */


#endif /* ARENA_RBITS_H */
