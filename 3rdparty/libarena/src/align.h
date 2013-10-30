/* ==========================================================================
 * arena/align.h - Custom Memory Allocator Interface
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
#ifndef ARENA_ALIGN_H
#define ARENA_ALIGN_H

#include <stddef.h>	/* offsetof() */


/*
 * Taken from glibc 2.3.5 ptmalloc2 implementation. Seems reasonable.
 * Possible alternative: (sizeof (union { void *p; double d; })).
 */
#ifndef ARENA_SYSTEM_ALIGNMENT
#define ARENA_SYSTEM_ALIGNMENT	(2 * sizeof (size_t))
#endif


/*
 * Calculates the adjustment needed to push `p' to boundary `align'. 
 *
 * NOTE: `align' MUST BE a power of 2.
 */
#ifndef ARENA_BOUNDARY_OFFSETOF
#define ARENA_BOUNDARY_OFFSETOF(p,align) \
	(((align) - ((uintptr_t)(p) % (align))) & ~(align))
#endif


/*
 * Calculates the adjustment needed to push `p' to boundary `align'.  This
 * version does not constrain the arguments, specifically it does not
 * require `align' to be a power of 2.
 */
#ifndef ARENA_XBOUNDARY_OFFSETOF
#define ARENA_XBOUNDARY_OFFSETOF(p,align) \
	((((align) - ((uintptr_t)(p) % (align))) != (align))? ((align) - ((uintptr_t)(p) % (align))) : (align))
#endif


/*
 * Determine a safe, but not necessarily smallest, boundary that type `t'
 * can be aligned on.
 */
#ifndef arena_alignof
#define arena_alignof(t)	offsetof(struct { char c; t x; }, x)
#endif


/*
 * Is `i' a power of 2?
 */
#ifndef arena_powerof2
#define arena_powerof2(i)	((((i) - 1) & (i)) == 0)
#endif

#endif /* ARENA_ALIGN_H */
