/* ==========================================================================
 * libarena/proto.h - Custom Memory Allocator Interface
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
#ifndef ARENA_PROTO_H
#define ARENA_PROTO_H

/*
 * Generic memory allocator interface.
 * .malloc
 * 	USAGE
 * 		malloc(<context>,<size>,<alignment>)
 *
 *	DESCRIPION
 * 		Allocates at least <size> bytes of memory, starting from
 * 		<alignment> boundary.
 *
 * 	RETURN VALUE
 * 		Returns a pointer to the allocated memory on success, or a
 * 		NULL pointer on failure.
 *
 * .realloc
 * 	USAGE
 * 		realloc(<context>,<pointer>,<size>,<alignment>)
 *
 * 	DESCRIPTION
 * 		Allocates a region of memory consisting of at least <size>
 * 		bytes with the same values as the previously allocated
 * 		region (up to the length of the lesser of the new and old
 * 		sizes) and aligned on an <alignment> byte boundary.
 *
 * 	RETURN VALUE
 * 		Returns a pointer to the allocated region on success, or a
 * 		NULL pointer on failure. On success use of <pointer> is
 * 		undefined, on failure the region named by <pointer> remains
 * 		unchanged.
 *
 * .free
 * 	USAGE
 * 		free(<context>,<pointer>)
 *
 * 	DESCRIPTION
 * 		Attempt to release all resources associated with the memory
 * 		region described by <pointer>.
 *
 * 	RETURN VALUE
 * 		None.
 *
 * .instanceof
 * 	USAGE
 * 		instanceof(<context>)
 *
 * 	DESCRIPTION
 * 		Returns a C string pointer which can be used by an allocator
 * 		specific interface or in an allocator specific manner to
 * 		uniquely determine the constructor of the arena instance. 
 * 		The pointer value, and not the contents of the C string,
 * 		identifies the instance. The C string can be used in
 * 		debugging output, etc., to describe the allocation strategy
 * 		and/or mechanism of the instance.
 *
 * 	RETURN VALUE
 * 		Returns a non-NULL, distinct pointer. However, there is no
 * 		guarantee any known allocator will be positively
 * 		identifiable.
 *
 * .strerror
 * 	USAGE
 * 		strerror(<context>)
 *
 * 	DESCRIPTION
 * 		Provides a descriptive string describing the last known
 * 		error encountered.
 *
 * 	RETURN VALUE
 * 		Immutable string describing last error, or ARENA_NOERROR if
 * 		no known error has occured.
 *
 * .clearerr
 * 	USAGE
 * 		clearerr(<context>)
 *
 * 	DESCRIPTION
 * 		Erase all previous error conditions.
 *
 * 	RETURN VALUE
 * 		None.
 */
typedef struct arena_prototype {
	void *(*malloc)(const struct arena_prototype *, size_t, size_t);

	void *(*realloc)(const struct arena_prototype *, void *, size_t, size_t);

	void (*free)(const struct arena_prototype *, void *);

	const char *(*instanceof)(const struct arena_prototype *);

	const char *(*strerror)(const struct arena_prototype *);

	void (*clearerr)(const struct arena_prototype *);
} arena_t;


/*
 * Special string pointer used to return "no error" from the strerror()
 * method of a prototyped arena interface.
 */
extern const char *ARENA_NOERROR;


/*
 * Pre-built arena interfaces built directly atop malloc(3), realloc(3) and
 * free(3).
 */
extern const struct arena_prototype *ARENA_STDLIB;
extern const struct arena_prototype *ARENA_STDLIB_ALIGNED;


/*
 * Debugging interface, the macro and function return non-zero if debugging
 * is enabled, zero otherwise.
 */
#define ARENA_DEBUG	arena_debug()

int arena_debug(void);


/*
 * Bundle all the important stuff with <arena/proto.h> for convenience.
 */
#ifndef LIBARENA_SOURCE
#include <arena/align.h>
#include <arena/arena.h>
#endif


#endif /* ARENA_PROTO_H */
