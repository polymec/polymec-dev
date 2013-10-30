/* ==========================================================================
 * libarena/src/arena.c - Custom Memory Allocator Interface
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
#include <stdlib.h>	/* size_t getenv(3) */
#include <stddef.h>	/* offsetof */
#include <stdarg.h>	/* va_list */

#ifdef _WIN32
#include <stddef.h>		/* ptrdiff_t */
#else
#include <inttypes.h>	/* ptrdiff_t */
#endif

#include <limits.h>	/* CHAR_BIT */

#include <string.h>	/* memcpy(3) memmove(3) */

#include <assert.h>	/* assert */

#include "align.h"
#include "proto.h"
#include "arena.h"
#include "rbits.h"
#include "util.h"
#include "queue.h"	/* SLIST_HEAD SLIST_ENTRY SLIST_INIT SLIST_FIRST */

#ifndef MIN
#define MIN(a, b)	(((a) < (b))? (a) : (b))
#endif

#ifndef MAX
#define MAX(a, b)	(((a) > (b))? (a) : (b))
#endif

#ifndef roundup
#define roundup(x, y)	((((x) + ((y) - 1)) / (y)) * (y))
#endif

static const struct arena_block {
	size_t size;

	struct {
		unsigned char *next;
	} pos;

	SLIST_ENTRY(arena_block) sle;

	unsigned char bytes[];
} arena_block_initializer;


const struct arena_options arena_defaults = {
	ARENA_SYSTEM_ALIGNMENT,
	1 << 15,	/* 32768 */
};



static const struct arena {
	struct arena_prototype interface;		/* Exportable interface */

	const struct arena_prototype *allocator;	/* Internal "core" allocator */

	SLIST_HEAD(,arena_block) blocks;
	unsigned int nblocks;

	struct arena_options options;
} arena_initializer;


static struct arena_block *arena_block_malloc(struct arena *a, size_t len, size_t align) {
	struct arena_block *b;
	size_t size;

	if (align == 0)
		align	= a->options.alignment;

	if (!ARENA_DEBUG)
		size	= MAX(a->options.blocklen,offsetof(struct arena_block, bytes) + len + align - 1 + RBITS_MAXLEN);
	else
		size	= MAX(a->options.blocklen,offsetof(struct arena_block, bytes) + len + align - 1 + rbits_len(len));

	if (!(b = a->allocator->malloc(a->allocator,size,0)))
		return 0;

	*b		= arena_block_initializer;
	b->size		= size - offsetof(struct arena_block, bytes);
	b->pos.next	= b->bytes;

	return b;
} /* arena_block_malloc() */


static int arena_block_regionof(struct arena_block *b, void *p) {
	unsigned char *x	= p;

	return (x >= b->bytes && x < &b->bytes[b->size]);
} /* arena_block_regionof() */


struct arena *arena_open(const struct arena_options *opts, const struct arena_prototype *m) {
	struct arena tmp, *a;
	struct arena_block *b;

	if (!opts)
		opts	= &arena_defaults;

	if (!m)
		m	= ARENA_STDLIB;

	tmp		= arena_initializer;
	tmp.allocator	= m;
	tmp.options	= *opts;

	if (ARENA_DEBUG)
		tmp.options.blocklen	= 0;

	if (!(b = arena_block_malloc(&tmp,sizeof *a,0)))
		return 0;

	a		= (void *)(b->pos.next + rbits_ptroffset(b->pos.next,sizeof *a,ARENA_SYSTEM_ALIGNMENT));

	rbits_put(b->pos.next, (unsigned char *)a - b->pos.next, sizeof *a, 0);

	b->pos.next	= (void *)(a + 1);

	*a		= arena_initializer;

	a->allocator	= m;

	SLIST_INIT(&a->blocks)

	SLIST_INSERT_HEAD(&a->blocks,b,sle);
	a->nblocks++;

	a->options	= *opts;

	if (ARENA_DEBUG)
		a->options.blocklen	= 0;

	return a;
} /* arena_open() */


void arena_close(struct arena *a) {
	struct arena_block *b, *nxt;

	if (a == 0)
		return /* void */;

	for (b = SLIST_FIRST(&a->blocks); b; b = nxt) {
		nxt	= SLIST_NEXT(b,sle);

		a->allocator->free(a->allocator,b);
	}

	return /* void */;
} /* arena_close() */


void *arena_malloc(struct arena *a, size_t size, size_t align) {
	struct arena_block *b	= SLIST_FIRST(&a->blocks);
	unsigned char *p;

	if (size == 0)
		return 0;

	if (align == 0)
		align	= a->options.alignment;

/*	printf("p getoffset(%p,%d,%d): %d\n",b->pos.next,size,align,getoffset(b->pos.next,size,align));*/

	p	= b->pos.next + rbits_ptroffset(b->pos.next,size,align);

	if (!(p + size <= &b->bytes[b->size])) {
		size_t want;

		/*
		 * FIXME: This is sort of a hack. Think of a better way to
		 * implement this behavior, if not better behavior.
		 *
		 * If somebody is sequentially and repeatedly growing a
		 * buffer beyond the default block length, try to keep
		 * apace.
		 */
		if (size > a->options.blocklen)
			want	= roundup(2 * size, MAX(1, a->options.blocklen));
		else
			want	= size;

		if (!(b = arena_block_malloc(a,want,align)))
			return 0;

		SLIST_INSERT_HEAD(&a->blocks,b,sle);
		a->nblocks++;

		p	= b->pos.next + rbits_ptroffset(b->pos.next,size,align);
	}

	rbits_put(b->pos.next,p - b->pos.next,size,0);

	b->pos.next	= p + size;

	return p;
} /* arena_malloc() */


void *arena_realloc(struct arena *a, void *p_, size_t n, size_t align) {
	struct { size_t old; size_t new; } len;
	unsigned char *p	= p_;
	unsigned char *hp	= 0;	/* Head of pointer */
	unsigned char *tp;		/* Tail of pointer */
	struct arena_block *top	= SLIST_FIRST(&a->blocks);
	size_t offset;

	if (align == 0)
		align	= a->options.alignment;

	if (p == 0)
		return arena_malloc(a,n,align);
	else if (n == 0)
		return arena_free(a,p), (void *)0;

	len.new	= n;

	len.old = rbits_get(p - 1,&hp);
	assert((len.old > 0) && (hp != 0));

	tp	= p + len.old;

	/*
	 * NOTE: Realignment could make a shorter reallocation consume more
	 * space than the original allocation. The code follows from that
	 * observation.
	 */
	offset	= rbits_ptroffset(hp,len.new,align);

	if (hp + offset + len.new <= tp) {
		ptrdiff_t diff	= (p - hp) - offset;

		if (diff)
			p	= memmove(hp + offset,p,MIN(len.new,len.old));

		/* If the tail of the original allocation was the top of our
		 * stack, drop the stack pointer so the space can be reused.
		 * Otherwise, maintain the stored region size so a reverse
		 * order free will work as advertised (rather than losing
		 * track of previous allocations and leaking).
		 */
		if (top->pos.next == tp) {
			top->pos.next	= p + len.new;

			rbits_put(hp,p - hp,len.new,0);
		} else {
			/* Store the full size of this region, not just the
			 * requested.
			 */
			rbits_put(hp,p - hp,tp - p,0);
		}
	} else {
		/* Are we at the top of the stack and the requested
		 * allocation still fits? If so, shift the stack pointer and
		 * memmove the data. Otherwise, allocate a new region and
		 * copy the data, possibly leaking the old region.
		 */
		if (top->pos.next == tp && hp + offset + len.new <= &top->bytes[top->size]) {
			p	= memmove(hp + offset,p,MIN(len.new,len.old));

			top->pos.next	= p + len.new;

			rbits_put(hp,p - hp,len.new,0);
		} else if ((p = arena_malloc(a,len.new,align))) {
			/* If we were top of stack shift the pointer down. 
			 * This allows reallocations at the top of the stack
			 * to overflow a block, but still provide leak free
			 * reverse deallocation later on.
			 */
			 if (top->pos.next == tp)
			 	top->pos.next	= hp;

			(void)memcpy(p,p_,MIN(len.new,len.old));
		}
	}

	return p;
} /* arena_realloc() */


void arena_free(struct arena *a, void *p_) {
	unsigned char *hp, *tp, *p = p_;
	struct arena_block *top;
	size_t len;

	if (!p)
		return /* void */;

	top = SLIST_FIRST(&a->blocks);
	assert(top);

	len = rbits_get(p - 1,&hp);
	assert((len > 0) && (hp != 0));

	tp	= p + len;

	if (top->pos.next != tp)
		return /* void */;

	/* If the block is empty, free it */
	if (top->bytes == (top->pos.next = hp)) {
		SLIST_REMOVE_HEAD(&a->blocks,sle);
		a->nblocks--;

		a->allocator->free(a->allocator,top);
	}

	return /* void */;
} /* arena_free() */


void arena_mark(struct arena *a, void **p) {
	struct arena_block *top;

	top = SLIST_FIRST(&a->blocks);
	assert(top);

	*p	= top->pos.next;

	return /* void */;
} /* arena_mark() */


void arena_reset(struct arena *a, void *p_) {
	unsigned char *p	= p_;
	struct arena_block *b;

	while ((b = SLIST_FIRST(&a->blocks))) {
		if (arena_block_regionof(b,p)
		||  arena_block_regionof(b,p - 1)) {
			b->pos.next	= p;

			return /* void */;
		}

		/* Don't free ourselves! */
		assert(a->nblocks > 1);

		SLIST_REMOVE_HEAD(&a->blocks,sle);
		a->nblocks--;

		a->allocator->free(a->allocator,b);
	}

	assert(0 && "Bad arena reset!");

	return /* void */;
} /* arena_mark() */


static char arena_name[]	= __FILE__;

const char *arena_instanceof(struct arena *a) {
	return &arena_name[0];
} /* arena_instanceof() */


struct arena *arena_import(const struct arena_prototype *ap) {
	return (ap->instanceof(ap) == &arena_name[0])? (struct arena *)ap : 0;
} /* arena_import() */


size_t arena_lengthof(struct arena *a, void *p) {
	unsigned char *hp;
	size_t len;

	len	= rbits_get((unsigned char *)p - 1,&hp);

/*	printf("len=%d; p - hp=%d; p=%p\n",len,(int)((unsigned char *)p - hp),p);*/

	return len;
} /* arena_lengthof() */


int arena_regionof(struct arena *a, void *p) {
	struct arena_block *b;

	SLIST_FOREACH(b,&a->blocks,sle) {
		if (arena_block_regionof(b,p))
			return 1;
	}

	return 0;
} /* arena_regionof() */


const char *arena_strerror(struct arena *a) {
	return a->allocator->strerror(a->allocator);
} /* arena_strerror() */


void arena_clearerr(struct arena *a) {
	(a->allocator->clearerr)(a->allocator);

	return /* void */;
} /* arena_clearerr() */


const struct arena_prototype *arena_export(struct arena *a) {

	if (!a->interface.malloc) {
		a->interface.malloc	= (void *(*)(const struct arena_prototype *, size_t, size_t))&arena_malloc;
		a->interface.realloc	= (void *(*)(const struct arena_prototype *, void *, size_t, size_t))&arena_realloc;
		a->interface.free	= (void (*)(const struct arena_prototype *, void *))&arena_free;
		a->interface.instanceof	= (const char *(*)(const struct arena_prototype *))&arena_instanceof;
		a->interface.strerror	= (const char *(*)(const struct arena_prototype *))&arena_strerror;
		a->interface.clearerr	= (void (*)(const struct arena_prototype *))&arena_clearerr;
	}

	return &a->interface;
} /* arena_export() */


char *arena_strdup(struct arena *P, const char *src) {
	return arena_util_strdup(arena_export(P),src);
} /* arena_strdup() */


char *arena_strndup(struct arena *P, const char *src, size_t n) {
	return arena_util_strndup(arena_export(P),src,n);
} /* arena_strndup() */


void *arena_memdup(struct arena *P, const void *p, size_t n) {
	return arena_util_memdup(arena_export(P),p,n);
} /* arena_memdup() */


int arena_vasprintf(struct arena *P, char **dstp, const char *fmt, va_list ap) {
	return arena_util_vasprintf(arena_export(P),dstp,fmt,ap);
} /* arena_vasprintf() */


int arena_asprintf(struct arena *P, char **dstp, const char *fmt, ...) {
	va_list ap;
	int n;

	va_start(ap,fmt);

	n	= arena_util_vasprintf(arena_export(P),dstp,fmt,ap);

	va_end(ap);

	return n;
} /* arena_asprintf() */


char *arena_vsprintf(struct arena *P, const char *fmt, va_list ap) {
	return arena_util_vsprintf(arena_export(P),fmt,ap);
} /* arena_vsprintf() */


char *arena_sprintf(struct arena *P, const char *fmt, ...) {
	va_list ap;
	char *s;

	va_start(ap,fmt);

	s	= arena_util_vsprintf(arena_export(P),fmt,ap);

	va_end(ap);

	return s;
} /* arena_sprintf() */


#if ARENA_MAIN


#include <stdio.h>


void test0(int argc, char *argv[]) {
	ARENA *a;
	struct arena_options opts = arena_defaults;
	char *list[256];
	int i, n;

	opts.blocklen = 64;

	if (!(a = arena_open(&opts,ARENA_STDLIB)))
		err(EXIT_FAILURE,"arena_open");

	for (n = 0; n < 8; n++) {
	for (i = 0; i < sizeof list / sizeof *list; i++) {
		/*if (!(list[i] = arena_malloc(a,2,1)))
			err(EXIT_FAILURE,"arena_malloc");

		list[i][0]	= '0' + i;
		list[i][1]	= '\0';*/
		if (!(list[i] = arena_util_sprintf(arena_export(a),"%d",i)))
			err(EXIT_FAILURE,"arena_util_sprintf");
	}

	for (i = 0; i < sizeof list / sizeof *list; i++) {
		printf("list[%d] = %s\n",i,list[i]);
		printf("regionof: %s\n",(arena_regionof(a,list[i]))? "yes" : "no");
		printf("lengthof: %d\n",(int)arena_lengthof(a,list[i]));
	}

	printf("nblocks = %d\n",(int)a->nblocks);

	for (i = (sizeof list / sizeof *list) - 1; i >= 0; i--) {
		arena_free(a,list[i]);
	}

	printf("nblocks = %d\n",(int)a->nblocks);
	}

	arena_close(a);
} /* test0() */


int main(int argc, char *argv[]) {
	int n;
	unsigned char len[RBITS_MAXLEN];
	size_t size;
	unsigned char *last;
	unsigned char *p;

#if 0
	printf("RBITS_MAXLEN = %d\n",RBITS_MAXLEN);

	size	= strtoul(argv[1],0,0);

	last	= rbits_put(len,sizeof len,size,0);

	printf("last - len: %d\n",(int)(last - len));

	/*printf("n = %d\n",n);*/

	p	= &len[sizeof len - 1];
	size	= rbits_get(p,&p);

	printf("size = %lu\n",(unsigned long)size);
	printf("p - len: %d\n",(int)(p - len));
#endif
	test0(argc,argv);

	return 0;
}

#endif /* ARENA_MAIN */
