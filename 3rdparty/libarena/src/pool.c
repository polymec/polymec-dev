/* ==========================================================================
 * libarena/src/pool.c - Custom Memory Allocator Interface
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
#include <stdio.h>
#include <stdlib.h>	/* size_t */

#ifdef _WIN32
typedef unsigned int uint32_t;
#else
#include <inttypes.h>	/* uint32_t */
#endif

#include <stdarg.h>	/* va_list */
#include <stddef.h>	/* offsetof */

#include <string.h>	/* memcpy(3) memmove(3) */

#include <assert.h>

#include "align.h"
#include "proto.h"
#include "pool.h"
#include "rbits.h"
#include "util.h"
#include "queue.h"

#ifndef MIN
#define MIN(a, b)	(((a) < (b))? (a) : (b))
#endif

#ifndef MAX
#define MAX(a, b)	(((a) > (b))? (a) : (b))
#endif


static const struct pool_bucket_options pool_bucket_defaults[] = {
	{ 32,   32 },
	{ 128,  32 },
	{ 512,  32 },
	{ 4096,  8 },
	{ 0,     0 },
};

const struct pool_options pool_defaults = {
	ARENA_SYSTEM_ALIGNMENT,
	pool_bucket_defaults,
};


struct pool_chunk {
	SLIST_ENTRY(pool_chunk) sle;

	unsigned char bytes[];
}; /* struct pool_chunk */


struct pool_block {
	SLIST_ENTRY(pool_block) sle;

	size_t nbytes;

	unsigned char *bytep;

	unsigned char bytes[];
}; /* struct pool_block */


static const struct pool_bucket {
	size_t bucketlen;
	size_t nperblock;

	size_t bytep_offset;
	size_t chunk_size;

	SLIST_HEAD(, pool_chunk) chunks;

	CIRCLEQ_ENTRY(pool_bucket) cqe;
} pool_bucket_initializer;


static const struct pool {
	struct arena_prototype interface;		/* Exportable interface */

	const struct arena_prototype *allocator;	/* Internal "core" allocator */

	size_t alignment;				/* Default alignment */

	SLIST_HEAD(, pool_block) blocks;

	size_t nbuckets;
	CIRCLEQ_HEAD(, pool_bucket) buckets;

	/*
	 * Must provide the range of indexes that pool_hibit() returns,
	 * currently 0 to 31.
	 */
	struct pool_bucket *bucket_index[sizeof(uint32_t) * CHAR_BIT];
} pool_initializer;


/*
 * Returns the highest bit "index" of the passed unsigned integer.
 *
 * Examples:
 *
 * 	1     -> 0  (00000000000000000000000000000001)
 * 	                                            ^ 0th
 *
 * 	65536 -> 15 (00000000000000001111111111111111)
 * 	                             ^ 15th
 *
 * 	65536 -> 16 (00000000000000010000000000000000)
 * 	                            ^ 16th
 */
static inline int pool_hibit(uint32_t n) {
	register int i	= (n & 0xffff0000)? 16 : 0;

	if ((n >>= i) & 0xff00)
		i |= 8, n >>= 8;

	if (n & 0xf0)
		i |= 4, n >>= 4;

	if (n & 0x0c)
		i |= 2, n >>= 2;

	return (i | (n >> 1));
} /* pool_hibit() */


/*
 * Allocate a new memory block and push it onto the block stack. Upon
 * returning block->bytep is suitably aligned and there is sufficient space
 * from block->bytep to provide at least `len' bytes of usable memory.
 */
static struct pool_block *pool_block_push(struct pool *P, size_t len) {
	struct pool_block *blk	= 0;
	/* We must satisfy both data and data structure alignment needs. */
	size_t align	= MAX(ARENA_SYSTEM_ALIGNMENT,P->alignment);
	size_t offset	= ARENA_BOUNDARY_OFFSETOF(offsetof(struct pool_block, bytes),align);

	if (!(blk = P->allocator->malloc(P->allocator,offsetof(struct pool_block, bytes) + offset + len,align)))
		return 0;

	blk->nbytes	= offset + len;
	blk->bytep	= blk->bytes + offset;

	SLIST_INSERT_HEAD(&P->blocks,blk,sle);

	return blk;
} /* pool_block_push() */


static inline int pool_bucket_indexof(struct pool_bucket *b) {
	return pool_hibit(b->bucketlen);
} /* pool_bucket_indexof() */


/*
 * Update the pool bucket lookup index at the slot derived from the
 * provided bucket. Then fill any lower slots.
 */
static void pool_bucket_reindex(struct pool *P, struct pool_bucket *b) {
	int i	= pool_bucket_indexof(b);

	if (!P->bucket_index[i]
	||  b->bucketlen < P->bucket_index[i]->bucketlen) {
		P->bucket_index[i]	= b;

		for (i--; i >= 0 && !P->bucket_index[i]; i--) {
			P->bucket_index[i]	= b;
		}
	}

	return /* void */;
} /* pool_bucket_reindex() */


/*
 * Insert the provided bucket in-order into the pool buckets list.
 */
static void pool_bucket_insert(struct pool *P, struct pool_bucket *b) {
	struct { struct pool_bucket *this; } i;

	/* Try to jump into the list, otherwise begin at the [high] end. */
	if (!(i.this = P->bucket_index[pool_bucket_indexof(b)]))
		i.this	= CIRCLEQ_LAST(&P->buckets);

	while (i.this != CIRCLEQ_END(&P->buckets)) {
		if (b->bucketlen >= i.this->bucketlen)
			break;

		i.this	= CIRCLEQ_PREV(i.this,cqe);
	}

	if (i.this != CIRCLEQ_END(&P->buckets))
		CIRCLEQ_INSERT_AFTER(&P->buckets,i.this,b,cqe);
	else
		CIRCLEQ_INSERT_HEAD(&P->buckets,b,cqe);

	P->nbuckets++;

	return /* void */;
} /* pool_bucket_insert() */


/*
 * Calculate the bytep offset required for maximum guaranteed alignment (the
 * alignment specified in the options) and the actual size to allocate per
 * chunk so chunk arrays are properly aligned.
 */
static void pool_bucket_analyze(struct pool *P, struct pool_bucket *b) {
	size_t align_data	= P->alignment;
	size_t align_struct	= MAX(ARENA_SYSTEM_ALIGNMENT,P->alignment);
	size_t rbits		= rbits_len(b->bucketlen);
	size_t offset		= rbits + ARENA_BOUNDARY_OFFSETOF(rbits + offsetof(struct pool_chunk, bytes),align_data);

#if 0
	printf("align_data: %d\n",(int)align_data);
	printf("align_struct: %d\n",(int)align_struct);
	printf("rbits: %d\n",(int)rbits);
	printf("offsetof(c,bytes): %d\n",(int)offsetof(__typeof__(struct pool_chunk),bytes));
	printf("offset: %d\n",(int)offset);
#endif

	b->bytep_offset	= offset;
	b->chunk_size	= offsetof(struct pool_chunk, bytes) + offset + b->bucketlen;
	b->chunk_size	+= ARENA_BOUNDARY_OFFSETOF(b->chunk_size,align_struct);

	return /* void */;
} /* pool_bucket_analyze() */


static struct pool_bucket *pool_bucket_add(struct pool *P, const struct pool_bucket_options *opts) {
	struct pool_block *blk;
	struct pool_bucket *b;

	if (!(blk = pool_block_push(P,sizeof *b)))
		return 0;

	b		= (void *)blk->bytep;
	blk->bytep	+= sizeof *b;

	b->bucketlen	= opts->bucketlen;
	b->nperblock	= (ARENA_DEBUG)? 1 : opts->nperblock;

	SLIST_INIT(&b->chunks);

	pool_bucket_analyze(P,b);

	pool_bucket_insert(P,b);

	pool_bucket_reindex(P,b);

	return b;
} /* pool_bucket_add() */


static inline size_t pool_power2(size_t i) {
#if defined SIZE_MAX
	i--;
	i |= i >> 1;
	i |= i >> 2;
	i |= i >> 4;
	i |= i >> 8;
	i |= i >> 16;
#if SIZE_MAX != 0xffffffffu
	i |= i >> 32;
#endif
	return ++i;
#else
#error No SIZE_MAX defined
#endif
} /* pool_power2() */


static inline size_t pool_roundup(size_t i) {
	if (i > ~(((size_t)-1) >> 1u))
		return (size_t)-1;
        else
		return pool_power2(i);
} /* pool_roundup() */


/*
 * For a given size find the lowest bucket which can fulfill an allocation
 * of that size. If no bucket is found, try to add a new one.
 */
static struct pool_bucket *pool_bucket_find(struct pool *P, size_t len, int tryhard) {
	struct pool_bucket_options opts;
	struct pool_bucket *b;

	if ((b = P->bucket_index[pool_hibit(len)])) {
		while (b != CIRCLEQ_END(&P->buckets) && b->bucketlen < len)
			b	= CIRCLEQ_NEXT(b,cqe);

		if (b != CIRCLEQ_END(&P->buckets))
			return b;
	}

	if (!tryhard)
		return 0;

	opts.bucketlen	= pool_roundup(len);
	opts.nperblock	= 1;

	return pool_bucket_add(P,&opts);
} /* pool_bucket_find() */


static struct pool_chunk *pool_bucket_grow(struct pool *P, struct pool_bucket *b) {
	struct pool_block *blk;
	struct { struct pool_chunk *pos; struct pool_chunk *end; } c;

	if (!(blk = pool_block_push(P,b->chunk_size * b->nperblock)))
		return 0;

	c.pos	= (void *)(blk->bytep);
	c.end	= (void *)(blk->bytep + (b->nperblock * b->chunk_size));

	while (c.pos < c.end) {
		SLIST_INSERT_HEAD(&b->chunks,c.pos,sle);

		c.pos	= (void *)(((unsigned char *)c.pos) + b->chunk_size);
	}

	blk->bytep	+= b->chunk_size * b->nperblock;

	return SLIST_FIRST(&b->chunks);
} /* pool_bucket_grow() */


static unsigned char *pool_recover(struct pool *P, struct pool_bucket **b, struct pool_chunk **c, unsigned char *q) {
	unsigned char *p;
	size_t bucketlen;

	bucketlen	= rbits_get(q - 1,&p);
	*c		= (void *)(p - offsetof(struct pool_chunk, bytes));

	*b = P->bucket_index[pool_hibit(bucketlen)];
	assert(*b);

	while (*b != CIRCLEQ_END(&P->buckets) && (*b)->bucketlen != bucketlen)
		*b	= CIRCLEQ_NEXT((*b),cqe);

	assert (*b != CIRCLEQ_END(&P->buckets));

	return p;
} /* pool_recover() */


void *pool_get(struct pool *P, size_t len, size_t align) {
	size_t bucketlen	= len;
	struct pool_bucket *b;
	struct pool_chunk *c;
	unsigned char *p;

	if (align == 0)
		align		= P->alignment;
	else if (align > P->alignment)
		bucketlen	+= align - P->alignment;

	if (!(b = pool_bucket_find(P,bucketlen,1)))
		return 0;

	if (!(c = SLIST_FIRST(&b->chunks))
	&&  !(c = pool_bucket_grow(P,b)))
		return 0;

	SLIST_REMOVE_HEAD(&b->chunks,sle);

	p	= &c->bytes[rbits_ptroffset(&c->bytes[0],b->bucketlen,align)];

	(void)rbits_put(&c->bytes[0],p - &c->bytes[0],b->bucketlen,0);

	return p;
} /* pool_get() */


void *pool_realloc(struct pool *P, void *q, size_t dstlen, size_t align) {
	struct pool_chunk *c;
	struct pool_bucket *b, *bx;
	unsigned char *p;
	size_t srcoff, dstoff;

	if (align == 0)
		align	= P->alignment;

	/*
	 * Make things easy on ourselves.
	 */
	if (dstlen == 0) {
		pool_put(P,q);

		return 0;
	}

	if (q == 0)
		return pool_get(P, dstlen, align);

	p	= pool_recover(P,&b,&c,q);

	if (align > P->alignment)
		bx	= pool_bucket_find(P,dstlen + (align - P->alignment),1);
	else
		bx	= pool_bucket_find(P,dstlen,1);

	if (!bx) {
		return 0;
	} else if (bx != b) {	/* Check if (bx->bucketlen > b->bucketlen)? */
		srcoff	= (unsigned char *)q - p;

		if (!(q = pool_get(P,dstlen,align)))
			return 0;

		(void)memcpy(q, p + srcoff, MIN(dstlen, (size_t)(&c->bytes[b->chunk_size - offsetof(struct pool_chunk, bytes)] - (p + srcoff))));

		SLIST_INSERT_HEAD(&b->chunks,c,sle);
	} else {
		srcoff	= (unsigned char *)q - p;
		dstoff	= rbits_ptroffset(p,b->bucketlen,align);

		/* If previous boundary is greater than current, then the
		 * previous is aligned for the new one (given power of 2
		 * arithmetic). If not, we need to shift the contents since
		 * the reverse relationship isn't necessarily true.
		 */
		if (dstoff > srcoff) {
			q	= memmove(p + dstoff, p + srcoff, MIN(dstlen, (size_t)(&c->bytes[b->chunk_size - offsetof(struct pool_chunk, bytes)] - (p + srcoff))));
		} else
			q	= p + srcoff;
	}

	return q;
} /* pool_realloc() */


void pool_put(struct pool *P, void *q) {
	struct pool_chunk *c;
	struct pool_bucket *b;

	if (q == 0)
		return /* void */;

	(void)pool_recover(P,&b,&c,q);

	SLIST_INSERT_HEAD(&b->chunks,c,sle);

	return /* void */;
} /* pool_put() */


static char pool_name[]	= __FILE__;

const char *pool_instanceof(struct pool *P) {
	return &pool_name[0];
} /* pool_instanceof() */


struct pool *pool_import(const struct arena_prototype *ap) {
	return (ap->instanceof(ap) == &pool_name[0])? (struct pool *)ap : 0;
} /* pool_import() */


const char *pool_strerror(struct pool *P) {
	return P->allocator->strerror(P->allocator);
} /* pool_strerror() */


void pool_clearerr(struct pool *P) {
	(P->allocator->clearerr)(P->allocator);

	return /* void */;
} /* pool_clearerr() */


const struct arena_prototype *pool_export(struct pool *P) {
	if (!P->interface.malloc) {
		P->interface.malloc	= (void *(*)(const struct arena_prototype *, size_t, size_t))&pool_get;
		P->interface.realloc	= (void *(*)(const struct arena_prototype *, void *, size_t, size_t))&pool_realloc;
		P->interface.free	= (void (*)(const struct arena_prototype *, void *))&pool_put;
		P->interface.instanceof	= (const char *(*)(const struct arena_prototype *))&pool_instanceof;
		P->interface.strerror	= (const char *(*)(const struct arena_prototype *))&pool_strerror;
		P->interface.clearerr	= (void (*)(const struct arena_prototype *))&pool_clearerr;
	}

	return &P->interface;
} /* pool_export() */


POOL *pool_open(const struct pool_options *opts, const struct arena_prototype *m) {
	struct pool *P	= 0;
	int i;

	if (!opts)
		opts	= &pool_defaults;

	if (!m)
		m	= ARENA_STDLIB;

	if (!(P = m->malloc(m,sizeof *P,0)))
		return 0;

	*P		= pool_initializer;
	P->allocator	= m;
	P->alignment	= opts->alignment;

	SLIST_INIT(&P->blocks);
	CIRCLEQ_INIT(&P->buckets);

	for (i = 0; opts->buckets[i].bucketlen > 0; i++) {
		if (!pool_bucket_add(P,&opts->buckets[i]))
			goto fail;
	}

	return P;
fail:
	pool_close(P);

	return 0;
} /* pool_open() */


void pool_close(POOL *P) {
	struct pool_block *b;

	/*
	 * Release everything in reverse order. Block list is a LIFO.
	 */
	if (P) {
		while ((b = SLIST_FIRST(&P->blocks))) {
			SLIST_REMOVE_HEAD(&P->blocks,sle);

			P->allocator->free(P->allocator,b);
		}

		P->allocator->free(P->allocator,P);
	}

	return /* void */;
} /* pool_close() */


char *pool_strdup(struct pool *P, const char *src) {
	return arena_util_strdup(pool_export(P),src);
} /* pool_strdup() */


char *pool_strndup(struct pool *P, const char *src, size_t n) {
	return arena_util_strndup(pool_export(P),src,n);
} /* pool_strndup() */


void *pool_memdup(struct pool *P, const void *p, size_t n) {
	return arena_util_memdup(pool_export(P),p,n);
} /* pool_memdup() */


int pool_vasprintf(struct pool *P, char **dstp, const char *fmt, va_list ap) {
	return arena_util_vasprintf(pool_export(P),dstp,fmt,ap);
} /* pool_vasprintf() */


int pool_asprintf(struct pool *P, char **dstp, const char *fmt, ...) {
	va_list ap;
	int n;

	va_start(ap,fmt);

	n	= arena_util_vasprintf(pool_export(P),dstp,fmt,ap);

	va_end(ap);

	return n;
} /* pool_asprintf() */


char *pool_vsprintf(struct pool *P, const char *fmt, va_list ap) {
	return arena_util_vsprintf(pool_export(P),fmt,ap);
} /* pool_vsprintf() */


char *pool_sprintf(struct pool *P, const char *fmt, ...) {
	va_list ap;
	char *s;

	va_start(ap,fmt);

	s	= arena_util_vsprintf(pool_export(P),fmt,ap);

	va_end(ap);

	return s;
} /* pool_sprintf() */


#if POOL_MAIN

#include <stdio.h>
#include <err.h>
#include <string.h>

int main(int argc, char *argv[]) {
	struct pool_options opts = pool_defaults;
	POOL *p;
	struct pool_bucket *b;
	int i;

	/*return printf("%d\n",rbits_ptroffset((unsigned char *)16,23,16));*/

	opts.alignment = (argc > 1)? atoi(argv[1]) : ARENA_SYSTEM_ALIGNMENT;

	if (!(p = pool_open(&opts,0)))
		err(1,"pool_open");

	CIRCLEQ_FOREACH(b,&p->buckets,cqe) {
		printf("bucket: %d\n",(int)b->bucketlen);
		printf("bytep_offset: %d\n",(int)b->bytep_offset);
		printf("chunk_size: %d\n",(int)b->chunk_size);
		printf("\n");
	}

	for (i = 0; i < 1024; i++) {
		void *ptr	= pool_get(p,12,16);

		strcpy(ptr,"abcdefghijk");

		pool_put(p,ptr);
	}

	pool_close(p);

	return 0;
}

#endif /* POOL_MAIN */
