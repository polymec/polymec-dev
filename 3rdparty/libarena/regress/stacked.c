#include <stdio.h>	/* vsnprintf(3) fprintf(3) */
#include <stdlib.h>	/* EXIT_FAILURE EXIT_SUCCESS srandom(3) random(3) exit(3) */
#include <stddef.h>	/* offsetof */

#include <string.h>	/* strerror(3) strlen(3) */

#include <errno.h>	/* errno */

#include <time.h>	/* time(2) */

#include <zlib.h>

#include "pool.h"
#include "arena.h"
#include "proto.h"
#include "queue.h"


struct data {
	unsigned long order;

	const struct arena_prototype *ap;

	unsigned long checksum;

	LIST_ENTRY(data) le;

	size_t nbytes;

	unsigned char bytes[1];
};


/*
 * Fake err(3) and errx(3) from err.h since Solaris doesn't provide it.
 */
static void err(int retcode, const char *fmt, ...) {
	int sys_errno	= errno;
	char buf[256];
	va_list ap;

	va_start(ap, fmt);

	(void)vsnprintf(buf, sizeof buf, fmt, ap);

	va_end(ap);

	(void)fprintf(stderr, "%s: %s\n", buf, strerror(sys_errno));

	exit(retcode);
} /* err() */


static void errx(int retcode, const char *fmt, ...) {
	int sys_errno	= errno;
	char buf[256];
	va_list ap;

	va_start(ap, fmt);

	(void)vsnprintf(buf, sizeof buf, fmt, ap);

	va_end(ap);

	(void)fprintf(stderr, "%s\n", buf);

	exit(retcode);
} /* errx() */


int main(int argc, char *argv[]) {
	LIST_HEAD(,data) list	= LIST_HEAD_INITIALIZER(&list);
	const struct arena_prototype *ap[6];
	ARENA *a0, *a1, *a2;
	POOL *p0, *p1, *p2;
	FILE *fp;
	struct data *d, *dn;
	unsigned long i, j, n;

	if (!(p0 = pool_open(&pool_defaults,ARENA_STDLIB)))
		err(EXIT_FAILURE,"pool_open(&pool_defaults,0)");

	if (!(a0 = arena_open(&arena_defaults,0)))
		err(EXIT_FAILURE,"arena_open(&arena_defaults,0)");

	if (!(p1 = pool_open(&pool_defaults,arena_export(a0))))
		err(EXIT_FAILURE,"pool_open(&pool_defaults,%p)",(void *)arena_export(a0));

	if (!(a1 = arena_open(&arena_defaults,pool_export(p0))))
		err(EXIT_FAILURE,"arena_open(&arena_defaults,%p)",(void *)pool_export(p0));

	if (!(p2 = pool_open(&pool_defaults,arena_export(a1))))
		err(EXIT_FAILURE,"pool_open(&pool_defaults,%p)",(void *)arena_export(a1));

	if (!(a2 = arena_open(&arena_defaults,pool_export(p1))))
		err(EXIT_FAILURE,"arena_open(&arena_defaults,%p)",(void *)pool_export(p1));

	ap[0]	= pool_export(p0);
	ap[1]	= arena_export(a0);
	ap[2]	= pool_export(p1);

	ap[3]	= arena_export(a1);
	ap[4]	= pool_export(p2);
	ap[5]	= arena_export(a2);

	if (!(fp = fopen("/dev/urandom","r")))
		err(EXIT_FAILURE,"fopen(/dev/urandom,r)");

	srandom(time(0));

	for (i = 0; i < sizeof ap / sizeof *ap; i++) {
		for (j = 1; j <= 1<<10; j++) {
			/* random size from 1 -> j */
			n	= j % (unsigned long)random();
	
			if (!(d = ap[i]->malloc(ap[i],offsetof(struct data,bytes) + n,0)))
				errx(EXIT_FAILURE,"ap->malloc(%p,%d,%d): %s",(void *)ap[i],offsetof(struct data,bytes) + n,0,(ap[i]->strerror)(ap[i]));

			d->ap		= ap[i];
			d->order	= random();
			d->nbytes	= n;
	
			if (!(n == fread(d->bytes,1,d->nbytes,fp)))
				err(EXIT_FAILURE,"fread(%p,1,%d,/dev/urandom)",(void *)d->bytes,1,(int)d->nbytes,(void *)fp);
	
			d->checksum	= crc32(0,d->bytes,d->nbytes);
	
			LIST_FOREACH(dn,&list,le) {
				if (d->order > dn->order) {
					LIST_INSERT_AFTER(dn,d,le);
	
					d	= 0;
	
					break;
				}
			}
	
			if (d) {
				LIST_INSERT_HEAD(&list,d,le);
	
				d	= 0;
			}
		}

		for (d = LIST_FIRST(&list); d; d = dn) {
			dn	= LIST_NEXT(d,le);

			if (random() & 0x01)
				continue;

			LIST_REMOVE(d,le);

			if (d->checksum != crc32(0,d->bytes,d->nbytes))
				errx(EXIT_FAILURE,"bad crc32");

			d->ap->free(d->ap,d);
		}
	}

	while ((d = LIST_FIRST(&list))) {
		LIST_REMOVE(d,le);

		if (d->checksum != crc32(0,d->bytes,d->nbytes))
			errx(EXIT_FAILURE,"bad crc32");
	}

	return 0;
}
