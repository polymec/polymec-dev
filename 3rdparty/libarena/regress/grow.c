#include <stddef.h> /* size_t */
#include <stdlib.h> /* strtoul(3) */
#include <unistd.h> /* getopt(3) */
#include <err.h>    /* err(3) */

#include "pool.h"
#include "proto.h"


int main(int argc, char **argv) {
	extern char *optarg;
	int optc;
	size_t align = 0;
	size_t length = 1UL << 6;
	struct pool *P;
	char *buf = NULL, *tmp;
	size_t i;

	while (-1 != (optc = getopt(argc, argv, "a:n:"))) {
		switch (optc) {
		case 'a':
			align = strtoul(optarg, NULL, 0);

			break;
		case 'n':
			length = strtoul(optarg, NULL, 0);

			break;
		}
	}

	P = pool_open(&pool_defaults, NULL);

	for (i = 1; i <= length; i++) {
		if (!(tmp = pool_realloc(P, buf, i, align)))
			err(1, "oops");

		buf = tmp;
		buf[i - 1] = 0x80;
	}

	pool_put(P, buf);

	pool_close(P);

	return 0;
} /* main() */

