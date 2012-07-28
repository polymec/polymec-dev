#ifndef ARBI_HEAP_H
#define ARBI_HEAP_H

#include "arbi.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct heap_t heap_t;
typedef int (*heap_item_cmp)(void*, void*);
typedef void (*heap_item_dtor)(void*);

heap_t* heap_new(heap_item_cmp cmp, heap_item_dtor dtor);
void    heap_free(heap_t* heap);
void*   heap_extract(heap_t* heap);
void*   heap_inspect(heap_t* heap);
void    heap_insert(heap_t* heap, void* item);
void    heap_delete(heap_t* heap);
void    heap_clear(heap_t* heap);
void    heap_merge(heap_t* heap, heap_t* other_heap);
void    heap_delinsert(heap_t* heap, void* item);

#ifdef __cplusplus
}
#endif

#endif
