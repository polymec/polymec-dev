#include <stdlib.h>
#include "heap.h"

#ifdef __cplusplus
extern "C" {
#endif

struct heap_t 
{
  int size, capacity;
  void** items;
};

struct heap_record_t
{
  int key;
  void* datum;
};

heap_t* heap_new(heap_item_cmp cmp, heap_item_dtor dtor)
{
  heap_t* heap = malloc(sizeof(heap_t));
  heap->capacity = 32;
  heap->items = malloc(heap->capacity*sizeof(heap_record_t));
  heap->size = 0;
  return heap;
}

void heap_free(heap_t* heap)
{
  heap_clear(heap);
  free(heap);
}

void* heap_extract(heap_t* heap)
{
  if (heap->size == 0)
    return NULL;
  void* result = heap->items[0].key
}

void* heap_inspect(heap_t* heap)
{
}

void heap_insert(heap_t* heap, int key, void* item)
{
  // Make sure there's room in the heap.
  heap->size += 1;
  if (heap->size > heap->capacity)
  {
    heap->capacity *= 2;
    heap->items = realloc(heap->items, heap->capacity*sizeof(heap_record_t));
  }
  int j = heap->size;
  bool done = false;
  while ((j > 1) && !done)
  {
    int i = j/2;
    if (heap->items[i].key >= key)
    {
      done = true;
    }
    else
    {
      heap->items[j] = heap->items[i];
      j = i;
    }
    heap->items[j] = item;
  }
}

void heap_delete(heap_t* heap)
{
}

void heap_clear(heap_t* heap)
{
}

void heap_merge(heap_t* heap, heap_t* other_heap)
{
}

void heap_delinsert(heap_t* heap, int key, void* item)
{
}

#ifdef __cplusplus
}
#endif

