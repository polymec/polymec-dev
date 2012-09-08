#ifndef ARBI_HASHMAP_H
#define ARBI_HASHMAP_H

#include "arbi.h"
#include "tommyhash.h"

#ifdef __cplusplus
extern "C" {
#endif

// hashmap_t is a hash table based off of the tommyds hashlin data structure,
// with a few niceties thrown in. For example, all elements are of the same 
// type in our hashmap, so we can register comparison and hash functions 
// upon creation of the map instead of on a per-element basis.
typedef struct hashmap_t hashmap_t;
typedef tommy_hash_t hashmap_hash_func;
typedef tommy_compare_func hashmap_comp_func;
typedef void (*hashmap_element_dtor)(void*);

// Create a new, empty hash map that uses the given comparison function for 
// keys and the given destructor for its elements.
hashmap_t* hashmap_new(hashmap_comp_func comp, hashmap_element_dtor dtor);

// Create a new, empty hash map with the given hashing function.
hashmap_t* hashmap_new_with_hash_func(hashmap_comp_func comp, hashmap_element_dtor dtor, hashmap_hash_func hash);

// Frees resources for the given hash map.
void hashmap_free(hashmap_t* hashmap);

// Empties the hashmap.
void hashmap_clear(hashmap_t* hashmap);

// Inserts the given datum into the hashmap.
void hashmap_insert(hashmap_t* hashmap, void* key, void* value);

// Removes the datum with the given key from the hashmap.
void hashmap_delete(hashmap_t* hashmap, void* key);

// Returns the value associated with the given key, or NULL if it is not found.
void* hashmap_find(hashmap_t* hashmap, void* key);

// Returns the value associated with the given key (or NULL), and deletes it 
// from the hash using a single lookup.
void* hashmap_find_and_delete(hashmap_t* hashmap, void* key);

// Returns the number of items in the hashmap.
int hashmap_size(hashmap_t* hashmap);

// Returns a list containing the keys in the hashmap. Do not free this list.
list_t* hashmap_keys(hashmap_t* hashmap);

#ifdef __cplusplus
}
#endif

#endif
