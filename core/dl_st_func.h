#ifndef POLYMEC_DL_ST_FUNC_H
#define POLYMEC_DL_ST_FUNC_H

#include "st_func.h"

// Construct a space-time function from a dynamically-loaded module
// that is either compiled in-place or loaded into Polymec's runtime database.
st_func_t* dl_st_func_new(const char* name);

// Sets the C compiler to use to build dynamically-loaded functions.
void dl_st_func_set_compiler(const char* cc,
                             const char* cflags);

// Sets the directory in which shared objects are stored.
void dl_st_func_set_so_dir(const char* path);

// Registers the dynamically-loaded function with the given name with 
// Polymec, building its shared object using the given source code.
void dl_st_func_register(const char* name, const char* source_code);

#endif

