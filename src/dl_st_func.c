#include <dlfcn.h>
#include <string.h>
#include <stdlib.h>
#include "dl_st_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// C compiler and flags to use.
static char* _cc     = NULL;
static char* _cflags = NULL;

// A database of compiled functions.
static char** _func_names   = NULL;
static char** _func_sources = NULL;
static char** _func_objects = NULL;
static int    _num_funcs = 0;
static int    _capacity = 0;

// The shared object directory.
static char* _so_dir = NULL;

// Called at exit.
static void dl_st_atexit()
{
  if (_func_names != NULL)
  {
    for (int i = 0; i < _num_funcs; ++i)
    {
      free(_func_names[i]);
      free(_func_sources[i]);
      free(_func_objects[i]);
    }
    free(_func_names);
    free(_func_sources);
    free(_func_objects);
  }
  if (_so_dir != NULL)
    free(_so_dir);
  if (_cc != NULL)
  {
    free(_cc);
    free(_cflags);
  }
}

static char* func_object(const char* func_name)
{
  int i = 0;
  while (i < _num_funcs)
  {
    if (!strcmp(_func_names[i], func_name))
      return _func_objects[i];
    ++i;
  }
  return NULL;
}

static char* compile_so(const char* name, const char* source_file)
{
  if (_so_dir == NULL)
  {
    // By default, use the CWD.
    // FIXME: This is lousy.
    _so_dir = strdup(".");
  }
  char obj_name[128];
  snprintf(obj_name, 128, "%s/%s.so", _so_dir, name);
  char cmd[1024];
  snprintf(cmd, 1024, "%s %s -o %s %s", _cc, _cflags, obj_name, source_file);
  system(cmd);
  return strdup(obj_name);
}

void dl_st_func_set_so_dir(const char* path)
{
  if (_so_dir != NULL)
    free(_so_dir);
  _so_dir = strdup(path);
}

void dl_st_func_register(const char* name, const char* source_file)
{
  arbi_atexit(&dl_st_atexit);
  if (_func_names == NULL)
  {
    _capacity = 32;
    _func_names = malloc(_capacity*sizeof(char*));
    _func_sources = malloc(_capacity*sizeof(char*));
    _func_objects = malloc(_capacity*sizeof(char*));
  }
  else if (_num_funcs == _capacity)
  {
    _capacity *= 2;
    _func_names = realloc(_func_names, _capacity*sizeof(char*));
    _func_sources = realloc(_func_sources, _capacity*sizeof(char*));
    _func_objects = realloc(_func_objects, _capacity*sizeof(char*));
  }

  // Try to compile the thing to a shared object.
  _func_objects[_num_funcs] = compile_so(name, source_file);
  if (_func_objects[_num_funcs] == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "dl_st_func_register: Could not compile %s into a shared object.", source_file);
    arbi_error(err);
    return;
  }

  // Register the function and its source file.
  _func_names[_num_funcs] = strdup(name);
  _func_sources[_num_funcs] = strdup(source_file);

  ++_num_funcs;
}

static void dl_st_dtor(void* handle)
{
  dlclose(handle);
}

st_func_t* dl_st_func_new(const char* name)
{
  if (_cc == NULL)
  {
    arbi_error("dl_st_func_new: No C compiler is set!");
    return NULL;
  }

  // Search for the object that implements this function.
  char* obj_path = func_object(name);
  if (obj_path == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "dl_st_func_new: Could not find dynamically-loadable st_func '%s'", name);
    arbi_error(err);
    return NULL;
  }

  // Okay--try to load it.
  void* handle = dlopen(obj_path, RTLD_LOCAL | RTLD_LAZY);
  if (handle == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "dl_st_func_new: Could not load dynamically-loadable st_func '%s'", name);
    arbi_error(err);
    return NULL;
  }

  // Look for the proper symbols.
  st_vtable vtable;
  vtable.eval = dlsym(handle, "eval");
  vtable.dtor = &dl_st_dtor;
  int* h = dlsym(handle, "homogeneous");
  int* c = dlsym(handle, "constant");
  int* n = dlsym(handle, "num_comp");
  st_func_homogeneity_t homogenity = ST_INHOMOGENEOUS;
  if ((h != NULL) && (*h == 1))
    homogenity = ST_HOMOGENEOUS;
  st_func_constancy_t constancy = ST_NONCONSTANT;
  if ((c != NULL) && (*c == 1))
    constancy = ST_CONSTANT;
  int num_comp = (n != NULL) ? *n : 1;

  return st_func_new(name, handle, vtable, homogenity, constancy, num_comp);

}

void dl_st_func_set_compiler(const char* cc,
                             const char* cflags)
{
  if (_cc != NULL)
  {
    free(_cc);
    free(_cflags);
  }
  _cc = strdup(cc);
  _cflags = strdup(cflags);
}

#ifdef __cplusplus
}
#endif

