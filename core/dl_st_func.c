#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <dlfcn.h>
#include "dl_st_func.h"

// Compiler defaults
#define xstr(s) str(s)
#define str(s) #s
#define DEFAULT_CC xstr(POLYMEC_C_COMPILER)
#define DEFAULT_CFLAGS "-shared"
#define INCLUDE_DIR xstr(POLYMEC_INCLUDE_DIR)
#define TP_INCLUDE_DIR xstr(POLYMEC_TP_INCLUDE_DIR)

#ifdef __cplusplus
extern "C" {
#endif

// C compiler and flags to use.
static char* _cc     = NULL;
static char* _cflags = NULL;

// A database of compiled functions.
static char** _func_names   = NULL;
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
      free(_func_objects[i]);
    }
    free(_func_names);
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

  if (_cc == NULL)
  {
    _cc = strdup(DEFAULT_CC);
    _cflags = strdup(DEFAULT_CFLAGS);
  }

  char obj_name[128];
  snprintf(obj_name, 128, "%s/%s.so", _so_dir, name);
  char cmd[1024];
  snprintf(cmd, 1024, "%s %s -o %s -I%s -I%s %s", _cc, _cflags, obj_name, INCLUDE_DIR, TP_INCLUDE_DIR, source_file);
  int status = system(cmd);
  if (status != 0)
  {
    char err[1024];
    snprintf(err, 1024, "compile_so: Compilation errors while building %s.", name);
    polymec_error(err);
    return NULL;
  }
  return strdup(obj_name);
}

static void preprocess_source(const char* source_code, char* pp_file)
{
  const char* header = 
    "#include <math.h>\n"
    "#include \"core/st_func.h\"\n\n"
    "#ifdef __cplusplus\n"
    "extern \"C\" {\n"
    "#endif\n\n";
  const char* footer = 
    "#ifdef __cplusplus\n"
    "}\n"
    "#endif\n";
  int fd = mkstemps(pp_file, 2);
  ASSERT(fd != -1);
  FILE* f = fdopen(fd, "w");
  fprintf(f, "%s%s%s", header, source_code, footer);
  fclose(f);
}

void dl_st_func_set_so_dir(const char* path)
{
  if (_so_dir != NULL)
    free(_so_dir);
  _so_dir = strdup(path);
}

void dl_st_func_register(const char* name, const char* source_code)
{
  polymec_atexit(&dl_st_atexit);
  if (_func_names == NULL)
  {
    _capacity = 32;
    _func_names = malloc(_capacity*sizeof(char*));
    _func_objects = malloc(_capacity*sizeof(char*));
  }
  else if (_num_funcs == _capacity)
  {
    _capacity *= 2;
    _func_names = realloc(_func_names, _capacity*sizeof(char*));
    _func_objects = realloc(_func_objects, _capacity*sizeof(char*));
  }

  // First, we preprocess the source to add some smarts.
  char pp_file[128];
  snprintf(pp_file, 128, "%sXXXXXX.c", name);
  preprocess_source(source_code, pp_file);

  // Now Try to compile the thing to a shared object.
  _func_objects[_num_funcs] = compile_so(name, pp_file);
  unlink(pp_file);
  if (_func_objects[_num_funcs] == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "dl_st_func_register: Could not compile %s into a shared object.", name);
    polymec_error(err);
    return;
  }

  // Register the function.
  _func_names[_num_funcs] = strdup(name);
  ++_num_funcs;
}

static void dl_st_dtor(void* handle)
{
  dlclose(handle);
}

st_func_t* dl_st_func_new(const char* name)
{
  // Search for the object that implements this function.
  char* obj_path = func_object(name);
  if (obj_path == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "dl_st_func_new: Could not find dynamically-loadable st_func '%s'", name);
    polymec_error(err);
    return NULL;
  }

  // Okay--try to load it.
  void* handle = dlopen(obj_path, RTLD_LOCAL | RTLD_LAZY);
  if (handle == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "dl_st_func_new: Could not load dynamically-loadable st_func '%s'", name);
    polymec_error(err);
    return NULL;
  }

  // Look for the proper symbols.
  st_vtable vtable;
  vtable.eval = dlsym(handle, "eval");
  if (vtable.eval == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "dl_st_func_new: Could not load eval function for st_func '%s'", name);
    polymec_error(err);
    return NULL;
  }
  vtable.dtor = &dl_st_dtor;
  bool* h = dlsym(handle, "homogeneous");
  bool* c = dlsym(handle, "constant");
  int* n = dlsym(handle, "num_comp");
  st_func_homogeneity_t homogenity = ST_INHOMOGENEOUS;
  if ((h != NULL) && (*h))
    homogenity = ST_HOMOGENEOUS;
  st_func_constancy_t constancy = ST_NONCONSTANT;
  if ((c != NULL) && (*c))
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

