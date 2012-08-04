#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "core/model.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct 
{
  char* name;
  char* usage;
  model_vtable vtable;
} model_metadata_t;

// Model database.
static model_metadata_t** model_db = NULL;
static int model_db_size = 0;
static int model_db_cap = 32;

struct model_t 
{
  void* context;
  char* name;
  char* usage;
  model_vtable vtable;
};

static void destroy_model_db()
{
  for (int i = 0; i < model_db_size; ++i)
  {
    free(model_db[i]->name);
    free(model_db[i]->usage);
    free(model_db[i]);
  }
  free(model_db);
}

void register_model(const char* name, 
                    const char* usage_string,
                    model_vtable vtable)
{
  ASSERT(name != NULL);
  ASSERT(usage_string != NULL);
  model_metadata_t* metadata = malloc(sizeof(model_metadata_t));
  metadata->name = strdup(name);
  metadata->usage = strdup(usage_string);
  metadata->vtable = vtable;

  // Stash the model metadata!
  if (model_db == NULL)
  {
    model_db = malloc(model_db_cap*sizeof(model_metadata_t));
    arbi_atexit(&destroy_model_db);
  }
  else
  {
    if (model_db_size == model_db_cap)
    {
      model_db_cap *= 2;
      model_db = realloc(model_db, model_db_cap*sizeof(model_metadata_t));
    }
  }
  model_db[model_db_size] = metadata;
  model_db_size++;
}
 
char** registered_models(int* num_models)
{
  *num_models = model_db_size;
  char** array = malloc(model_db_size*sizeof(char*));
  for (int i = 0; i < model_db_size; ++i)
    array[i] = model_db[i]->name;
  return array;
}

bool model_exists(const char* name)
{
  for (int i = 0; i < model_db_size; ++i)
  {
    if (!strcmp(model_db[i]->name, name))
      return true;
  }
  return false;
}

model_t* model_new(const char* name, options_t* options)
{
  // Dig up the model.
  model_metadata_t* metadata = NULL;
  for (int i = 0; i < model_db_size; ++i)
  {
    if (!strcmp(model_db[i]->name, name))
      metadata = model_db[i];
  }
  if (metadata == NULL) return NULL;

  model_t* model = malloc(sizeof(model_t));
  model->context = metadata->vtable.ctor(options);
  model->name = metadata->name;     // Borrowed from model_db
  model->usage = metadata->usage;   // Borrowed from model_db
  model->vtable = metadata->vtable;
  return model;
}

void model_free(model_t* model)
{
  if ((model->context != NULL) && (model->vtable.dtor != NULL))
    model->vtable.dtor(model->context);
  free(model);
}

char* model_name(model_t* model)
{
  return model->name;
}

// Print usage information for the model to the given file stream.
void model_usage(model_t* model, FILE* stream)
{
  fprintf(stream, "%s\n", model->usage);
}

void model_run_benchmark(model_t* model, const char* benchmark, options_t* opts)
{
  if (model->vtable.run_benchmark != NULL)
    model->vtable.run_benchmark(model->context, benchmark, opts);
  else
  {
    char err[1024];
    snprintf(err, 1024, "No benchmarks are defined for the model '%s'.", model->name);
    arbi_error(err);
  }
}

// Initialize the model at the given time.
void model_init(model_t* model, double t)
{
  model->vtable.init(model->context, t);
}

// Returns the largest permissible time step that can be taken by the model
// starting at time t.
double model_max_dt(model_t* model, double t, char* reason)
{
  if (model->vtable.max_dt != NULL)
    return model->vtable.max_dt(model->context, t, reason);
  else
  {
    strcpy(reason, "No time step constraints.");
    return FLT_MAX;
  }
}

void model_advance(model_t* model, double t, double dt)
{
  model->vtable.advance(model->context, t, dt);
}

void model_load(model_t* model, io_interface_t* io, double* t, int* step)
{
  model->vtable.load(model->context, io, t, step);
}

void model_dump(model_t* model, io_interface_t* io, double t, int step)
{
  model->vtable.dump(model->context, io, t, step);
}

void model_plot(model_t* model, plot_interface_t* plot, double t, int step)
{
  model->vtable.plot(model->context, plot, t, step);
}

#ifdef __cplusplus
}
#endif

