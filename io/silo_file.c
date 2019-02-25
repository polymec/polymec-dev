// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <sys/stat.h>
#include <dirent.h>
#include "core/polymec.h"

#if POLYMEC_HAVE_MPI
// We use this MPI tag to pass batons back and forth between processes.
static int SILO_FILE_MPI_TAG = 1111;

// Poor Man's Parallel I/O stuff.
#include "pmpio.h"
#endif

#include "silo.h"
#include "core/arch.h"
#include "core/logging.h"
#include "core/array.h"
#include "core/array_utils.h"
#include "core/timer.h"
#include "io/silo_file.h"

#if POLYMEC_HAVE_DOUBLE_PRECISION
#define SILO_FLOAT_TYPE DB_DOUBLE
#else
#define SILO_FLOAT_TYPE DB_FLOAT
#endif

//-------------------------------------------------------------------------
// The following unpublished functions are not part of the formal silo_file 
// API, and offer access to low-level Silo internals that can be used 
// directly in other polymec-based libraries.
//-------------------------------------------------------------------------

// Access the underlying SILO file MPI communicator.
MPI_Comm silo_file_comm(silo_file_t* file);

// Access the underlying SILO file descriptor.
DBfile* silo_file_dbfile(silo_file_t* file);

// Creates a SILO options list from the given metadata object.
DBoptlist* optlist_from_metadata(field_metadata_t* metadata, int component);

// Clones a metadata-related SILO options list.
DBoptlist* optlist_clone(DBoptlist* optlist);

// Frees memory associated with a metadata-related SILO options list.
void optlist_free(DBoptlist* optlist);

// Writes metadata identifying a subdomain of a global polymesh of the given 
// type to the given SILO file. This is used for objects that appear on more 
// than 1 process via domain decomposition.
void silo_file_add_subdomain_mesh(silo_file_t* file,
                                  const char* mesh_name, 
                                  int silo_mesh_type,
                                  DBoptlist* optlist);

// Writes metadata identifying a subdomain of a global field of the given 
// type to the given SILO file. This is used for objects that appear on more 
// than 1 process via domain decomposition.
void silo_file_add_subdomain_field(silo_file_t* file,
                                   const char* mesh_name, 
                                   const char* field_name,
                                   int silo_var_type,
                                   DBoptlist* optlist);

// This pushes the given directory onto the silo file's directory stack.
void silo_file_push_dir(silo_file_t* file, const char* dir);

// This pushes the domain directory onto the silo file's directory stack.
void silo_file_push_domain_dir(silo_file_t* file);

// This pops the last directory off the silo file's directory stack.
void silo_file_pop_dir(silo_file_t* file);

// This provides access to a "scratch space" in memory that allows named 
// temporary data to be stored and passed between methods. 
string_ptr_unordered_map_t* silo_file_scratch(silo_file_t* file);

// Writes field metadata to a file.
void silo_file_write_field_metadata(silo_file_t* file,
                                    const char* md_name,
                                    field_metadata_t* md);
// Reads field metadata from a file. 
void silo_file_read_field_metadata(silo_file_t* file,
                                   const char* md_name,
                                   field_metadata_t* md);

// Checks for the existence of field metadata.
bool silo_file_contains_field_metadata(silo_file_t* file,
                                       const char* md_name);

//-------------------------------------------------------------------------
// End unpublished functions
//-------------------------------------------------------------------------

// SILO Compression -- global parameter. Disabled by default.
static int _silo_compression_level = -1; 
static bool _silo_compression_level_changed = false;

static void silo_atexit(void)
{
  DBSetCompression(NULL);
}

static void silo_set_compression(void)
{
  if (_silo_compression_level_changed)
  {
    char options[1024];
    snprintf(options, 1024, "ERRMODE=FALLBACK,METHOD=GZIP,LEVEL=%d", _silo_compression_level);
    DBSetCompression(options);
    _silo_compression_level_changed = false;
    polymec_atexit(silo_atexit);
  }
}

void silo_enable_compression(int level)
{
  ASSERT(level >= 0);
  ASSERT(level <= 9);
  if (level != _silo_compression_level)
  {
    _silo_compression_level = level; 
    _silo_compression_level_changed = true;
  }
}

DBoptlist* optlist_from_metadata(field_metadata_t* md, int component)
{
  DBoptlist* optlist = NULL;
  optlist = DBMakeOptlist(4);
  const char* name = field_metadata_name(md, component);
  if (name != NULL)
    DBAddOption(optlist, DBOPT_LABEL, string_dup(name));
  const char* units = field_metadata_units(md, component);
  if (units != NULL)
    DBAddOption(optlist, DBOPT_UNITS, string_dup(units));
  int* conserved = polymec_malloc(sizeof(int));
  *conserved = field_metadata_conserved(md, component);
  int* extensive = polymec_malloc(sizeof(int));
  *extensive = field_metadata_extensive(md, component);
  DBAddOption(optlist, DBOPT_CONSERVED, conserved);
  DBAddOption(optlist, DBOPT_EXTENSIVE, extensive);
  return optlist;
}

// These helpers are used in the read/write operations.
DBoptlist* optlist_clone(DBoptlist* optlist)
{
  DBoptlist* clone = NULL;
  if (optlist != NULL)
  {
    clone = DBMakeOptlist(4);
    char* label = DBGetOption(optlist, DBOPT_LABEL);
    if (label != NULL)
      DBAddOption(clone, DBOPT_LABEL, string_dup(label));
    char* units = DBGetOption(optlist, DBOPT_UNITS);
    if (units != NULL)
      DBAddOption(clone, DBOPT_UNITS, string_dup(units));
    int* conserved = DBGetOption(optlist, DBOPT_CONSERVED);
    if (conserved != NULL)
    {
      int* my_conserved = polymec_malloc(sizeof(int));
      *my_conserved = *conserved;
      DBAddOption(clone, DBOPT_CONSERVED, my_conserved);
    }
    int* extensive = DBGetOption(optlist, DBOPT_EXTENSIVE);
    if (extensive != NULL)
    {
      int* my_extensive = polymec_malloc(sizeof(int));
      *my_extensive = *extensive;
      DBAddOption(clone, DBOPT_EXTENSIVE, my_extensive);
    }
  }
  return clone;
}

void optlist_free(DBoptlist* optlist)
{
  if (optlist != NULL)
  {
    char* label = DBGetOption(optlist, DBOPT_LABEL);
    if (label != NULL)
      string_free(label);
    char* units = DBGetOption(optlist, DBOPT_UNITS);
    if (units != NULL)
      string_free(units);
    int* conserved = DBGetOption(optlist, DBOPT_CONSERVED);
    if (conserved != NULL)
      polymec_free(conserved);
    int* extensive = DBGetOption(optlist, DBOPT_EXTENSIVE);
    if (extensive != NULL)
      polymec_free(extensive);
    DBFreeOptlist(optlist);
  }
}

#if POLYMEC_HAVE_MPI
static void* pmpio_create_file(const char* filename,
                               const char* dir_name,
                               void* user_data)
{
  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, DB_CLOBBER, DB_LOCAL, NULL, driver);
  if (strcmp(dir_name, "/") != 0)
    DBMkDir(file, dir_name);
  return (void*)file;
}

static void* pmpio_open_file(const char* filename, 
                             const char* dir_name,
                             PMPIO_iomode_t iomode, 
                             void* user_data)
{
  int driver = DB_HDF5;
  DBfile* file;
  if (iomode == PMPIO_WRITE)
  { 
    FILE* f = fopen(filename, "r");
    if (f == NULL)
      file = DBCreate(filename, DB_CLOBBER, DB_LOCAL, NULL, driver);
    else
    {
      fclose(f);
      file = DBOpen(filename, driver, DB_APPEND);
    }
    if (strcmp(dir_name, "/") != 0)
      DBMkDir(file, dir_name);
  }
  else
    file = DBOpen(filename, driver, DB_READ);
  return (void*)file;
}

static void pmpio_close_file(void* file, void* user_data)
{
  DBClose((DBfile*)file);
}

// Object representing data in a subdomain.
typedef struct
{
  char* name;
  int type;
  DBoptlist* optlist;
} subdomain_mesh_t;

static subdomain_mesh_t* subdomain_mesh_new(const char* mesh_name, 
                                            int mesh_type, 
                                            DBoptlist* optlist)
{
  subdomain_mesh_t* mesh = polymec_malloc(sizeof(subdomain_mesh_t));
  mesh->name = string_dup(mesh_name);
  mesh->type = mesh_type;
  mesh->optlist = optlist_clone(optlist);
  return mesh;
}

static void subdomain_mesh_free(subdomain_mesh_t* mesh)
{
  polymec_free(mesh->name);
  optlist_free(mesh->optlist);
  polymec_free(mesh);
}

// Object representing data in a multi-mesh.
typedef struct
{
  char* mesh_name;
  char* name;
  int type;
  DBoptlist* optlist;
} subdomain_field_t;

// Constructors for various subdomain-related objects.
static subdomain_field_t* subdomain_field_new(const char* mesh_name,
                                              const char* field_name,
                                              int field_type,
                                              DBoptlist* optlist)
{
  subdomain_field_t* field = polymec_malloc(sizeof(subdomain_field_t));
  field->mesh_name = string_dup(mesh_name);
  field->name = string_dup(field_name);
  field->type = field_type;
  field->optlist = optlist_clone(optlist);
  return field;
}

static void subdomain_field_free(subdomain_field_t* field)
{
  polymec_free(field->mesh_name);
  polymec_free(field->name);
  optlist_free(field->optlist);
  polymec_free(field);
}

#endif

bool silo_file_query(const char* file_prefix,
                     const char* directory,
                     int* num_files,
                     int* num_mpi_processes,
                     int_slist_t* steps)
{
  START_FUNCTION_TIMER();

  // No blank strings allowed for queries.
  ASSERT(strlen(file_prefix) > 0);
  ASSERT(strlen(directory) > 0);

  // Rank 0 does all the dirty work.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
  {
    // Inspect the given directory's contents.
    if (!directory_exists(directory))
    {
      log_urgent("Can't query non-existent directory: %s", directory);
      STOP_FUNCTION_TIMER();
      return false;
    }
    string_slist_t* files_in_dir = files_within_directory(directory);

    // Try to find a master file or single data file.
    char data_file[FILENAME_MAX+1];
    data_file[0] = '\0';
    bool found_steps = false;
    {
      string_slist_node_t* node = files_in_dir->front;
      while (node != NULL)
      {
        char path[FILENAME_MAX+1];
        snprintf(path, FILENAME_MAX, "%s", file_prefix); 
        if (strstr(node->value, path) && strstr(node->value, "silo"))
        {
          strncpy(data_file, node->value, FILENAME_MAX);
          snprintf(data_file, FILENAME_MAX, "%s/%s", directory, node->value);
          snprintf(path, FILENAME_MAX, "%s-", file_prefix); 
          if (strstr(node->value, path))
            found_steps = true;
          break;
        }
        node = node->next;
      }
    }
    if (strlen(data_file) == 0)
    {
      log_urgent("silo_file_query: Could not find %s/%s-*.silo.", directory, file_prefix);
      string_slist_free(files_in_dir);
      STOP_FUNCTION_TIMER();
      return false;
    }

    // Open up the file and see whether it's a master file or a serial data file.
    int driver = DB_HDF5;
    DBfile* file = DBOpen(data_file, driver, DB_READ);

    // What's in there?
    bool is_master = DBInqVarExists(file, "POLYMEC_SILO_MASTER_FILE");
    if (is_master)
    {
      // How many MPI processes were used to construct the data set?
      DBtoc* toc = DBGetToc(file); 
      int my_num_mpi_procs = -1;
      for (int f = 0; f < toc->nmultimesh; ++f)
      {
        DBmultimesh* subdomain_mesh = DBGetMultimesh(file, toc->multimesh_names[f]);
        if (my_num_mpi_procs == -1)
          my_num_mpi_procs = subdomain_mesh->nblocks;
        ASSERT(my_num_mpi_procs == subdomain_mesh->nblocks);
        DBFreeMultimesh(subdomain_mesh);
      }
      *num_mpi_processes = my_num_mpi_procs;

      // How many files are in the data set?
      if (DBInqVarExists(file, "num_files"))
        DBReadVar(file, "num_files", num_files);
      else
      {
        log_urgent("silo_file_query: Could not read number of files in set.");
        DBClose(file);
        STOP_FUNCTION_TIMER();
        return false;
      }
    }
    else
    {
      if (!DBInqVarExists(file, "POLYMEC_SILO_FILE"))
      {
        log_urgent("silo_file_query: Invalid Silo file.");
        DBClose(file);
        STOP_FUNCTION_TIMER();
        return false;
      }
      ASSERT(DBInqVarExists(file, "num_mpi_procs"));

      // This dataset consists of a single file.
      *num_files = 1;
      DBReadVar(file, "num_mpi_procs", num_mpi_processes);
    }

    DBClose(file);

    // Search for available steps.
    if ((steps != NULL) && found_steps)
    {
      int_slist_clear(steps);
      string_slist_node_t* node = files_in_dir->front;
      while (node != NULL)
      {
        char path[FILENAME_MAX+1];
        snprintf(path, FILENAME_MAX, "%s-", file_prefix); 
        char* p1 = strstr(node->value, path);
        char* p2 = strstr(node->value, ".silo");
        if ((p1 != NULL) && (p2 != NULL))
        {
          char* c = p1 + strlen(path);
          char num[p2-c+1];
          strncpy(num, c, p2-c);
          num[p2-c] = '\0';
          if (string_is_number(num))
            int_slist_append(steps, atoi(num));
        }
        node = node->next;
      }
    }

    // Clean up.
    string_slist_free(files_in_dir);
  }

  // Now spread the word to other processes.
  int num_steps = (steps != NULL) ? (int)steps->size : 0;
  int data[3] = {*num_files, *num_mpi_processes, num_steps};
  MPI_Bcast(data, 3, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank != 0)
  {
    *num_files = data[0];
    *num_mpi_processes = data[1];
    num_steps = data[2];
  }

  if (steps != NULL)
  {
    int step_buffer[num_steps];
    if (rank == 0)
    {
      // Spit the steps into the array and sort them.
      int_slist_node_t* node = steps->front;
      int i = 0;
      while (node != NULL)
      {
        step_buffer[i++] = node->value;
        node = node->next;
      }
      ASSERT(i == num_steps);
      int_qsort(step_buffer, num_steps);
    }
    MPI_Bcast(step_buffer, num_steps, MPI_INT, 0, MPI_COMM_WORLD);

    // Shuffle them back into our linked list.
    int_slist_clear(steps);
    for (int i = 0; i < num_steps; ++i)
      int_slist_append(steps, step_buffer[i]);
  }
  STOP_FUNCTION_TIMER();
  return true;
}

struct silo_file_t 
{
  // File data.
  DBfile* dbfile;

  // Metadata.
  char prefix[FILENAME_MAX+1], directory[FILENAME_MAX+1], filename[FILENAME_MAX+1];
  int step;
  real_t time;
  int mode; // Open for reading (DB_READ) or writing (DB_CLOBBER)? 
  string_ptr_unordered_map_t* expressions;

  // Directory stack.
  string_slist_t* dirs;

  // Scratch space for storing named temporary data.
  string_ptr_unordered_map_t* scratch;

  MPI_Comm comm;
#if POLYMEC_HAVE_MPI
  // Stuff for poor man's parallel I/O.
  PMPIO_baton_t* baton;
  int num_files, mpi_tag, nproc, rank, group_rank, rank_in_group;

  // Data appearing on more than one proc within a domain decomposition.
  ptr_array_t* subdomain_meshes;
  ptr_array_t* subdomain_fields;
#endif
};

// Expression struct.
typedef struct
{
  int type;
  char* definition;
} silo_expression_t;

static void write_expressions_to_file(silo_file_t* file, DBfile* dbfile)
{
  ASSERT(file->mode == DB_CLOBBER);
  ASSERT(file->expressions != NULL);

  if (file->expressions->size > 0)
  {
    // Write out the expressions.
    int pos = 0, k = 0;
    char* name;
    void* val;
    char* names[file->expressions->size];
    int types[file->expressions->size];
    char* defs[file->expressions->size];
    while (string_ptr_unordered_map_next(file->expressions, &pos, &name, &val))
    {
      names[k] = name;
      silo_expression_t* exp = val;
      types[k] = exp->type;
      defs[k] = exp->definition;
      ++k;
    }
    DBPutDefvars(dbfile, "expressions", file->expressions->size,
                 (const char* const*)names, (const int*)types, 
                 (const char* const*)defs, NULL);
  }
}

static void write_provenance_to_file(silo_file_t* file)
{
  // We harvest provenance information from polymec.
  char* provenance_str;
  size_t provenance_strlen;
  FILE* stream = open_memstream(&provenance_str, &provenance_strlen);
  polymec_provenance_fprintf(stream);
  fclose(stream);

  silo_file_push_dir(file, "/");
  int provenance_len = (int)(strlen(provenance_str));
  DBWrite(file->dbfile, "provenance_string", provenance_str, &provenance_len, 1, DB_CHAR);
  silo_file_pop_dir(file);

  // Note we have to use free here instead of string_free or polymec_free, 
  // since open_memstream uses vanilla malloc and realloc.
  free(provenance_str);
}

#if POLYMEC_HAVE_MPI
static void write_subdomains_to_file(silo_file_t* file)
{
  ASSERT(file->mode == DB_CLOBBER);
  START_FUNCTION_TIMER();

  if ((file->nproc == 1) || (file->rank_in_group != 0))
  {
    STOP_FUNCTION_TIMER();
    return;
  }
  int num_chunks = file->nproc / file->num_files;

  // Stick in step/time information if needed.
  DBoptlist* optlist = DBMakeOptlist(2);
  if (file->step >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &file->step);
  if (reals_equal(file->time, -REAL_MAX))
  {
    double t = (double)file->time;
    DBAddOption(optlist, DBOPT_DTIME, &t);
  }

  // Write out subdomain meshes.
  for (int i = 0; i < file->subdomain_meshes->size; ++i)
  {
    subdomain_mesh_t* mesh = file->subdomain_meshes->data[i];

    char* mesh_names[num_chunks];
    int mesh_types[num_chunks];
    for (int j = 0; j < num_chunks; ++j)
    {
      char mesh_name[FILENAME_MAX+1];
      snprintf(mesh_name, FILENAME_MAX, "domain_%d/%s", j, mesh->name);
      mesh_names[j] = string_dup(mesh_name);
      mesh_types[j] = mesh->type;
    }

    // Write the point mesh and variable data.
    DBPutMultimesh(file->dbfile, mesh->name, num_chunks, 
                   (char const* const*)mesh_names, mesh_types, optlist);

    // Clean up.
    for (int j = 0; j < num_chunks; ++j)
      polymec_free(mesh_names[j]);
  }
  DBFreeOptlist(optlist);

  // Subdomain fields.
  for (int i = 0; i < file->subdomain_fields->size; ++i)
  {
    subdomain_field_t* field = file->subdomain_fields->data[i];

    // Fields and associated meshes.
    char* field_names[num_chunks];
    int field_types[num_chunks];
    for (int j = 0; j < num_chunks; ++j)
    {
      // Field name.
      char field_name[FILENAME_MAX+1];
      snprintf(field_name, FILENAME_MAX, "domain_%d/%s", j, field->name);
      field_names[j] = string_dup(field_name);
      field_types[j] = field->type;
    }

    // Write the field data.
    DBPutMultivar(file->dbfile, field->name, num_chunks, 
                  (char const* const*)field_names, field_types, field->optlist); 

    // Clean up.
    for (int j = 0; j < num_chunks; ++j)
      polymec_free(field_names[j]);
  }
  STOP_FUNCTION_TIMER();
}

static void write_master_file(silo_file_t* file)
{
  ASSERT(file->mode == DB_CLOBBER);
  if (file->num_files == 1) return;

  START_FUNCTION_TIMER();

  // FIXME: Should change this to use Silo's name schemes for multi-block 
  // FIXME: objects when we start to Get Real Parallel.

  char master_file_name[FILENAME_MAX+1];
  if (file->step == -1)
    snprintf(master_file_name, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
  else
    snprintf(master_file_name, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, file->step);
  PMPIO_baton_t* baton = PMPIO_Init(1, PMPIO_WRITE, file->comm, 
                                    file->mpi_tag+1, pmpio_create_file, 
                                    pmpio_open_file, pmpio_close_file, NULL);
  log_debug("write_master_file: Waiting for baton...");
  DBfile* master = (DBfile*)PMPIO_WaitForBaton(baton, master_file_name, "/");
  log_debug("write_master_file: Got the baton...");

  // Write our stamp of approval.
  int one = 1;
  DBWrite(master, "POLYMEC_SILO_MASTER_FILE", &one, &one, 1, DB_INT);

  // Stick in step/time information if needed.
  DBoptlist* optlist = DBMakeOptlist(2);
  if (file->step >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &file->step);
  if (reals_equal(file->time, -REAL_MAX))
  {
    double t = (double)file->time;
    DBAddOption(optlist, DBOPT_DTIME, &t);
  }

  int num_files = file->num_files;
  int num_chunks = file->nproc / num_files;

  // Meshes.
  for (int i = 0; i < file->subdomain_meshes->size; ++i)
  {
    subdomain_mesh_t* mesh = file->subdomain_meshes->data[i];

    // Mesh.
    char* mesh_names[file->num_files*num_chunks];
    int mesh_types[file->num_files*num_chunks];
    for (int j = 0; j < file->num_files; ++j)
    {
      for (int c = 0; c < num_chunks; ++c)
      {
        char mesh_name[FILENAME_MAX+1];
        mesh_types[num_chunks*j+c] = mesh->type;
        if (file->step == -1)
          snprintf(mesh_name, FILENAME_MAX, "%d/%s.silo:/domain_%d/%s", j, file->prefix, c, mesh->name);
        else
          snprintf(mesh_name, FILENAME_MAX, "%d/%s-%d.silo:/domain_%d/%s", j, file->prefix, file->step, c, mesh->name);
        mesh_names[num_chunks*j+c] = string_dup(mesh_name);
      }
    }

    // Write the subdomain_mesh.
    int stat = DBPutMultimesh(master, mesh->name, file->num_files*num_chunks, 
                              (char const* const*)mesh_names, mesh_types, 
                              optlist);
    if (stat == -1)
      polymec_error("Error writing multi-mesh to Silo master file %s.", master_file_name);

    // Clean up.
    for (int j = 0; j < num_files*num_chunks; ++j)
      polymec_free(mesh_names[j]);
  }

  // Variables.
  for (int i = 0; i < file->subdomain_fields->size; ++i)
  {
    subdomain_field_t* field = file->subdomain_fields->data[i];

    // Fields.
    char* field_names[file->num_files*num_chunks];
    int field_types[num_files*num_chunks];
    for (int j = 0; j < file->num_files; ++j)
    {
      for (int c = 0; c < num_chunks; ++c)
      {
        char field_name[FILENAME_MAX+1];
        if (file->step == -1)
          snprintf(field_name, FILENAME_MAX, "%d/%s.silo:/domain_%d/%s", j, file->prefix, c, field->name);
        else
          snprintf(field_name, FILENAME_MAX, "%d/%s-%d.silo:/domain_%d/%s", j, file->prefix, file->step, c, field->name);
        field_names[num_chunks*j+c] = string_dup(field_name);
        field_types[num_chunks*j+c] = field->type;
      }
    }

    // Write the subdomain_fieldiable data.
    DBPutMultivar(master, field->name, num_files*num_chunks, 
                  (char const* const *)field_names, field_types, optlist);

    // Clean up.
    for (int j = 0; j < num_files*num_chunks; ++j)
      polymec_free(field_names[j]);
  }

  // Finally, write the number of files to the master file.
  DBWrite(master, "num_files", &num_files, &one, 1, DB_INT);

  // Don't forget to write expressions to the master file.
  write_expressions_to_file(file, master);

  // NOTE: Master files don't get provenance information.

  DBFreeOptlist(optlist);

  log_debug("write_master_file: Handing off baton.");
  PMPIO_HandOffBaton(baton, (void*)master);
  PMPIO_Finish(baton);

  STOP_FUNCTION_TIMER();
}
#endif

// Sets the prefix for the file, stripping .silo off if it's there.
static void set_prefix(silo_file_t* file,
                       const char* prefix)
{
  char pre[FILENAME_MAX+1];
  strncpy(pre, prefix, FILENAME_MAX);
  char* suffix = strstr(pre, ".silo");
  if (suffix != NULL)
    suffix[0] = '\0';
  strcpy(file->prefix, pre);
}

silo_file_t* silo_file_new(MPI_Comm comm,
                           const char* file_prefix,
                           const char* directory,
                           int num_files,
                           int step,
                           real_t time)
{
  START_FUNCTION_TIMER();
  ASSERT((num_files == -1) || (num_files > 0));

  // Set compression if needed.
  silo_set_compression();

  silo_file_t* file = polymec_malloc(sizeof(silo_file_t));
  file->expressions = string_ptr_unordered_map_new();

  set_prefix(file, file_prefix);

  file->dirs = NULL;
  file->comm = comm;
#if POLYMEC_HAVE_MPI
  file->mpi_tag = SILO_FILE_MPI_TAG;
  MPI_Comm_size(file->comm, &file->nproc);
  MPI_Comm_rank(file->comm, &file->rank);
  if (num_files == -1)
    file->num_files = file->nproc;
  else
    file->num_files = num_files;
  ASSERT(file->num_files <= file->nproc);

  if (file->nproc > 1)
  {
    // We put the entire data set into a directory named after the 
    // prefix, and every process gets its own subdirectory therein.

    // Create the master directory if we need to.
    if ((strcmp(directory, ".") == 0) && (file->nproc > 1))
      polymec_error("silo_file_new: Multi-process filesets cannot be generated in the current working directory.");
    else if (strlen(directory) == 0)
      snprintf(file->directory, FILENAME_MAX, "%s_%dprocs", file->prefix, file->nproc);
    else
      strncpy(file->directory, directory, FILENAME_MAX);
    if (file->rank == 0)
    {
      if (strcmp(file->directory, ".") != 0)
        create_directory(file->directory, S_IRWXU | S_IRWXG);
      MPI_Barrier(file->comm);
    }
    else
      MPI_Barrier(file->comm);

    // Initialize poor man's I/O and figure out group ranks.
    file->baton = PMPIO_Init(file->num_files, PMPIO_WRITE, file->comm, file->mpi_tag, 
                             pmpio_create_file, pmpio_open_file, 
                             pmpio_close_file, 0);
    file->group_rank = PMPIO_GroupRank(file->baton, file->rank);
    file->rank_in_group = PMPIO_RankInGroup(file->baton, file->rank);

    // Determine a file name and directory name.
    if (file->num_files > 1)
    {
      // Create a subdirectory for each group.
      char group_dir_name[FILENAME_MAX+1];
      snprintf(group_dir_name, FILENAME_MAX, "%s/%d", file->directory, file->group_rank);
      if (file->rank_in_group == 0)
        create_directory(group_dir_name, S_IRWXU | S_IRWXG);

      if (step == -1)
        snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", group_dir_name, file->prefix);
      else
        snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", group_dir_name, file->prefix, step);
    }
    else
    {
      ASSERT(file->group_rank == 0);
      ASSERT(file->rank_in_group == file->rank);
      if (step == -1)
        snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
      else
        snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, step);
    }
    char silo_dir_name[FILENAME_MAX+1];
    snprintf(silo_dir_name, FILENAME_MAX, "domain_%d", file->rank_in_group);
    log_debug("silo_file_new: Opening %s for writing and waiting for baton...", file->filename);
    file->dbfile = (DBfile*)PMPIO_WaitForBaton(file->baton, file->filename, silo_dir_name);
    log_debug("silo_file_new: Got the baton...");

    file->subdomain_meshes = ptr_array_new();
    file->subdomain_fields = ptr_array_new();
  }
  else
  {
    file->group_rank = 0;
    file->rank_in_group = 0;

    if (strlen(directory) == 0)
      strncpy(file->directory, ".", FILENAME_MAX);
    else
      strncpy(file->directory, directory, FILENAME_MAX);

    if (step == -1)
      snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
    else
      snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, step);

    int driver = DB_HDF5;
    if (strcmp(file->directory, ".") != 0)
      create_directory(file->directory, S_IRWXU | S_IRWXG);
    log_debug("silo_file_new: Opening %s for writing...", file->filename);
    file->dbfile = DBCreate(file->filename, DB_CLOBBER, DB_LOCAL, NULL, driver);
  }
#else
  if (strlen(directory) == 0)
    strncpy(file->directory, ".", FILENAME_MAX);
  else
    strncpy(file->directory, directory, FILENAME_MAX);

  if (step == -1)
    snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
  else
    snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, step);

  int driver = DB_HDF5;
  create_directory(file->directory, S_IRWXU | S_IRWXG);
  log_debug("silo_file_new: Opening %s for writing...", file->filename);
  file->dbfile = DBCreate(file->filename, DB_CLOBBER, DB_LOCAL, NULL, driver);
#endif

  silo_file_push_dir(file, "/");
  file->mode = DB_CLOBBER;
  file->step = step;
  file->time = time;

  // Write the step and cycle.
  {
    int one = 1;
    DBWrite(file->dbfile, "dtime", &(file->time), &one, 1, DB_DOUBLE);
    DBWrite(file->dbfile, "cycle", &(file->step), &one, 1, DB_INT);
  }

  // Write our stamp of approval.
  int one = 1;
#if POLYMEC_HAVE_MPI
  if (file->rank_in_group == 0)
#endif
    DBWrite(file->dbfile, "POLYMEC_SILO_FILE", &one, &one, 1, DB_INT);

  // If we're writing to a single file, encode our number of processes.
#if POLYMEC_HAVE_MPI
  if ((file->num_files == 1) && (file->rank == 0))
#endif
  {
#if POLYMEC_HAVE_MPI
    DBWrite(file->dbfile, "num_mpi_procs", &file->nproc, &one, 1, DB_INT);
#else
    DBWrite(file->dbfile, "num_mpi_procs", &one, &one, 1, DB_INT);
#endif
  }

  // Initialize scratch space.
  file->scratch = string_ptr_unordered_map_new();

  STOP_FUNCTION_TIMER();
  return file;
}

static void show_provenance_on_debug_log(silo_file_t* file)
{
  // If we're printing debug messages, show provenance of the file.
  if (log_level() == LOG_DEBUG)
  {
    int provenance_len = DBGetVarLength(file->dbfile, "provenance_string");
    char* provenance = polymec_malloc(sizeof(char) * (provenance_len + 1));
    if (provenance_len > 0)
    {
      DBReadVar(file->dbfile, "provenance_string", provenance);
      provenance[provenance_len] = '\0';
    }

    // We may need to temporarily bump up the debug log's buffness.
    int max_message_size, flush_every;
    get_log_buffering(LOG_DEBUG, &max_message_size, &flush_every);
    set_log_buffering(LOG_DEBUG, (int)strlen(provenance)+10, 1); // SUPER DUPER!!!
    log_debug_literal(provenance); // No formatting.
    set_log_buffering(LOG_DEBUG, max_message_size, flush_every); // Back to normal

    string_free(provenance);
  }
}

silo_file_t* silo_file_open(MPI_Comm comm,
                            const char* file_prefix,
                            const char* directory,
                            int step, 
                            real_t* time)
{
  START_FUNCTION_TIMER();

  // Set compression if needed.
  silo_set_compression();

  silo_file_t* file = polymec_calloc(1, sizeof(silo_file_t));
  file->mode = DB_READ;
  file->step = -1;
  file->time = -REAL_MAX;
  file->expressions = NULL;

  set_prefix(file, file_prefix);

  int nproc = 1;
  file->comm = comm;
  file->dirs = NULL;
#if POLYMEC_HAVE_MPI
  file->mpi_tag = SILO_FILE_MPI_TAG;
  MPI_Comm_size(file->comm, &nproc); 
#endif

  // Provide a default for the directory for querying.
  if (strlen(directory) == 0)
  {
    if (nproc > 1)
      snprintf(file->directory, FILENAME_MAX, "%s_%dprocs", file->prefix, nproc);
    else
      strncpy(file->directory, ".", FILENAME_MAX);
  }
  else
    strncpy(file->directory, directory, FILENAME_MAX);

  // Query the dataset for the number of files and MPI processes and step.
  int num_files, num_mpi_procs;
  int_slist_t* steps = int_slist_new();
  if (!silo_file_query(file_prefix, file->directory, &num_files, &num_mpi_procs, steps))
  {
    int_slist_free(steps);
    log_info("silo_file_open: Invalid file.");
    STOP_FUNCTION_TIMER();
    return NULL;
  }

  log_debug("silo_file_open: Found file written by %d MPI processes.", num_mpi_procs);

  // For now, we only support reading files that were written with the same 
  // number of processes.
  if (nproc != num_mpi_procs)
  {
    log_urgent("silo_file_open: Cannot read file written by %d MPI processes "
               "into communicator with %d processes.", nproc, num_mpi_procs);
    STOP_FUNCTION_TIMER();
    return NULL;
  }

  // Check to see whether the requested step is available, or whether the 
  // latest one is requested (with -1).
  if (step >= 0)
  {
    bool step_found = false;
    int_slist_node_t* node = steps->front;
    while (node != NULL)
    {
      if (node->value == step)
      {
        step_found = true;
        break;
      }
      else if (node->value > step) // steps are sorted
        break;
      node = node->next;
    }
    if (!step_found)
    {
      log_urgent("silo_file_open: Step %d was not found for prefix '%s' in directory %s.", step, file->prefix, directory);
      int_slist_free(steps);
      polymec_free(file);
      STOP_FUNCTION_TIMER();
      return NULL;
    }
  }
  int_slist_free(steps);

#if POLYMEC_HAVE_MPI
  // The way these things are defined for a file has to do with how the 
  // file was generated, not how we are currently running.
  MPI_Comm_rank(file->comm, &file->rank); 
  file->num_files = num_files; // number of files in the data set.
  file->nproc = num_mpi_procs; // number of MPI procs used to write the thing.

  if (file->nproc > 1)
  {
    // Look for the master directory.
    if (strlen(directory) == 0)
      snprintf(file->directory, FILENAME_MAX, "%s_%dprocs", file->prefix, file->nproc);
    else
      strncpy(file->directory, directory, FILENAME_MAX);
    int dir_exists;
    if (file->rank == 0)
      dir_exists = (int)(directory_exists(file->directory));
    MPI_Bcast(&dir_exists, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!dir_exists)
    {
      log_urgent("silo_file_open: Master directory %s does not exist for file prefix %s.",
                 file->directory, file->prefix);
      polymec_free(file);
      STOP_FUNCTION_TIMER();
      return NULL;
    }

    // Initialize poor man's I/O and figure out group ranks.
    MPI_Barrier(file->comm);
    file->baton = PMPIO_Init(file->num_files, PMPIO_READ, file->comm, file->mpi_tag, 
                             pmpio_create_file, pmpio_open_file, pmpio_close_file, 0);
    file->group_rank = PMPIO_GroupRank(file->baton, file->rank);
    file->rank_in_group = PMPIO_RankInGroup(file->baton, file->rank);

    if (file->num_files > 1)
    {
      // Make sure a subdirectory exists for each group.
      char group_dir_name[FILENAME_MAX+1];
      snprintf(group_dir_name, FILENAME_MAX, "%s/%d", file->directory, file->group_rank);
      if (file->rank_in_group == 0)
      {
        DIR* group_dir = opendir(group_dir_name);
        if (group_dir == NULL)
        {
          log_urgent("silo_file_open: Group directory %s does not exist for file prefix %s.",
                     group_dir_name, file->prefix);
          polymec_free(file);
          STOP_FUNCTION_TIMER();
          return NULL;
        }
        else
          closedir(group_dir);
      }

      // Determine a file name and directory name.
      if (step == -1)
        snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", group_dir_name, file->prefix);
      else
        snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", group_dir_name, file->prefix, step);
    }
    else
    {
      ASSERT(file->group_rank == 0);
      ASSERT(file->rank_in_group == file->rank);
      if (step == -1)
        snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
      else
        snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, step);
    }

    char silo_dir_name[FILENAME_MAX+1];
    snprintf(silo_dir_name, FILENAME_MAX, "domain_%d", file->rank_in_group);
    log_debug("silo_file_open: Opening %s for reading and waiting for baton...", file->filename);
    file->dbfile = (DBfile*)PMPIO_WaitForBaton(file->baton, file->filename, silo_dir_name);
    log_debug("silo_file_open: Got the baton...");
  }
  else
  {
    if (strlen(directory) == 0)
      strncpy(file->directory, ".", FILENAME_MAX);
    else
      strncpy(file->directory, directory, FILENAME_MAX);

    if (step == -1)
      snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
    else
      snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, step);

    int driver = DB_HDF5;
    log_debug("silo_file_open: Opening %s for reading...", file->filename);
    file->dbfile = DBOpen(file->filename, driver, file->mode);
  }
#else
  if (strlen(directory) == 0)
    strncpy(file->directory, ".", FILENAME_MAX);
  else
    strncpy(file->directory, directory, FILENAME_MAX);

  if (step == -1)
    snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
  else
    snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, step);

  int driver = DB_HDF5;
  log_debug("silo_file_open: Opening %s for reading...", file->filename);
  file->dbfile = DBOpen(file->filename, driver, file->mode);
#endif

  silo_file_push_dir(file, "/");
  show_provenance_on_debug_log(file);

  // Get step/time information.
  if (DBInqVarExists(file->dbfile, "dtime"))
  {
    double dtime;
    DBReadVar(file->dbfile, "dtime", &dtime);
    file->time = (real_t)dtime;
  }
  else
    file->time = 0.0;
  if (DBInqVarExists(file->dbfile, "cycle"))
    DBReadVar(file->dbfile, "cycle", &file->step);
  else
    file->step = -1;

  if (time != NULL)
    *time = file->time;

  // Initialize scratch space.
  file->scratch = string_ptr_unordered_map_new();

  STOP_FUNCTION_TIMER();
  return file;
}

void silo_file_close(silo_file_t* file)
{
  START_FUNCTION_TIMER();
#if POLYMEC_HAVE_MPI
  if (file->nproc > 1)
  {
    // Finish working on this process.
    if (file->mode == DB_CLOBBER)
    {
      // Write multi-block objects to the file if needed.
      write_subdomains_to_file(file);
      write_expressions_to_file(file, file->dbfile);
      write_provenance_to_file(file);
    }

    log_debug("silo_file_close: Handing off baton.");
    PMPIO_HandOffBaton(file->baton, (void*)file->dbfile);
    PMPIO_Finish(file->baton);

    if (file->mode == DB_CLOBBER)
    {
      // Write the uber-master file containing any multiobjects if need be.
      write_master_file(file);
    }

    MPI_Barrier(file->comm);

    if (file->subdomain_meshes != NULL)
      ptr_array_free(file->subdomain_meshes);
    if (file->subdomain_fields != NULL)
      ptr_array_free(file->subdomain_fields);
  }
  else
  {
    if (file->mode == DB_CLOBBER)
    {
      write_expressions_to_file(file, file->dbfile);
      write_provenance_to_file(file);
    }
    DBClose(file->dbfile);
  }
#else
  // Write the file.
  if (file->mode == DB_CLOBBER)
  {
    write_expressions_to_file(file, file->dbfile);
    write_provenance_to_file(file);
  }
  DBClose(file->dbfile);
#endif

  log_debug("silo_file_close: Closed file.");

  // Clean up.
  if (file->dirs != NULL)
    string_slist_free(file->dirs);
  if (file->scratch != NULL)
    string_ptr_unordered_map_free(file->scratch);
  if (file->expressions != NULL)
    string_ptr_unordered_map_free(file->expressions);
  polymec_free(file);
  STOP_FUNCTION_TIMER();
}

static void silo_file_write_tags(silo_file_t* file, tagger_t* tagger, const char* tag_list_name)
{
  ASSERT(file->mode == DB_CLOBBER);

  // Pack the tags into a compound array.
  int_array_t* elem_lengths = int_array_new();
  string_array_t* elem_names = string_array_new();
  int_array_t* tag_data = int_array_new();

  int pos = 0, *tag;
  size_t tag_size;
  char* tag_name;
  while (tagger_next_tag(tagger, &pos, &tag_name, &tag, &tag_size))
  {
    int_array_append(elem_lengths, (int)tag_size);
    string_array_append(elem_names, tag_name);
    for (int i = 0; i < tag_size; ++i)
      int_array_append(tag_data, tag[i]);
  }

  // Write the compound array.
  if (elem_names->size > 0)
  {
    DBPutCompoundarray(file->dbfile, tag_list_name, 
                       (char const* const*)elem_names->data, elem_lengths->data,
                       (int)elem_names->size, tag_data->data, (int)tag_data->size, DB_INT, 0);
  }

  // Clean up.
  int_array_free(elem_lengths);
  string_array_free(elem_names);
  int_array_free(tag_data);
}

static void silo_file_read_tags(silo_file_t* file, const char* tag_list_name, tagger_t* tagger)
{
  if (DBInqVarExists(file->dbfile, tag_list_name))
  {
    DBcompoundarray* var = DBGetCompoundarray(file->dbfile, (char*)tag_list_name);
    ASSERT(var != NULL);
    int num_tags = var->nelems;
    char** tag_names = var->elemnames;
    int* tag_sizes = var->elemlengths;
    int* array = var->values;
    int j = 0;
    for (int i = 0; i < num_tags; ++i)
    {
      int* tag = tagger_create_tag(tagger, (const char*)tag_names[i], tag_sizes[i]);
      if (tag == NULL) // Does the tag already exist?
      {
        size_t size;
        tag = tagger_tag(tagger, (const char*)tag_names[i], &size);
        ASSERT(size == tag_sizes[i]);
      }
      memcpy(tag, &array[j], sizeof(int) * tag_sizes[i]);
      j += tag_sizes[i];
    }
    ASSERT(j == var->nvalues);
    DBFreeCompoundarray(var);
  }
}

void silo_file_write_exchanger(silo_file_t* file, const char* exchanger_name, exchanger_t* ex)
{
  START_FUNCTION_TIMER();

  // Collapse the exchanger into a set of integers.
  // Format is: [nprocs rank send_map receive_map]
  // where send_map and receive_map are sets of integers encoding mappings
  // of processes to sets of indices. Such a map has the format
  // [num_procs [process num_indices i0 i1 ... iN] ... ]
  int_array_t* array = int_array_new();

  int pos = 0, proc, *indices, num_indices;
  int_array_append(array, exchanger_num_sends(ex));
  while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
  {
    int_array_append(array, proc);
    int_array_append(array, num_indices);
    for (int i = 0; i < num_indices; ++i)
      int_array_append(array, indices[i]);
  }
  pos = 0;
  int_array_append(array, exchanger_num_receives(ex));
  while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
  {
    int_array_append(array, proc);
    int_array_append(array, num_indices);
    for (int i = 0; i < num_indices; ++i)
      int_array_append(array, indices[i]);
  }

  // Write the exchanger array to the file.
  silo_file_write_int_array(file, exchanger_name, array->data, array->size);

  // Clean up.
  int_array_free(array);
  STOP_FUNCTION_TIMER();
}

exchanger_t* silo_file_read_exchanger(silo_file_t* file, const char* exchanger_name, MPI_Comm comm)
{
  START_FUNCTION_TIMER();

  // Read the exchanger array in.
  size_t size;
  int* array = silo_file_read_int_array(file, exchanger_name, &size);

  // Create the exchanger.
  exchanger_t* ex = exchanger_new(comm);
  int i = 0;
  int num_sends = array[i++];
  for (int j = 0; j < num_sends; ++j)
  {
    int proc = array[i++];
    int num_indices = array[i++];
    exchanger_set_send(ex, proc, &array[i], num_indices, true);
    i += num_indices;
  }
  int num_receives = array[i++];
  for (int j = 0; j < num_receives; ++j)
  {
    int proc = array[i++];
    int num_indices = array[i++];
    exchanger_set_receive(ex, proc, &array[i], num_indices, true);
    i += num_indices;
  }
  ASSERT(i == size);

  // Clean up.
  polymec_free(array);

  STOP_FUNCTION_TIMER();
  return ex;
}

void silo_file_write_polymesh(silo_file_t* file,
                              const char* mesh_name,
                              polymesh_t* mesh)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_CLOBBER);

  silo_file_push_domain_dir(file);

  // This is optional for now, but we'll give it anyway.
  char *coordnames[3];
  coordnames[0] = (char*)"xcoords";
  coordnames[1] = (char*)"ycoords";
  coordnames[2] = (char*)"zcoords";

  // Node coordinates.
  int num_nodes = mesh->num_nodes;
  real_t* x = polymec_malloc(sizeof(real_t) * num_nodes);
  real_t* y = polymec_malloc(sizeof(real_t) * num_nodes);
  real_t* z = polymec_malloc(sizeof(real_t) * num_nodes);
  for (int i = 0; i < num_nodes; ++i)
  {
    x[i] = mesh->nodes[i].x;
    y[i] = mesh->nodes[i].y;
    z[i] = mesh->nodes[i].z;
  }
  real_t* coords[3];
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;

  // The polyhedral zone list is referred to in the options list.
  DBoptlist* optlist = DBMakeOptlist(10);
  char zonelist_name[FILENAME_MAX+1];
  snprintf(zonelist_name, FILENAME_MAX, "%s_zonelist", mesh_name);
  DBAddOption(optlist, DBOPT_PHZONELIST, zonelist_name);

  // Stick in step/time information if needed.
  if (file->step >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &file->step);
  if (reals_equal(file->time, -REAL_MAX))
  {
    double t = (double)file->time;
    DBAddOption(optlist, DBOPT_DTIME, &t);
  }

  // Write out the 3D polyhedral mesh.
  int num_cells = mesh->num_cells;
  int result = DBPutUcdmesh(file->dbfile, (char*)mesh_name, 3, 
                            (char const* const*)coordnames, coords,
                            num_nodes, num_cells, 0, 0, SILO_FLOAT_TYPE, optlist);
  if (result == -1)
    polymec_error("silo_file_write_polymesh: Could not write mesh '%s'.", mesh_name);

  // Partial cleanup.
  polymec_free(x);
  polymec_free(y);
  polymec_free(z);

  // Construct the silo face-node info.  We rely on the mesh having
  // the faces nodes arranged counter-clockwise around the face.
  int num_faces = mesh->num_faces;
  int* face_node_counts = polymec_malloc(sizeof(int) * num_faces);
  char* ext_faces = polymec_malloc(sizeof(char) * num_faces);
  for (int i = 0; i < num_faces; ++i)
  {
    face_node_counts[i] = mesh->face_node_offsets[i+1] - mesh->face_node_offsets[i];
    if (mesh->face_cells[2*i+1] == -1)
      ext_faces[i] = 0x1;
    else
      ext_faces[i] = 0x0;
  }

  // Construct the silo cell-face info.  Silo uses the same 1's complement
  // convention we use for indicating face orientation, so we can
  // simply copy our faces.
  int* cell_face_counts = polymec_calloc(num_cells + mesh->num_ghost_cells,
                                         sizeof(int));
  for (int i = 0; i < num_cells; ++i)
    cell_face_counts[i] = mesh->cell_face_offsets[i+1] - mesh->cell_face_offsets[i];

  // Write the connectivity information.
  result = DBPutPHZonelist(file->dbfile, zonelist_name, num_faces, face_node_counts,
                           mesh->face_node_offsets[num_faces], mesh->face_nodes,
                           ext_faces, num_cells + mesh->num_ghost_cells, cell_face_counts,
                           mesh->cell_face_offsets[num_cells], mesh->cell_faces,
                           0, 0, num_cells-1, optlist);
  if (result == -1)
    polymec_error("silo_file_write_polymesh: Could not write connectivity data for mesh '%s'.", mesh_name);

  // Partial cleanup.
  polymec_free(face_node_counts);
  polymec_free(ext_faces);
  polymec_free(cell_face_counts);

  // Finally, write out the face_cells array.
  {
    char name[FILENAME_MAX+1];
    snprintf(name, FILENAME_MAX, "%s_face_cells", mesh_name);
    silo_file_write_int_array(file, name, mesh->face_cells, 2*mesh->num_faces);
  }

  // Write out tag information.
  {
    char tag_name[FILENAME_MAX+1];
    snprintf(tag_name, FILENAME_MAX, "%s_node_tags", mesh_name);
    silo_file_write_tags(file, mesh->node_tags, tag_name);
    snprintf(tag_name, FILENAME_MAX, "%s_edge_tags", mesh_name);
    silo_file_write_tags(file, mesh->edge_tags, tag_name);
    snprintf(tag_name, FILENAME_MAX, "%s_face_tags", mesh_name);
    silo_file_write_tags(file, mesh->face_tags, tag_name);
    snprintf(tag_name, FILENAME_MAX, "%s_cell_tags", mesh_name);
    silo_file_write_tags(file, mesh->cell_tags, tag_name);
  }

  // Write out exchanger information. 
  {
    char ex_name[FILENAME_MAX+1];
    snprintf(ex_name, FILENAME_MAX, "%s_cell_exchanger", mesh_name);
    silo_file_write_exchanger(file, ex_name, polymesh_exchanger(mesh, POLYMESH_CELL));
  }

  // Write out the number of mesh cells/faces/nodes/edges to special variables.
  char num_cells_var[FILENAME_MAX+1], num_faces_var[FILENAME_MAX+1],
       num_edges_var[FILENAME_MAX+1], num_nodes_var[FILENAME_MAX+1];
  snprintf(num_cells_var, FILENAME_MAX, "%s_mesh_num_cells", mesh_name);
  snprintf(num_faces_var, FILENAME_MAX, "%s_mesh_num_faces", mesh_name);
  snprintf(num_edges_var, FILENAME_MAX, "%s_mesh_num_edges", mesh_name);
  snprintf(num_nodes_var, FILENAME_MAX, "%s_mesh_num_nodes", mesh_name);
  int one = 1;
  DBWrite(file->dbfile, num_cells_var, &mesh->num_cells, &one, 1, DB_INT);
  DBWrite(file->dbfile, num_faces_var, &mesh->num_faces, &one, 1, DB_INT);
  DBWrite(file->dbfile, num_edges_var, &mesh->num_edges, &one, 1, DB_INT);
  DBWrite(file->dbfile, num_nodes_var, &mesh->num_nodes, &one, 1, DB_INT);
  
  // Clean up.
  DBFreeOptlist(optlist);

#if POLYMEC_HAVE_MPI
  // For parallel environments, add a subdomain entry.
  if (file->nproc > 1)
    silo_file_add_subdomain_mesh(file, mesh_name, DB_UCDMESH, NULL);
#endif

  silo_file_pop_dir(file);

  STOP_FUNCTION_TIMER();
}

extern void polymesh_set_exchanger(polymesh_t* mesh, 
                                   polymesh_centering_t centering,
                                   exchanger_t* exchanger);
polymesh_t* silo_file_read_polymesh(silo_file_t* file,
                                    const char* mesh_name)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_READ);

  silo_file_push_domain_dir(file);

  DBucdmesh* ucd_mesh = DBGetUcdmesh(file->dbfile, mesh_name);
  if (ucd_mesh == NULL)
  {
    log_urgent("No mesh named '%s' was found within the Silo file.", mesh_name);
    silo_file_pop_dir(file);
    STOP_FUNCTION_TIMER();
    return NULL;
  }
  ASSERT(ucd_mesh->ndims == 3);

  // Also get the polyhedral zone list.
  char phzl_name[FILENAME_MAX+1];
  snprintf(phzl_name, FILENAME_MAX, "%s_zonelist", mesh_name);
  DBphzonelist* ph_zonelist = DBGetPHZonelist(file->dbfile, phzl_name);
  if (ph_zonelist == NULL)
  {
    log_urgent("Mesh '%s' is not a polymec polyhedral mesh.", mesh_name);
    silo_file_pop_dir(file);
    STOP_FUNCTION_TIMER();
    return NULL;
  }

  // Decipher the mesh object.
  int num_cells = ph_zonelist->hi_offset + 1;
  int num_ghost_cells = ph_zonelist->nzones - num_cells;
  int num_faces = ph_zonelist->nfaces;
  int num_nodes = ucd_mesh->nnodes;
#if POLYMEC_HAVE_MPI
  MPI_Comm comm = file->comm;
#else
  MPI_Comm comm = MPI_COMM_WORLD;
#endif
  polymesh_t* mesh = polymesh_new(comm, num_cells, num_ghost_cells, 
                                  num_faces, num_nodes);

  // Set node positions.
  real_t* x = ucd_mesh->coords[0];
  real_t* y = ucd_mesh->coords[1];
  real_t* z = ucd_mesh->coords[2];
  for (int n = 0; n < num_nodes; ++n)
  {
    mesh->nodes[n].x = x[n];
    mesh->nodes[n].y = y[n];
    mesh->nodes[n].z = z[n];
  }

  // Set up cell face counts and face node counts.
  mesh->cell_face_offsets[0] = 0;
  for (int c = 0; c < num_cells; ++c)
    mesh->cell_face_offsets[c+1] = mesh->cell_face_offsets[c] + ph_zonelist->facecnt[c];
  mesh->face_node_offsets[0] = 0;
  for (int f = 0; f < num_faces; ++f)
    mesh->face_node_offsets[f+1] = mesh->face_node_offsets[f] + ph_zonelist->nodecnt[f];
  polymesh_reserve_connectivity_storage(mesh);

  // Read in the face_cells array.
  {
    char name[FILENAME_MAX+1];
    snprintf(name, FILENAME_MAX, "%s_face_cells", mesh_name);
    size_t num_face_cells;
    int* face_cells = silo_file_read_int_array(file, name, &num_face_cells);
    ASSERT(num_face_cells == 2*mesh->num_faces);
    memcpy(mesh->face_cells, face_cells, sizeof(int) * 2 * mesh->num_faces);
    polymec_free(face_cells);
  }

  // Fill in cell faces and face nodes.
  memcpy(mesh->cell_faces, ph_zonelist->facelist, sizeof(int) * mesh->cell_face_offsets[mesh->num_cells]);
  memcpy(mesh->face_nodes, ph_zonelist->nodelist, sizeof(int) * mesh->face_node_offsets[mesh->num_faces]);

  // Finish constructing the mesh.
  polymesh_construct_edges(mesh);
  polymesh_compute_geometry(mesh);

  // Read in tag information.
  {
    char tag_name[FILENAME_MAX+1];
    snprintf(tag_name, FILENAME_MAX, "%s_node_tags", mesh_name);
    silo_file_read_tags(file, tag_name, mesh->node_tags);
    snprintf(tag_name, FILENAME_MAX, "%s_edge_tags", mesh_name);
    silo_file_read_tags(file, tag_name, mesh->edge_tags);
    snprintf(tag_name, FILENAME_MAX, "%s_face_tags", mesh_name);
    silo_file_read_tags(file, tag_name, mesh->face_tags);
    snprintf(tag_name, FILENAME_MAX, "%s_cell_tags", mesh_name);
    silo_file_read_tags(file, tag_name, mesh->cell_tags);
  }

  // Read in exchanger information.
  {
    char ex_name[FILENAME_MAX+1];
    snprintf(ex_name, FILENAME_MAX, "%s_cell_exchanger", mesh_name);
    polymesh_set_exchanger(mesh, POLYMESH_CELL, silo_file_read_exchanger(file, ex_name, mesh->comm));
  }

  // Clean up.
  DBFreeUcdmesh(ucd_mesh);
  DBFreePHZonelist(ph_zonelist);

  silo_file_pop_dir(file);

  STOP_FUNCTION_TIMER();
  return mesh;
}

bool silo_file_contains_polymesh(silo_file_t* file, const char* mesh_name)
{
  bool result = false;
  silo_file_push_domain_dir(file);
  result = (DBInqVarExists(file->dbfile, mesh_name) && 
            (DBInqVarType(file->dbfile, mesh_name) == DB_UCDMESH));
  silo_file_pop_dir(file);
  return result;
}

static void silo_file_write_polymesh_field_comp(silo_file_t* file,
                                                const char* field_name,
                                                const char* mesh_name,
                                                int component,
                                                polymesh_centering_t centering,
                                                real_t* field_data,
                                                field_metadata_t* md)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_CLOBBER);

  silo_file_push_domain_dir(file);

  // How many elements does our mesh have?
  char num_elems_var[FILENAME_MAX+1];
  int cent;
  switch(centering)
  {
    case POLYMESH_CELL: cent = DB_ZONECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_cells", mesh_name); break;
    case POLYMESH_FACE: cent = DB_FACECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_faces", mesh_name); break;
    case POLYMESH_EDGE: cent = DB_EDGECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_edges", mesh_name); break;
    case POLYMESH_NODE: cent = DB_NODECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_nodes", mesh_name);
  }
  ASSERT(DBInqVarExists(file->dbfile, num_elems_var));
  int num_elems;
  DBReadVar(file->dbfile, num_elems_var, &num_elems);

  // Create an optlist from our field metadata.
  DBoptlist* optlist = optlist_from_metadata(md, component);

  char var_name[FILENAME_MAX+1];
  snprintf(var_name, FILENAME_MAX, "%s_%d_%s", field_name, component, mesh_name);
  DBPutUcdvar1(file->dbfile, var_name, mesh_name, field_data, num_elems, NULL, 0, SILO_FLOAT_TYPE, cent, optlist);

#if POLYMEC_HAVE_MPI
  // Add a subdomain entry for parallel environments.
  if (file->nproc > 1)
    silo_file_add_subdomain_field(file, mesh_name, field_name, DB_UCDVAR, optlist);
#endif

  optlist_free(optlist);
  silo_file_pop_dir(file);

  STOP_FUNCTION_TIMER();
}

void silo_file_write_polymesh_field(silo_file_t* file,
                                    const char* field_name,
                                    const char* mesh_name,
                                    polymesh_field_t* field)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_CLOBBER);

  // How many elements does our mesh have?
  silo_file_push_domain_dir(file);
  char num_elems_var[FILENAME_MAX+1];
  int cent;
  switch(field->centering)
  {
    case POLYMESH_CELL: cent = DB_ZONECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_cells", mesh_name); break;
    case POLYMESH_FACE: cent = DB_FACECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_faces", mesh_name); break;
    case POLYMESH_EDGE: cent = DB_EDGECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_edges", mesh_name); break;
    case POLYMESH_NODE: cent = DB_NODECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_nodes", mesh_name);
  }
  ASSERT(DBInqVarExists(file->dbfile, num_elems_var));
  int num_elems;
  DBReadVar(file->dbfile, num_elems_var, &num_elems);
  silo_file_pop_dir(file);

  // Write the field metadata.
  field_metadata_t* md = polymesh_field_metadata(field);
  char md_name[FILENAME_MAX+1];
  snprintf(md_name, FILENAME_MAX, "%s_%s_md", field_name, mesh_name);
  silo_file_write_field_metadata(file, md_name, md);

  DECLARE_POLYMESH_FIELD_ARRAY(field_data, field);
  real_t* comp_data = polymec_malloc(sizeof(real_t) * num_elems);
  for (int c = 0; c < field->num_components; ++c)
  {
    for (int i = 0; i < num_elems; ++i)
      comp_data[i] = field_data[i][c];
    silo_file_write_polymesh_field_comp(file, field_name, mesh_name, c, 
                                        field->centering, comp_data, md);
  }
  polymec_free(comp_data);
  STOP_FUNCTION_TIMER();
}

static bool silo_file_read_polymesh_field_comp(silo_file_t* file,
                                               const char* field_name,
                                               const char* mesh_name,
                                               int component,
                                               polymesh_centering_t centering,
                                               real_t* field_data)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_READ);

  silo_file_push_domain_dir(file);

  // How many elements does our mesh have?
  char num_elems_var[FILENAME_MAX+1];
  int cent = DB_ZONECENT;
  switch(centering)
  {
    case POLYMESH_CELL: cent = DB_ZONECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_cells", mesh_name); break;
    case POLYMESH_FACE: cent = DB_FACECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_faces", mesh_name); break;
    case POLYMESH_EDGE: cent = DB_EDGECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_edges", mesh_name); break;
    case POLYMESH_NODE: cent = DB_NODECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_nodes", mesh_name);
  }
  ASSERT(DBInqVarExists(file->dbfile, num_elems_var));
  int num_elems;
  DBReadVar(file->dbfile, num_elems_var, &num_elems);

  char var_name[FILENAME_MAX+1];
  snprintf(var_name, FILENAME_MAX, "%s_%d_%s", field_name, component, mesh_name);
  DBucdvar* var = DBGetUcdvar(file->dbfile, var_name);
  if (var == NULL)
  {
    log_urgent("Component %d of field '%s' was not found in the Silo file.", component, field_name);
    silo_file_pop_dir(file);
    STOP_FUNCTION_TIMER();
    return false;
  }
  if (var->centering != cent)
  {
    log_urgent("Field '%s' has the incorrect centering.", field_name);
    silo_file_pop_dir(file);
    STOP_FUNCTION_TIMER();
    return false;
  }
  memcpy(field_data, var->vals[0], sizeof(real_t) * num_elems);
  DBFreeUcdvar(var);

  silo_file_pop_dir(file);
  STOP_FUNCTION_TIMER();
  return true;
}

void silo_file_read_polymesh_field(silo_file_t* file,
                                   const char* field_name,
                                   const char* mesh_name,
                                   polymesh_field_t* field)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_READ);

  silo_file_push_domain_dir(file);

  // How many elements does our mesh have?
  char num_elems_var[FILENAME_MAX+1];
  int cent = DB_ZONECENT;
  switch(field->centering)
  {
    case POLYMESH_CELL: cent = DB_ZONECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_cells", mesh_name); break;
    case POLYMESH_FACE: cent = DB_FACECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_faces", mesh_name); break;
    case POLYMESH_EDGE: cent = DB_EDGECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_edges", mesh_name); break;
    case POLYMESH_NODE: cent = DB_NODECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_nodes", mesh_name);
  }
  ASSERT(DBInqVarExists(file->dbfile, num_elems_var));
  int num_elems;
  DBReadVar(file->dbfile, num_elems_var, &num_elems);

  silo_file_pop_dir(file);

  // Read the field metadata.
  field_metadata_t* md = polymesh_field_metadata(field);
  char md_name[FILENAME_MAX+1];
  snprintf(md_name, FILENAME_MAX, "%s_%s_md", field_name, mesh_name);
  silo_file_read_field_metadata(file, md_name, md);

  DECLARE_POLYMESH_FIELD_ARRAY(field_data, field);
  real_t* comp_data = polymec_malloc(sizeof(real_t) * num_elems);
  for (int c = 0; c < field->num_components; ++c)
  {
    silo_file_read_polymesh_field_comp(file, field_name, mesh_name, c, field->centering, comp_data);
    for (int i = 0; i < num_elems; ++i)
      field_data[i][c] = comp_data[i];
  }
  polymec_free(comp_data);
  STOP_FUNCTION_TIMER();
}

bool silo_file_contains_polymesh_field(silo_file_t* file, 
                                       const char* field_name, 
                                       const char* mesh_name,
                                       polymesh_centering_t centering)
{
  START_FUNCTION_TIMER();
  bool result = false;
  if (silo_file_contains_polymesh(file, mesh_name)) // mesh exists...
  {
    // Look for the field's metadata array.
    char md_name[FILENAME_MAX+1];
    snprintf(md_name, FILENAME_MAX, "%s_%s_md", field_name, mesh_name);
    result = silo_file_contains_field_metadata(file, md_name);
  }

  STOP_FUNCTION_TIMER();
  return result;
}

void silo_file_write_point_cloud(silo_file_t* file,
                                 const char* cloud_name,
                                 point_cloud_t* cloud)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_CLOBBER);

  silo_file_push_domain_dir(file);

  // Point coordinates.
  size_t num_points = cloud->num_points;
  point_t* points = cloud->points;
  real_t* x = polymec_malloc(sizeof(real_t) * num_points);
  real_t* y = polymec_malloc(sizeof(real_t) * num_points);
  real_t* z = polymec_malloc(sizeof(real_t) * num_points);
  for (size_t i = 0; i < num_points; ++i)
  {
    x[i] = points[i].x;
    y[i] = points[i].y;
    z[i] = points[i].z;
  }
  real_t* coords[3];
  coords[0] = &(x[0]);
  coords[1] = &(y[0]);
  coords[2] = &(z[0]);

  // Write out the point mesh.
  DBPutPointmesh(file->dbfile, (char*)cloud_name, 3, coords, (int)num_points, SILO_FLOAT_TYPE, NULL); 
  polymec_free(x);
  polymec_free(y);
  polymec_free(z);

  // Write out the number of points to a special variable.
  char num_points_var[FILENAME_MAX+1];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", cloud_name);
  int one = 1;
  DBWrite(file->dbfile, num_points_var, &num_points, &one, 1, DB_INT);
  
  // Write out tag information.
  {
    char tag_name[FILENAME_MAX+1];
    snprintf(tag_name, FILENAME_MAX, "%s_node_tags", cloud_name);
    silo_file_write_tags(file, cloud->tags, tag_name);
  }

#if POLYMEC_HAVE_MPI
  if (file->nproc > 1)
    silo_file_add_subdomain_mesh(file, cloud_name, DB_POINTMESH, NULL);
#endif

  silo_file_pop_dir(file);

  STOP_FUNCTION_TIMER();
}

point_cloud_t* silo_file_read_point_cloud(silo_file_t* file,
                                          const char* cloud_name)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_READ);

  silo_file_push_domain_dir(file);

  // How many points does our cloud have?
  int num_points;
  char num_points_var[FILENAME_MAX+1];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", cloud_name);
  ASSERT(DBInqVarExists(file->dbfile, num_points_var));
  DBReadVar(file->dbfile, num_points_var, &num_points);

  DBpointmesh* pm = DBGetPointmesh(file->dbfile, (char*)cloud_name);
  if (pm == NULL)
  {
    log_urgent("Point mesh '%s' was not found in the Silo file.", cloud_name);
    silo_file_pop_dir(file);
    STOP_FUNCTION_TIMER();
    return NULL;
  }
  ASSERT(num_points == pm->nels);
  point_t* points = polymec_malloc(sizeof(point_t) * num_points);
  real_t* x = pm->coords[0];
  real_t* y = pm->coords[1];
  real_t* z = pm->coords[2];
  for (int p = 0; p < num_points; ++p)
  {
    points[p].x = x[p];
    points[p].y = y[p];
    points[p].z = z[p];
  }
#if POLYMEC_HAVE_MPI
  point_cloud_t* cloud = point_cloud_from_points(file->comm, points, num_points);
#else
  point_cloud_t* cloud = point_cloud_from_points(MPI_COMM_WORLD, points, num_points);
#endif
  DBFreePointmesh(pm);
  polymec_free(points);

  // Read in tag information.
  {
    char tag_name[FILENAME_MAX+1];
    snprintf(tag_name, FILENAME_MAX, "%s_node_tags", cloud_name);
    silo_file_read_tags(file, tag_name, cloud->tags);
  }

  silo_file_pop_dir(file);

  STOP_FUNCTION_TIMER();
  return cloud;
}

void silo_file_write_planar_polymesh(silo_file_t* file, 
                                     const char* mesh_name,
                                     planar_polymesh_t* mesh)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_CLOBBER);

  silo_file_push_domain_dir(file);

  // This is optional for now, but we'll give it anyway.
  char *coordnames[2];
  coordnames[0] = (char*)"xcoords";
  coordnames[1] = (char*)"ycoords";

  // Node coordinates.
  int num_nodes = mesh->num_nodes;
  real_t* x = polymec_malloc(sizeof(real_t) * num_nodes);
  real_t* y = polymec_malloc(sizeof(real_t) * num_nodes);
  for (int i = 0; i < num_nodes; ++i)
  {
    x[i] = mesh->nodes[i].x;
    y[i] = mesh->nodes[i].y;
  }
  real_t* coords[2];
  coords[0] = x;
  coords[1] = y;

  // Stick in step/time information if needed.
  DBoptlist* optlist = DBMakeOptlist(10);
  if (file->step >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &file->step);
  if (reals_equal(file->time, -REAL_MAX))
  {
    double t = (double)file->time;
    DBAddOption(optlist, DBOPT_DTIME, &t);
  }

  // Write cells to the file.
  char zonelist_name[FILENAME_MAX+1];
  snprintf(zonelist_name, FILENAME_MAX, "%s_zonelist", mesh_name);
  {
    int num_cell_edges[mesh->num_cells];
    int total_num_cell_nodes = mesh->cell_edge_offsets[mesh->num_cells];
    int cell_nodes[total_num_cell_nodes];
    int shapes[mesh->num_cells], ones[mesh->num_cells];
    int l = 0;
    for (int c = 0; c < mesh->num_cells; ++c)
    {
      shapes[c] = DB_ZONETYPE_POLYGON;
      ones[c] = 1;
      num_cell_edges[c] = mesh->cell_edge_offsets[c+1] - mesh->cell_edge_offsets[c];
      for (int n = 0; n < num_cell_edges[c]; ++n, ++l)
      {
        int edge = mesh->cell_edges[l];
        if (edge < 0)
          cell_nodes[l] = mesh->edge_nodes[2*(~edge)+1];
        else
          cell_nodes[l] = mesh->edge_nodes[2*edge];
      }
    }
    DBPutZonelist2(file->dbfile, zonelist_name, mesh->num_cells, 2, 
                   cell_nodes, total_num_cell_nodes, 0, 0, 0,
                   shapes, num_cell_edges, ones, mesh->num_cells, NULL);
  }

  // Write edges (2D faces) to the file.
  int num_edge_nodes = 2;
  char facelist_name[FILENAME_MAX+1];
  snprintf(facelist_name, FILENAME_MAX, "%s_facelist", mesh_name);
  DBPutFacelist(file->dbfile, facelist_name, mesh->num_edges, 2, mesh->edge_nodes, 
                2*mesh->num_edges, 0, NULL, &num_edge_nodes, &(mesh->num_edges), 1, 
                NULL, NULL, 0);

  // Write out the (2D) planar poly(gonal) mesh.
  int num_cells = mesh->num_cells;
  int result = DBPutUcdmesh(file->dbfile, (char*)mesh_name, 2, 
                            (char const* const*)coordnames, coords,
                            num_nodes, num_cells, zonelist_name, facelist_name, SILO_FLOAT_TYPE, optlist);
  if (result == -1)
    polymec_error("silo_file_write_planar_polymesh: Could not write mesh '%s'.", mesh_name);

  // Partial cleanup.
  polymec_free(x);
  polymec_free(y);

  // Finally, write out the cell_edges, edge_cells, edge_nodes arrays.
  {
    char ce_name[FILENAME_MAX+1];
    snprintf(ce_name, FILENAME_MAX, "%s_cell_edges", mesh_name);
    silo_file_write_int_array(file, ce_name, mesh->cell_edges, mesh->cell_edge_offsets[mesh->num_cells]);

    char ec_name[FILENAME_MAX+1];
    snprintf(ec_name, FILENAME_MAX, "%s_edge_cells", mesh_name);
    silo_file_write_int_array(file, ec_name, mesh->edge_cells, 2*mesh->num_edges);

    char en_name[FILENAME_MAX+1];
    snprintf(en_name, FILENAME_MAX, "%s_edge_nodes", mesh_name);
    silo_file_write_int_array(file, en_name, mesh->edge_nodes, 2*mesh->num_edges);
  }

  // Write out tag information.
  {
    char tag_name[FILENAME_MAX+1];
    snprintf(tag_name, FILENAME_MAX, "%s_node_tags", mesh_name);
    silo_file_write_tags(file, mesh->node_tags, tag_name);
    snprintf(tag_name, FILENAME_MAX, "%s_edge_tags", mesh_name);
    silo_file_write_tags(file, mesh->edge_tags, tag_name);
    snprintf(tag_name, FILENAME_MAX, "%s_cell_tags", mesh_name);
    silo_file_write_tags(file, mesh->cell_tags, tag_name);
  }

  // Write out the number of mesh cells/faces/nodes/edges to special variables.
  char num_cells_var[FILENAME_MAX+1], num_edges_var[FILENAME_MAX+1], 
       num_nodes_var[FILENAME_MAX+1];
  snprintf(num_cells_var, FILENAME_MAX, "%s_mesh_num_cells", mesh_name);
  snprintf(num_edges_var, FILENAME_MAX, "%s_mesh_num_edges", mesh_name);
  snprintf(num_nodes_var, FILENAME_MAX, "%s_mesh_num_nodes", mesh_name);
  int one = 1;
  DBWrite(file->dbfile, num_cells_var, &mesh->num_cells, &one, 1, DB_INT);
  DBWrite(file->dbfile, num_edges_var, &mesh->num_edges, &one, 1, DB_INT);
  DBWrite(file->dbfile, num_nodes_var, &mesh->num_nodes, &one, 1, DB_INT);
  
  // Clean up.
  DBFreeOptlist(optlist);

  silo_file_pop_dir(file);

  STOP_FUNCTION_TIMER();
}

planar_polymesh_t* silo_file_read_planar_polymesh(silo_file_t* file, 
                                                  const char* mesh_name)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_READ);

  silo_file_push_domain_dir(file);

  DBucdmesh* ucd_mesh = DBGetUcdmesh(file->dbfile, mesh_name);
  if (ucd_mesh == NULL)
  {
    log_urgent("No mesh named '%s' was found within the Silo file.", mesh_name);
    silo_file_pop_dir(file);
    return NULL;
  }
  if (ucd_mesh->ndims != 2)
  {
    log_urgent("Mesh '%s' is not a polymec planar polygonal mesh.", mesh_name);
    silo_file_pop_dir(file);
    return NULL;
  }

  // Decipher the mesh object.
  int num_cells = ucd_mesh->zones->nzones;
  int num_edges = ucd_mesh->faces->nfaces;
  int num_nodes = ucd_mesh->nnodes;
  planar_polymesh_t* mesh = planar_polymesh_new(num_cells, num_edges, num_nodes);

  // Set node positions.
  real_t* x = ucd_mesh->coords[0];
  real_t* y = ucd_mesh->coords[1];
  for (int n = 0; n < num_nodes; ++n)
  {
    mesh->nodes[n].x = x[n];
    mesh->nodes[n].y = y[n];
  }

  // Set up cell face counts and face node counts.
  mesh->cell_edge_offsets[0] = 0;
  for (int c = 0; c < num_cells; ++c)
    mesh->cell_edge_offsets[c+1] = mesh->cell_edge_offsets[c] + ucd_mesh->zones->shapesize[c];
  planar_polymesh_reserve_connectivity_storage(mesh);

  // Read in the cell_edges array.
  {
    char name[FILENAME_MAX+1];
    snprintf(name, FILENAME_MAX, "%s_cell_edges", mesh_name);
    size_t num_cell_edges;
    int* cell_edges = silo_file_read_int_array(file, name, &num_cell_edges);
    ASSERT(num_cell_edges == mesh->cell_edge_offsets[mesh->num_cells]);
    memcpy(mesh->cell_edges, cell_edges, sizeof(int) * num_cell_edges);
    polymec_free(cell_edges);
  }

  // Read in the edge_cells array.
  {
    char name[FILENAME_MAX+1];
    snprintf(name, FILENAME_MAX, "%s_edge_cells", mesh_name);
    size_t num_edge_cells;
    int* edge_cells = silo_file_read_int_array(file, name, &num_edge_cells);
    ASSERT(num_edge_cells == 2*mesh->num_edges);
    memcpy(mesh->edge_cells, edge_cells, sizeof(int) * 2 * mesh->num_edges);
    polymec_free(edge_cells);
  }

  // Read in the edge_nodes array.
  {
    char name[FILENAME_MAX+1];
    snprintf(name, FILENAME_MAX, "%s_edge_nodes", mesh_name);
    size_t num_edge_nodes;
    int* edge_nodes = silo_file_read_int_array(file, name, &num_edge_nodes);
    ASSERT(num_edge_nodes == 2*mesh->num_edges);
    memcpy(mesh->edge_nodes, edge_nodes, sizeof(int) * 2 * mesh->num_edges);
    polymec_free(edge_nodes);
  }

  // Read in tag information.
  {
    char tag_name[FILENAME_MAX+1];
    snprintf(tag_name, FILENAME_MAX, "%s_node_tags", mesh_name);
    silo_file_read_tags(file, tag_name, mesh->node_tags);
    snprintf(tag_name, FILENAME_MAX, "%s_edge_tags", mesh_name);
    silo_file_read_tags(file, tag_name, mesh->edge_tags);
    snprintf(tag_name, FILENAME_MAX, "%s_cell_tags", mesh_name);
    silo_file_read_tags(file, tag_name, mesh->cell_tags);
  }

  // Clean up.
  DBFreeUcdmesh(ucd_mesh);
  silo_file_pop_dir(file);

  STOP_FUNCTION_TIMER();
  return mesh;
}

bool silo_file_contains_planar_polymesh(silo_file_t* file, 
                                        const char* mesh_name)
{
  bool result = false;
  silo_file_push_domain_dir(file);
  result = (DBInqVarExists(file->dbfile, mesh_name) && 
            (DBInqVarType(file->dbfile, mesh_name) == DB_UCDMESH));
  silo_file_pop_dir(file);
  return result;
}

bool silo_file_contains_point_cloud(silo_file_t* file, const char* cloud_name)
{
  bool result = false;
  silo_file_push_domain_dir(file);
  result = (DBInqVarExists(file->dbfile, cloud_name) && 
            (DBInqVarType(file->dbfile, cloud_name) == DB_POINTMESH));
  silo_file_pop_dir(file);
  return result;
}

static void silo_file_write_point_field_comp(silo_file_t* file,
                                             const char* field_name,
                                             const char* cloud_name,
                                             int component,
                                             real_t* field_data,
                                             field_metadata_t* md)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_CLOBBER);

  silo_file_push_domain_dir(file);

  // How many points does our mesh have?
  char num_points_var[FILENAME_MAX+1];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", cloud_name);
  ASSERT(DBInqVarExists(file->dbfile, num_points_var));
  int num_points;
  DBReadVar(file->dbfile, num_points_var, &num_points);

  DBoptlist* optlist = optlist_from_metadata(md, component);

  // Write the field.
  char var_name[FILENAME_MAX+1];
  snprintf(var_name, FILENAME_MAX, "%s_%d_%s", field_name, component, cloud_name);
  DBPutPointvar1(file->dbfile, var_name, cloud_name, field_data, num_points, SILO_FLOAT_TYPE, optlist);

#if POLYMEC_HAVE_MPI
  if (file->nproc > 1)
    silo_file_add_subdomain_field(file, cloud_name, field_name, DB_POINTVAR, optlist);
#endif
  optlist_free(optlist);

  silo_file_pop_dir(file);
  STOP_FUNCTION_TIMER();
}

static bool silo_file_read_point_field_comp(silo_file_t* file,
                                            const char* field_name,
                                            const char* cloud_name,
                                            int component, 
                                            real_t* field_data)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_READ);

  silo_file_push_domain_dir(file);

  char var_name[FILENAME_MAX+1];
  snprintf(var_name, FILENAME_MAX, "%s_%d_%s", field_name, component, cloud_name);
  DBmeshvar* var = DBGetPointvar(file->dbfile, var_name);
  if (var == NULL)
  {
    log_urgent("Component %d of field '%s' was not found in the Silo file.", 
               component, field_name);
    silo_file_pop_dir(file);
    STOP_FUNCTION_TIMER();
    return false;
  }
  memcpy(field_data, var->vals[0], sizeof(real_t) * var->nels);
  DBFreePointvar(var);
  silo_file_pop_dir(file);
  STOP_FUNCTION_TIMER();
  return true;
}

void silo_file_write_point_field(silo_file_t* file,
                                 const char* field_name,
                                 const char* cloud_name,
                                 point_cloud_field_t* field)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_CLOBBER);

  // How many points does our mesh have?
  silo_file_push_domain_dir(file);

  char num_points_var[FILENAME_MAX+1];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", cloud_name);
  ASSERT(DBInqVarExists(file->dbfile, num_points_var));
  int num_points;
  DBReadVar(file->dbfile, num_points_var, &num_points);

  silo_file_pop_dir(file);

  // Write the field metadata.
  field_metadata_t* md = point_cloud_field_metadata(field);
  char md_name[FILENAME_MAX+1];
  snprintf(md_name, FILENAME_MAX, "%s_%s_md", field_name, cloud_name);
  silo_file_write_field_metadata(file, md_name, md);

  DECLARE_POINT_CLOUD_FIELD_ARRAY(field_data, field);
  for (int c = 0; c < field->num_components; ++c)
  {
    silo_file_write_point_field_comp(file, field_name, cloud_name, c, 
                                     (real_t*)(field_data[c]), md);
  }
  STOP_FUNCTION_TIMER();
}

void silo_file_read_point_field(silo_file_t* file,
                                const char* field_name,
                                const char* cloud_name,
                                point_cloud_field_t* field)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_READ);

  silo_file_push_domain_dir(file);

  // How many points does our mesh have?
  char num_points_var[FILENAME_MAX+1];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", cloud_name);
  ASSERT(DBInqVarExists(file->dbfile, num_points_var));
  int num_points;
  DBReadVar(file->dbfile, num_points_var, &num_points);

  // Read the field metadata.
  field_metadata_t* md = point_cloud_field_metadata(field);
  char md_name[FILENAME_MAX+1];
  snprintf(md_name, FILENAME_MAX, "%s_%s_md", field_name, cloud_name);
  silo_file_read_field_metadata(file, md_name, md);

  DECLARE_POINT_CLOUD_FIELD_ARRAY(field_data, field);
  for (int c = 0; c < field->num_components; ++c)
    silo_file_read_point_field_comp(file, field_name, cloud_name, c, (real_t*)(field_data[c]));

  silo_file_pop_dir(file);

  STOP_FUNCTION_TIMER();
}

bool silo_file_contains_point_field(silo_file_t* file, 
                                    const char* field_name, 
                                    const char* cloud_name)
{
  bool result = false;
  if (silo_file_contains_point_cloud(file, cloud_name))  // point cloud exists...
  {
    // Look for the field's metadata array.
    char md_name[FILENAME_MAX+1];
    snprintf(md_name, FILENAME_MAX, "%s_%s_md", field_name, cloud_name);
    result = silo_file_contains_field_metadata(file, md_name);
  }
  return result;
}


static void expression_dtor(char* key, void* val)
{
  string_free(key);
  silo_expression_t* exp = val;
  string_free(exp->definition);
  polymec_free(exp);
}

void silo_file_write_scalar_expression(silo_file_t* file,
                                       const char* expression_name,
                                       const char* definition)

{
  ASSERT(file->mode == DB_CLOBBER);
  silo_expression_t* exp = polymec_malloc(sizeof(silo_expression_t));
  exp->type = DB_VARTYPE_SCALAR;
  exp->definition = string_dup(definition);
  string_ptr_unordered_map_insert_with_kv_dtor(file->expressions,
                                               string_dup(expression_name), exp,
                                               expression_dtor);
}

void silo_file_write_vector_expression(silo_file_t* file,
                                       const char* expression_name,
                                       const char* definition)
{
  ASSERT(file->mode == DB_CLOBBER);
  silo_expression_t* exp = polymec_malloc(sizeof(silo_expression_t));
  exp->type = DB_VARTYPE_VECTOR;
  exp->definition = string_dup(definition);
  string_ptr_unordered_map_insert_with_kv_dtor(file->expressions,
                                               string_dup(expression_name), exp,
                                               expression_dtor);
}

void silo_file_write_tensor_expression(silo_file_t* file,
                                       const char* expression_name,
                                       const char* definition)
{
  ASSERT(file->mode == DB_CLOBBER);
  silo_expression_t* exp = polymec_malloc(sizeof(silo_expression_t));
  exp->type = DB_VARTYPE_TENSOR;
  exp->definition = string_dup(definition);
  string_ptr_unordered_map_insert_with_kv_dtor(file->expressions,
                                               string_dup(expression_name), exp,
                                               expression_dtor);
}

bool silo_file_contains_string(silo_file_t* file,
                               const char* string_name)
{
  silo_file_push_domain_dir(file);
  char name[FILENAME_MAX+1];
  snprintf(name, FILENAME_MAX, "%s_string", string_name);
  bool result = DBInqVarExists(file->dbfile, name);
  silo_file_pop_dir(file);
  return result;
}

void silo_file_write_string(silo_file_t* file,
                            const char* string_name,
                            char* string_data)
{
  ASSERT(file->mode == DB_CLOBBER);
  ASSERT(string_data != NULL);

  silo_file_push_domain_dir(file);

  int string_size = (int)strlen(string_data);
  if (string_size > 0)
  {
    char string_data_name[FILENAME_MAX+1];
    snprintf(string_data_name, FILENAME_MAX, "%s_string", string_name);
    int result = DBWrite(file->dbfile, string_data_name, string_data, &string_size, 1, DB_CHAR);
    if (result != 0)
      polymec_error("silo_file_write_string: write of string '%s' failed.", string_name);
  }
  silo_file_pop_dir(file);
}

char* silo_file_read_string(silo_file_t* file,
                            const char* string_name)
{
  ASSERT(file->mode == DB_READ);

  silo_file_push_domain_dir(file);

  char string_data_name[FILENAME_MAX+1];
  snprintf(string_data_name, FILENAME_MAX, "%s_string", string_name);
  if (!DBInqVarExists(file->dbfile, string_data_name))
  {
    log_urgent("silo_file_read_string: Could not read string '%s'.", string_name);
    silo_file_pop_dir(file);
    return NULL;
  }
  int string_size = DBGetVarLength(file->dbfile, string_data_name);
  char* string = polymec_malloc(sizeof(char) * (string_size + 1));
  if (string_size > 0)
    DBReadVar(file->dbfile, string_data_name, string);
  string[string_size] = '\0';
  silo_file_pop_dir(file);
  return string;
}

void silo_file_write_real_array(silo_file_t* file,
                                const char* array_name,
                                real_t* array_data,
                                size_t array_size)
{
  ASSERT(file->mode == DB_CLOBBER);
  ASSERT(array_data != NULL);

  silo_file_push_domain_dir(file);

  if (array_size > 0)
  {
    char real_array_name[FILENAME_MAX+1];
    snprintf(real_array_name, FILENAME_MAX, "%s_real_array", array_name);
    int asize = (int)array_size;
    int result = DBWrite(file->dbfile, real_array_name, array_data, &asize, 1, SILO_FLOAT_TYPE);
    if (result != 0)
      polymec_error("silo_file_write_real_array: write of array '%s' failed.", array_name);
  }
  silo_file_pop_dir(file);
}

real_t* silo_file_read_real_array(silo_file_t* file,
                                  const char* array_name,
                                  size_t* array_size)
{
  ASSERT(file->mode == DB_READ);
  ASSERT(array_size != NULL);

  silo_file_push_domain_dir(file);

  char real_array_name[FILENAME_MAX+1];
  snprintf(real_array_name, FILENAME_MAX, "%s_real_array", array_name);
  real_t* result = NULL;
  if (!DBInqVarExists(file->dbfile, real_array_name))
  {
    log_urgent("silo_file_read_real_array: Could not read array '%s'.", array_name);
    silo_file_pop_dir(file);
    return NULL;
  }
  *array_size = (size_t)DBGetVarLength(file->dbfile, real_array_name);
  if (*array_size > 0)
  {
    real_t* array = polymec_malloc(sizeof(real_t) * *array_size);
    DBReadVar(file->dbfile, real_array_name, array);
    result = array;
  }
  silo_file_pop_dir(file);
  return result;
}

void silo_file_write_int_array(silo_file_t* file,
                               const char* array_name,
                               int* array_data,
                               size_t array_size)
{
  ASSERT(file->mode == DB_CLOBBER);
  ASSERT(array_data != NULL);

  silo_file_push_domain_dir(file);

  if (array_size > 0)
  {
    char int_array_name[FILENAME_MAX+1];
    snprintf(int_array_name, FILENAME_MAX, "%s_int_array", array_name);
    int asize = (int)array_size;
    int result = DBWrite(file->dbfile, int_array_name, array_data, &asize, 1, DB_INT);
    if (result != 0)
      polymec_error("silo_file_write_int_array: write of array '%s' failed.", array_name);
  }
  silo_file_pop_dir(file);
}

int* silo_file_read_int_array(silo_file_t* file,
                              const char* array_name,
                              size_t* array_size)
{
  ASSERT(file->mode == DB_READ);
  ASSERT(array_size != NULL);

  silo_file_push_domain_dir(file);

  char int_array_name[FILENAME_MAX+1];
  snprintf(int_array_name, FILENAME_MAX, "%s_int_array", array_name);
  int* result = NULL;
  if (!DBInqVarExists(file->dbfile, int_array_name))
  {
    log_urgent("silo_file_read_int_array: Could not read array '%s'.", array_name);
    silo_file_pop_dir(file);
    return NULL;
  }
  *array_size = (size_t)DBGetVarLength(file->dbfile, int_array_name);
  if (*array_size > 0)
  {
    int* array = polymec_malloc(sizeof(int) * *array_size);
    DBReadVar(file->dbfile, int_array_name, array);
    result = array;
  }
  silo_file_pop_dir(file);
  return result;
}

MPI_Comm silo_file_comm(silo_file_t* file)
{
  return file->comm;
}

DBfile* silo_file_dbfile(silo_file_t* file)
{
  return file->dbfile;
}

void silo_file_add_subdomain_mesh(silo_file_t* file,
                                  const char* mesh_name, 
                                  int silo_mesh_type,
                                  DBoptlist* optlist)
{
#if POLYMEC_HAVE_MPI
  ASSERT(file->mode == DB_CLOBBER);

  if (file->nproc > 1)
  {
    subdomain_mesh_t* mesh = subdomain_mesh_new(mesh_name, silo_mesh_type, optlist);
    ptr_array_append_with_dtor(file->subdomain_meshes, mesh, DTOR(subdomain_mesh_free));
  }
#endif
}

void silo_file_add_subdomain_field(silo_file_t* file,
                                   const char* mesh_name, 
                                   const char* field_name,
                                   int silo_var_type,
                                   DBoptlist* optlist)
{
#if POLYMEC_HAVE_MPI
  ASSERT((file->mode == DB_CLOBBER) || (file->mode == DB_APPEND));

  if (file->nproc > 1)
  {
    subdomain_field_t* field = subdomain_field_new(mesh_name, field_name, silo_var_type, optlist);
    ptr_array_append_with_dtor(file->subdomain_fields, field, DTOR(subdomain_field_free));
  }
#endif
}

void silo_file_push_dir(silo_file_t* file, const char* dir)
{
  if (file->dirs == NULL)
    file->dirs = string_slist_new();

  // Set the new directory and push it onto the stack.
  DBSetDir(file->dbfile, dir);
  string_slist_push_with_dtor(file->dirs, string_dup(dir), string_free);
}

void silo_file_push_domain_dir(silo_file_t* file)
{
#if POLYMEC_HAVE_MPI
  // Determine the domain directory.
  char domain_dir[FILENAME_MAX+1];
  if (file->nproc > 1)
  {
    snprintf(domain_dir, FILENAME_MAX, "/domain_%d", file->rank_in_group);
    DBSetDir(file->dbfile, domain_dir);
  }
  else
    sprintf(domain_dir, "/");
  silo_file_push_dir(file, domain_dir);
#else
  silo_file_push_dir(file, "/");
#endif
}

void silo_file_pop_dir(silo_file_t* file)
{
  ASSERT(!string_slist_empty(file->dirs));

  // Pop the current directory off the stack and delete it.
  char* dir = string_slist_pop(file->dirs, NULL);
  string_free(dir);

  // Set the directory to the front of the stack.
  DBSetDir(file->dbfile, file->dirs->front->value);
}

string_ptr_unordered_map_t* silo_file_scratch(silo_file_t* file)
{
  return file->scratch;
}

bool silo_file_contains_stencil(silo_file_t* file, const char* stencil_name)
{
  char name[FILENAME_MAX+1];
  snprintf(name, FILENAME_MAX, "%s_stencil_name", stencil_name);
  return silo_file_contains_string(file, name);
}

void silo_file_write_stencil(silo_file_t* file,
                             const char* stencil_name,
                             stencil_t* stencil)
{
  START_FUNCTION_TIMER();
  char name_name[FILENAME_MAX+1];
  snprintf(name_name, FILENAME_MAX, "%s_stencil_name", stencil_name);
  silo_file_write_string(file, name_name, stencil->name);

  char offsets_name[FILENAME_MAX+1];
  snprintf(offsets_name, FILENAME_MAX, "%s_stencil_offsets", stencil_name);
  silo_file_write_int_array(file, offsets_name, stencil->offsets, stencil->num_indices+1);

  char indices_name[FILENAME_MAX+1];
  snprintf(indices_name, FILENAME_MAX, "%s_stencil_indices", stencil_name);
  silo_file_write_int_array(file, indices_name, stencil->indices, stencil->offsets[stencil->num_indices]);

  if (stencil->ex != NULL)
  {
    char ex_name[FILENAME_MAX+1];
    snprintf(ex_name, FILENAME_MAX, "%s_stencil_ex", stencil_name);
    silo_file_write_exchanger(file, ex_name, stencil->ex);
  }
  STOP_FUNCTION_TIMER();
}

stencil_t* silo_file_read_stencil(silo_file_t* file,
                                  const char* stencil_name,
                                  MPI_Comm comm)
{
  START_FUNCTION_TIMER();
  stencil_t* s = polymec_malloc(sizeof(stencil_t));
  char name_name[FILENAME_MAX+1];
  snprintf(name_name, FILENAME_MAX, "%s_stencil_name", stencil_name);
  s->name = silo_file_read_string(file, name_name);

  char offsets_name[FILENAME_MAX+1];
  snprintf(offsets_name, FILENAME_MAX, "%s_stencil_offsets", stencil_name);
  size_t size;
  s->offsets = silo_file_read_int_array(file, offsets_name, &size);
  s->num_indices = (int)(size) - 1;

  if (s->offsets[s->num_indices] > 0)
  {
    char indices_name[FILENAME_MAX+1];
    snprintf(indices_name, FILENAME_MAX, "%s_stencil_indices", stencil_name);
    s->indices = silo_file_read_int_array(file, indices_name, &size);
    ASSERT((int)size == s->offsets[s->num_indices]);
  }
  else
    s->indices = NULL;

  char ex_name[FILENAME_MAX+1];
  snprintf(ex_name, FILENAME_MAX, "%s_stencil_ex", stencil_name);
  s->ex = silo_file_read_exchanger(file, ex_name, comm);
  STOP_FUNCTION_TIMER();
  return s;
}

bool silo_file_contains_neighbor_pairing(silo_file_t* file, 
                                         const char* neighbors_name)
{
  char name[FILENAME_MAX+1];
  snprintf(name, FILENAME_MAX, "%s_neighbor_pairing_name", neighbors_name);
  return silo_file_contains_string(file, name);
}

void silo_file_write_neighbor_pairing(silo_file_t* file,
                                      const char* neighbors_name,
                                      neighbor_pairing_t* neighbors)
{
  char name_name[FILENAME_MAX];
  snprintf(name_name, FILENAME_MAX, "%s_neighbor_pairing_name", neighbors_name);
  silo_file_write_string(file, name_name, neighbors->name);
  char pairs_name[FILENAME_MAX];
  snprintf(pairs_name, FILENAME_MAX, "%s_neighbor_pairing_pairs", neighbors_name);
  silo_file_write_int_array(file, pairs_name, neighbors->pairs, 2*neighbors->num_pairs);

  if (neighbors->ex != NULL)
  {
    char ex_name[FILENAME_MAX];
    snprintf(ex_name, FILENAME_MAX, "%s_neighbor_pairing_ex", neighbors_name);
    silo_file_write_exchanger(file, ex_name, neighbors->ex);
  }
}

neighbor_pairing_t* silo_file_read_neighbor_pairing(silo_file_t* file,
                                                    const char* neighbors_name,
                                                    MPI_Comm comm)
{
  neighbor_pairing_t* p = polymec_malloc(sizeof(neighbor_pairing_t));
  char name_name[FILENAME_MAX];
  snprintf(name_name, FILENAME_MAX, "%s_neighbor_pairing_name", neighbors_name);
  p->name = silo_file_read_string(file, name_name);
  char pairs_name[FILENAME_MAX];
  snprintf(pairs_name, FILENAME_MAX, "%s_neighbor_pairing_pairs", neighbors_name);
  size_t size;
  p->pairs = silo_file_read_int_array(file, pairs_name, &size);
  ASSERT((size % 2) == 0);
  p->num_pairs = (int)size/2;

  char ex_name[FILENAME_MAX];
  snprintf(ex_name, FILENAME_MAX, "%s_neighbor_pairing_ex", neighbors_name);
  p->ex = silo_file_read_exchanger(file, ex_name, comm);
  return p;
}

void silo_file_write_field_metadata(silo_file_t* file,
                                    const char* md_name,
                                    field_metadata_t* md)
{
  silo_file_push_domain_dir(file);

  // Pack any metadata into an array.
  size_t md_size = 1;
  int num_comps = field_metadata_num_components(md);
  for (int c = 0; c < num_comps; ++c)
  {
    const char* name = field_metadata_name(md, c);
    size_t name_size = (name != NULL) ? strlen(name) : 0;
    const char* units = field_metadata_units(md, c);
    size_t units_size = (units != NULL) ? strlen(units) : 0;
    md_size += 4 + name_size + units_size;
  }

  int pos = 0, comp, num_vectors = 0;
  ++md_size;
  while (field_metadata_next_vector(md, &pos, &comp)) 
  {
    ++num_vectors;
    ++md_size;
  }

  int num_tensor2s = 0;
  pos = 0;
  ++md_size;
  while (field_metadata_next_tensor2(md, &pos, &comp)) 
  {
    ++num_tensor2s;
    ++md_size;
  }

  int num_symtensor2s = 0;
  pos = 0;
  ++md_size;
  while (field_metadata_next_symtensor2(md, &pos, &comp)) 
  {
    ++num_symtensor2s;
    ++md_size;
  }

  int mda[md_size];
  memset(mda, 0, sizeof(int) * md_size);
  mda[0] = num_comps;
  int offset = 1;
  for (int c = 0; c < num_comps; ++c)
  {
    const char* name = field_metadata_name(md, c);
    int name_size = (int)((name != NULL) ? strlen(name) : 0);
    mda[offset++] = name_size; 
    for (int cc = 0; cc < name_size; ++cc, ++offset)
      mda[offset] = (int)(name[cc]);

    const char* units = field_metadata_units(md, c);
    int units_size = (int)((units != NULL) ? strlen(units) : 0);
    mda[offset++] = units_size; 
    for (int cc = 0; cc < units_size; ++cc, ++offset)
      mda[offset] = (int)(units[cc]);

    mda[offset++] = (int)(field_metadata_conserved(md, c));
    mda[offset++] = (int)(field_metadata_extensive(md, c));
  }
  mda[offset++] = num_vectors;
  pos = 0;
  while (field_metadata_next_vector(md, &pos, &comp))
    mda[offset++] = comp;
  mda[offset++] = num_tensor2s;
  pos = 0;
  while (field_metadata_next_tensor2(md, &pos, &comp))
    mda[offset++] = comp;
  mda[offset++] = num_symtensor2s;
  pos = 0;
  while (field_metadata_next_symtensor2(md, &pos, &comp))
    mda[offset++] = comp;
  ASSERT(offset == md_size);
  silo_file_write_int_array(file, md_name, mda, md_size);

  silo_file_pop_dir(file);
}

void silo_file_read_field_metadata(silo_file_t* file,
                                   const char* md_name,
                                   field_metadata_t* md)
{
  silo_file_push_domain_dir(file);

  size_t mda_size;
  int* mda = silo_file_read_int_array(file, md_name, &mda_size);
  if (mda != NULL)
  {
    int num_comps = mda[0];
    if (num_comps != field_metadata_num_components(md))
    {
      polymec_free(mda);
      polymec_error("silo_file_read_field_metadata: Inconsistent metadata size: %d components read, %d expected", 
                    num_comps, field_metadata_num_components(md));
    }
    int offset = 1;
    for (int c = 0; c < num_comps; ++c)
    {
      int name_size = mda[offset++];
      if (name_size > 0)
      {
        char name[name_size+1];
        for (int i = 0; i < name_size; ++i, ++offset)
          name[i] = (char)mda[offset];
        name[name_size] = '\0';
        field_metadata_set_name(md, c, name);
      }
      else
        field_metadata_set_name(md, c, NULL);

      int units_size = mda[offset++];
      if (name_size > 0)
      {
        char units[units_size+1];
        for (int i = 0; i < units_size; ++i, ++offset)
          units[i] = (char)mda[offset];
        units[units_size] = '\0';
        field_metadata_set_units(md, c, units);
      }
      else
        field_metadata_set_units(md, c, NULL);

      field_metadata_set_conserved(md, c, mda[offset++]);
      field_metadata_set_extensive(md, c, mda[offset++]);
    }
    int num_vectors = mda[offset++];
    for (int i = 0; i < num_vectors; ++i)
      field_metadata_set_vector(md, mda[offset++]);
    int num_tensor2s = mda[offset++];
    for (int i = 0; i < num_tensor2s; ++i)
      field_metadata_set_tensor2(md, mda[offset++]);
    int num_symtensor2s = mda[offset++];
    for (int i = 0; i < num_symtensor2s; ++i)
      field_metadata_set_symtensor2(md, mda[offset++]);
    polymec_free(mda);
  }

  silo_file_pop_dir(file);
}

bool silo_file_contains_field_metadata(silo_file_t* file,
                                       const char* md_name)
{
  silo_file_push_domain_dir(file);
  char arr_name[FILENAME_MAX+1];
  snprintf(arr_name, FILENAME_MAX, "%s_int_array", md_name);
  int result = DBInqVarExists(file->dbfile, arr_name);
  silo_file_pop_dir(file);
  return (result != 0);
}
