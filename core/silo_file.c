// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <sys/stat.h>
#include <dirent.h>
#include <gc/gc.h>
#include "core/polymec.h"
#if POLYMEC_HAVE_MPI
#include "mpi.h"
#include "pmpio.h"
#endif
#include "silo.h"
#include "core/arch.h"
#include "core/logging.h"
#include "core/silo_file.h"
#include "core/array.h"
#include "core/array_utils.h"
#include "core/timer.h"

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

// Access the underlying SILO file descriptor.
DBfile* silo_file_dbfile(silo_file_t* file);

// Creates a SILO options list from the given metadata object.
DBoptlist* optlist_from_metadata(silo_field_metadata_t* metadata);

// Clones a metadata-related SILO options list.
DBoptlist* optlist_clone(DBoptlist* optlist);

// Frees memory associated with a metadata-related SILO options list.
void optlist_free(DBoptlist* optlist);

// Writes metadata identifying a subdomain of a global mesh of the given 
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

//-------------------------------------------------------------------------
// End unpublished functions
//-------------------------------------------------------------------------

void silo_enable_compression(int level)
{
  ASSERT(level >= 0);
  ASSERT(level <= 9);
  // Global compression level for Silo.
  static int silo_compression_level = -1; // Disabled by default.
  if (level != silo_compression_level)
  {
    char options[1024];
    snprintf(options, 1024, "METHOD=GZIP,LEVEL=%d", level);
    DBSetCompression(options);
    silo_compression_level = level; 
  }
}

struct silo_field_metadata_t 
{
  char* label;
  char* units;
  bool conserved;
  bool extensive;
  int vector_component;
};

static void silo_field_metadata_free(void* ctx, void* dummy)
{
  silo_field_metadata_t* metadata = ctx;
  if (metadata->label)
    string_free(metadata->label);
  if (metadata->units)
    string_free(metadata->units);
}

silo_field_metadata_t* silo_field_metadata_new()
{
  silo_field_metadata_t* metadata = GC_MALLOC(sizeof(silo_field_metadata_t));
  metadata->label = NULL;
  metadata->units = NULL;
  metadata->conserved = false;
  metadata->extensive = false;
  metadata->vector_component = -1;
  GC_register_finalizer(metadata, silo_field_metadata_free, metadata, NULL, NULL);
  return metadata;
}

void silo_field_metadata_set_label(silo_field_metadata_t* metadata,
                                   const char* label)
{
  if (metadata->label != NULL)
    string_free(metadata->label);
  if (label != NULL)
    metadata->label = string_dup(label);
  else
    metadata->label = NULL;
}

void silo_field_metadata_set_units(silo_field_metadata_t* metadata,
                                   const char* units)
{
  if (metadata->units != NULL)
    string_free(metadata->units);
  if (units != NULL)
    metadata->units = string_dup(units);
  else
    metadata->units = NULL;
}

void silo_field_metadata_set_conserved(silo_field_metadata_t* metadata,
                                       bool is_conserved)
{
  metadata->conserved = is_conserved;
}

void silo_field_metadata_set_extensive(silo_field_metadata_t* metadata,
                                       bool is_extensive)
{
  metadata->extensive = is_extensive;
}

void silo_field_metadata_set_vector_component(silo_field_metadata_t* metadata,
                                              int component)
{
  ASSERT(component >= -1);
  ASSERT(component < 3);
  metadata->vector_component = component;
}

char* silo_field_metadata_label(silo_field_metadata_t* metadata)
{
  return metadata->label;
}

char* silo_field_metadata_units(silo_field_metadata_t* metadata)
{
  return metadata->units;
}

bool silo_field_metadata_conserved(silo_field_metadata_t* metadata)
{
  return metadata->conserved;
}

bool silo_field_metadata_extensive(silo_field_metadata_t* metadata)
{
  return metadata->extensive;
}

bool silo_field_metadata_is_vector_component(silo_field_metadata_t* metadata)
{
  return (metadata->vector_component != -1);
}

int silo_field_metadata_vector_component(silo_field_metadata_t* metadata)
{
  return metadata->vector_component;
}

// These helpers are used in the read/write operations.
static void read_mesh_metadata(DBucdvar* var, silo_field_metadata_t* metadata)
{
  if (metadata != NULL)
  {
    silo_field_metadata_set_label(metadata, var->label);
    silo_field_metadata_set_units(metadata, var->units);
    silo_field_metadata_set_conserved(metadata, var->conserved);
    silo_field_metadata_set_extensive(metadata, var->extensive);
  }
}

static void read_point_metadata(DBmeshvar* var, silo_field_metadata_t* metadata)
{
  if (metadata != NULL)
  {
    silo_field_metadata_set_label(metadata, var->label);
    silo_field_metadata_set_units(metadata, var->units);
    silo_field_metadata_set_conserved(metadata, var->conserved);
    silo_field_metadata_set_extensive(metadata, var->extensive);
  }
}

DBoptlist* optlist_from_metadata(silo_field_metadata_t* metadata)
{
  DBoptlist* optlist = NULL;
  if (metadata != NULL)
  {
    optlist = DBMakeOptlist(4);
    if (metadata->label != NULL)
      DBAddOption(optlist, DBOPT_LABEL, string_dup(metadata->label));
    if (metadata->units != NULL)
      DBAddOption(optlist, DBOPT_UNITS, string_dup(metadata->units));
    int* conserved = polymec_malloc(sizeof(int));
    *conserved = metadata->conserved;
    int* extensive = polymec_malloc(sizeof(int));
    *extensive = metadata->extensive;
    DBAddOption(optlist, DBOPT_CONSERVED, conserved);
    DBAddOption(optlist, DBOPT_EXTENSIVE, extensive);
  }
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
  DBSetDir(file, dir_name);
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
    DBSetDir(file, dir_name);
  }
  else
  {
    file = DBOpen(filename, driver, DB_READ);
    DBSetDir(file, dir_name);
  }
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
  polymec_free(mesh);
  optlist_free(mesh->optlist);
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
                     int_slist_t* cycles)
{
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
      return false;
    }
    string_slist_t* files_in_dir = files_within_directory(directory);

    // Try to find a master file or single data file.
    char data_file[FILENAME_MAX];
    data_file[0] = '\0';
    bool found_cycles = false;
    {
      string_slist_node_t* node = files_in_dir->front;
      while (node != NULL)
      {
        char path[FILENAME_MAX];
        snprintf(path, FILENAME_MAX, "%s", file_prefix); 
        if (strstr(node->value, path) && strstr(node->value, "silo"))
        {
          strncpy(data_file, node->value, FILENAME_MAX);
          snprintf(data_file, FILENAME_MAX, "%s/%s", directory, node->value);
          snprintf(path, FILENAME_MAX, "%s-", file_prefix); 
          if (strstr(node->value, path))
            found_cycles = true;
          break;
        }
        node = node->next;
      }
    }
    if (strlen(data_file) == 0)
    {
      log_urgent("silo_file_query: Could not find %s/%s-*.silo.", directory, file_prefix);
      string_slist_free(files_in_dir);
      return false;
    }

    // Open up the file and see whether it's a master file or a serial data file.
    bool is_master = false;
    int driver = DB_HDF5;
    DBfile* file = DBOpen(data_file, driver, DB_READ);

    // What's in there?
    DBtoc* toc = DBGetToc(file); 
    if ((toc->nucdmesh == 0) && (toc->nptmesh == 0) && (toc->nmultimesh > 0))
      is_master = true;
    if (is_master)
    {
      if (!DBInqVarExists(file, "POLYMEC_SILO_MASTER_FILE"))
      {
        log_urgent("silo_query_file: invalid Silo master file.");
        DBClose(file);
        return false;
      }

      // How many MPI processes were used to construct the data set?
      int my_num_mpi_procs = -1;
      for (int f = 0; f < toc->nmultimesh; ++f)
      {
        DBmultimesh* subdomain_mesh = DBGetMultimesh(file, toc->multimesh_names[f]);
        if (my_num_mpi_procs == -1)
          my_num_mpi_procs = subdomain_mesh->nblocks;
        ASSERT(my_num_mpi_procs == subdomain_mesh->nblocks);
      }
      *num_mpi_processes = my_num_mpi_procs;

      // How many files are in the data set?
      if (DBInqVarExists(file, "num_files"))
        DBReadVar(file, "num_files", num_files);
      else
      {
        log_urgent("silo_file_query: Could not read number of files in set.");
        DBClose(file);
        return false;
      }
    }
    else
    {
      if (!DBInqVarExists(file, "POLYMEC_SILO_FILE"))
      {
        log_urgent("silo_query_file: Invalid Silo file.");
        DBClose(file);
        return false;
      }

      // A single data file can only be written for a serial run.
      *num_files = 1;
      *num_mpi_processes = 1;
    }

    DBClose(file);

    // Search for available cycles.
    if ((cycles != NULL) && found_cycles)
    {
      int_slist_clear(cycles);
      string_slist_node_t* node = files_in_dir->front;
      while (node != NULL)
      {
        char path[FILENAME_MAX];
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
            int_slist_append(cycles, atoi(num));
        }
        node = node->next;
      }
    }

    // Clean up.
    string_slist_free(files_in_dir);
  }

  // Now spread the word to other processes.
  int num_cycles = (cycles != NULL) ? cycles->size : 0;
  int data[3] = {*num_files, *num_mpi_processes, num_cycles};
  MPI_Bcast(data, 3, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank != 0)
  {
    *num_files = data[0];
    *num_mpi_processes = data[1];
    num_cycles = data[2];
  }

  if (cycles != NULL)
  {
    int cycles_buffer[num_cycles];
    if (rank == 0)
    {
      // Spit the cycles into the array and sort them.
      int_slist_node_t* node = cycles->front;
      int i = 0;
      while (node != NULL)
      {
        cycles_buffer[i++] = node->value;
        node = node->next;
      }
      ASSERT(i == num_cycles);
      int_qsort(cycles_buffer, num_cycles);
    }
    MPI_Bcast(cycles_buffer, num_cycles, MPI_INT, 0, MPI_COMM_WORLD);

    // Shuffle them back into our linked list.
    int_slist_clear(cycles);
    for (int i = 0; i < num_cycles; ++i)
      int_slist_append(cycles, cycles_buffer[i]);
  }
  return true;
}

struct silo_file_t 
{
  // File data.
  DBfile* dbfile;

  // Metadata.
  char prefix[FILENAME_MAX], directory[FILENAME_MAX], filename[FILENAME_MAX];
  int cycle;
  real_t time;
  int mode; // Open for reading (DB_READ) or writing (DB_CLOBBER)? 
  string_ptr_unordered_map_t* expressions;

#if POLYMEC_HAVE_MPI
  // Stuff for poor man's parallel I/O.
  PMPIO_baton_t* baton;
  MPI_Comm comm;
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

  silo_file_write_string(file, "provenance", provenance_str);
  string_free(provenance_str);
}

#if POLYMEC_HAVE_MPI
static void write_subdomains_to_file(silo_file_t* file)
{
  ASSERT(file->mode == DB_CLOBBER);

  if (file->nproc == 1) return;
  if (file->rank_in_group != 0) return;
  int num_chunks = file->nproc / file->num_files;

  // Stick in cycle/time information if needed.
  DBoptlist* optlist = DBMakeOptlist(2);
  if (file->cycle >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &file->cycle);
  if (file->time != -FLT_MAX)
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
      char mesh_name[FILENAME_MAX];
      snprintf(mesh_name, FILENAME_MAX, "domain_%d/%s", j, mesh->name);
      mesh_names[j] = string_dup(mesh_name);
      mesh_types[j] = mesh->type;
    }

    // Write the point mesh and variable data.
    DBSetDir(file->dbfile, "/");
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
      char field_name[FILENAME_MAX];
      snprintf(field_name, FILENAME_MAX, "domain_%d/%s", j, field->name);
      field_names[j] = string_dup(field_name);
      field_types[j] = field->type;
    }

    // Write the field data.
    DBSetDir(file->dbfile, "/");
    DBPutMultivar(file->dbfile, field->name, num_chunks, 
                  (char const* const*)field_names, field_types, field->optlist); 

    // Clean up.
    for (int j = 0; j < num_chunks; ++j)
      polymec_free(field_names[j]);
  }
}

static void write_master_file(silo_file_t* file)
{
  ASSERT(file->mode == DB_CLOBBER);

  // FIXME: Should change this to use Silo's name schemes for multi-block 
  // FIXME: objects when we start to Get Real Parallel.

  char master_file_name[FILENAME_MAX];
  if (file->cycle == -1)
    snprintf(master_file_name, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
  else
    snprintf(master_file_name, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, file->cycle);
  PMPIO_baton_t* baton = PMPIO_Init(file->num_files, PMPIO_WRITE, file->comm, file->mpi_tag+1, 
                                    pmpio_create_file, pmpio_open_file, 
                                    pmpio_close_file, 0);
  DBfile* master = (DBfile*)PMPIO_WaitForBaton(baton, master_file_name, "/");

  // Write our stamp of approval.
  int one = 1;
  DBWrite(master, "POLYMEC_SILO_MASTER_FILE", &one, &one, 1, DB_INT);

  // Stick in cycle/time information if needed.
  DBoptlist* optlist = DBMakeOptlist(2);
  if (file->cycle >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &file->cycle);
  if (file->time != -FLT_MAX)
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
        char mesh_name[FILENAME_MAX];
        mesh_types[num_chunks*j+c] = mesh->type;
        if (file->cycle == -1)
          snprintf(mesh_name, FILENAME_MAX, "%d/%s.silo:/domain_%d/%s", j, file->prefix, c, mesh->name);
        else
          snprintf(mesh_name, FILENAME_MAX, "%d/%s-%d.silo:/domain_%d/%s", j, file->prefix, file->cycle, c, mesh->name);
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
        char field_name[FILENAME_MAX];
        if (file->cycle == -1)
          snprintf(field_name, FILENAME_MAX, "%d/%s.silo:/domain_%d/%s", j, file->prefix, c, field->name);
        else
          snprintf(field_name, FILENAME_MAX, "%d/%s-%d.silo:/domain_%d/%s", j, file->prefix, file->cycle, c, field->name);
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

  PMPIO_HandOffBaton(baton, (void*)master);
  PMPIO_Finish(baton);
}
#endif

silo_file_t* silo_file_new(MPI_Comm comm,
                           const char* file_prefix,
                           const char* directory,
                           int num_files,
                           int mpi_tag,
                           int cycle,
                           real_t time)
{
  START_FUNCTION_TIMER();
  silo_file_t* file = polymec_malloc(sizeof(silo_file_t));
  file->expressions = string_ptr_unordered_map_new();

  // Strip .silo off of the prefix if it's there.
  {
    char prefix[FILENAME_MAX];
    strncpy(prefix, file_prefix, FILENAME_MAX);
    char* suffix = strstr(prefix, ".silo");
    if (suffix != NULL)
      suffix[0] = '\0';
    strcpy(file->prefix, prefix);
  }

#if POLYMEC_HAVE_MPI
  file->comm = comm;
  MPI_Comm_size(file->comm, &file->nproc);
  MPI_Comm_rank(file->comm, &file->rank);
  if (num_files == -1)
    file->num_files = file->nproc;
  else
    file->num_files = num_files;
  ASSERT(file->num_files <= file->nproc);
  file->mpi_tag = mpi_tag;

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

    // Create a subdirectory for each group.
    char group_dir_name[FILENAME_MAX];
    snprintf(group_dir_name, FILENAME_MAX, "%s/%d", file->directory, file->group_rank);
    if (file->rank_in_group == 0)
      create_directory(group_dir_name, S_IRWXU | S_IRWXG);

    // Determine a file name and directory name.
    if (cycle == -1)
      snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", group_dir_name, file->prefix);
    else
      snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", group_dir_name, file->prefix, cycle);
    char silo_dir_name[FILENAME_MAX];
    snprintf(silo_dir_name, FILENAME_MAX, "domain_%d", file->rank_in_group);
    file->dbfile = (DBfile*)PMPIO_WaitForBaton(file->baton, file->filename, silo_dir_name);

    file->subdomain_meshes = ptr_array_new();
    file->subdomain_fields = ptr_array_new();
  }
  else
  {
    if (strlen(directory) == 0)
      strncpy(file->directory, ".", FILENAME_MAX);
    else
      strncpy(file->directory, directory, FILENAME_MAX);

    if (cycle == -1)
      snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
    else
      snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, cycle);

    int driver = DB_HDF5;
    if (strcmp(file->directory, ".") != 0)
      create_directory(file->directory, S_IRWXU | S_IRWXG);
    file->dbfile = DBCreate(file->filename, DB_CLOBBER, DB_LOCAL, NULL, driver);
    DBSetDir(file->dbfile, "/");
  }
#else
  if (strlen(directory) == 0)
    strncpy(file->directory, ".", FILENAME_MAX);
  else
    strncpy(file->directory, directory, FILENAME_MAX);

  if (cycle == -1)
    snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
  else
    snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, cycle);

  int driver = DB_HDF5;
  create_directory(file->directory, S_IRWXU | S_IRWXG);
  file->dbfile = DBCreate(file->filename, DB_CLOBBER, DB_LOCAL, NULL, driver);
  DBSetDir(file->dbfile, "/");
#endif

  file->mode = DB_CLOBBER;
  file->cycle = cycle;
  file->time = time;

  // Write our stamp of approval.
  int one = 1;
  DBWrite(file->dbfile, "POLYMEC_SILO_FILE", &one, &one, 1, DB_INT);

  STOP_FUNCTION_TIMER();
  return file;
}

static void show_provenance_on_debug_log(silo_file_t* file)
{
  // If we're printing debug messages, show provenance of the file.
  if (log_level() == LOG_DEBUG)
  {
    char* provenance = silo_file_read_string(file, "provenance");

    // We may need to temporarily bump up the debug log's buffness.
    int max_message_size, flush_every;
    get_log_buffering(LOG_DEBUG, &max_message_size, &flush_every);
    set_log_buffering(LOG_DEBUG, strlen(provenance)+10, 1); // SUPER DUPER!!!
    log_debug_literal(provenance); // No formatting.
    set_log_buffering(LOG_DEBUG, max_message_size, flush_every); // Back to normal

    string_free(provenance);
  }
}

silo_file_t* silo_file_open(MPI_Comm comm,
                            const char* file_prefix,
                            const char* directory,
                            int mpi_tag,
                            int cycle, 
                            real_t* time)
{
  START_FUNCTION_TIMER();
  silo_file_t* file = polymec_malloc(sizeof(silo_file_t));
  file->mode = DB_READ;
  file->cycle = -1;
  file->time = -FLT_MAX;
  file->expressions = NULL;

  // Strip .silo off of the prefix if it's there.
  {
    char prefix[FILENAME_MAX];
    strncpy(prefix, file_prefix, FILENAME_MAX);
    char* suffix = strstr(prefix, ".silo");
    if (suffix != NULL)
      suffix[0] = '\0';
    strcpy(file->prefix, prefix);
  }

  // Query the dataset for the number of files and MPI processes and cycles.
  int num_files, num_mpi_procs;
  int_slist_t* cycles = int_slist_new();
  if (!silo_file_query(file_prefix, directory, &num_files, &num_mpi_procs, cycles))
    polymec_error("silo_file_open: Invalid file.");

  log_debug("silo_file_open: Opened file written by %d MPI processes.", num_mpi_procs);

  // For now, we only support reading files that were written with the same 
  // number of processes.
  int nproc = 1;
#if POLYMEC_HAVE_MPI
  file->comm = comm;
  MPI_Comm_size(file->comm, &nproc); 
#endif
  if (nproc != num_mpi_procs)
    polymec_not_implemented("silo_file_open: reading files written with different\n"
                            "number of MPI processes is not yet supported.");

  // Check to see whether the requested cycle is available, or whether the 
  // latest one is requested (with -1).
  if (cycle >= 0)
  {
    bool cycle_found = false;
    int_slist_node_t* node = cycles->front;
    while (node != NULL)
    {
      if (node->value == cycle)
      {
        cycle_found = true;
        break;
      }
      else if (node->value > cycle) // cycles are sorted
        break;
      node = node->next;
    }
    if (!cycle_found)
      polymec_error("silo_file_open: Cycle %d was not found for prefix '%s' in directory %s.", cycle, file->prefix, directory);
  }

#if POLYMEC_HAVE_MPI
  // The way these things are defined for a file has to do with how the 
  // file was generated, not how we are currently running.
  file->comm = comm; // ...for lack of a better value. Plus, might be useful.
  MPI_Comm_rank(file->comm, &file->rank); // ...also might be useful.
  file->num_files = num_files; // number of files in the data set.
  file->nproc = num_mpi_procs; // number of MPI procs used to write the thing.
  file->mpi_tag = mpi_tag; // this is fine.

  if (file->nproc > 1)
  {
    // Look for the master directory.
    if (strlen(directory) == 0)
      snprintf(file->directory, FILENAME_MAX, "%s_%dprocs", file->prefix, file->nproc);
    else
      strncpy(file->directory, directory, FILENAME_MAX);
    if (file->rank == 0)
    {
      if (!directory_exists(file->directory))
      {
        polymec_error("silo_file_open: Master directory %s does not exist for file prefix %s.",
            file->directory, file->prefix);
      }
    }

    // Initialize poor man's I/O and figure out group ranks.
    MPI_Barrier(file->comm);
    file->baton = PMPIO_Init(file->num_files, PMPIO_READ, file->comm, file->mpi_tag, 
                             pmpio_create_file, pmpio_open_file, pmpio_close_file, 0);
    file->group_rank = PMPIO_GroupRank(file->baton, file->rank);
    file->rank_in_group = PMPIO_RankInGroup(file->baton, file->rank);

    // Make sure a subdirectory exists for each group.
    char group_dir_name[FILENAME_MAX];
    snprintf(group_dir_name, FILENAME_MAX, "%s/%d", file->directory, file->group_rank);
    if (file->rank_in_group == 0)
    {
      DIR* group_dir = opendir(group_dir_name);
      if (group_dir == NULL)
      {
        polymec_error("silo_file_open: Group directory %s does not exist for file prefix %s.",
            group_dir_name, file->prefix);
      }
      else
        closedir(group_dir);
    }

    // Determine a file name and directory name.
    if (cycle == -1)
      snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", group_dir_name, file->prefix);
    else
      snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", group_dir_name, file->prefix, cycle);
    char silo_dir_name[FILENAME_MAX];
    snprintf(silo_dir_name, FILENAME_MAX, "domain_%d", file->rank_in_group);
    file->dbfile = (DBfile*)PMPIO_WaitForBaton(file->baton, file->filename, silo_dir_name);

    show_provenance_on_debug_log(file);

    file->subdomain_meshes = ptr_array_new();
    file->subdomain_fields = ptr_array_new();
  }
  else
  {
    if (strlen(directory) == 0)
      strncpy(file->directory, ".", FILENAME_MAX);
    else
      strncpy(file->directory, directory, FILENAME_MAX);

    if (cycle == -1)
      snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
    else
      snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, cycle);

    int driver = DB_HDF5;
    file->dbfile = DBOpen(file->filename, driver, file->mode);
    DBSetDir(file->dbfile, "/");

    show_provenance_on_debug_log(file);
  }
#else
  if (strlen(directory) == 0)
    strncpy(file->directory, ".", FILENAME_MAX);
  else
    strncpy(file->directory, directory, FILENAME_MAX);

  if (cycle == -1)
    snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
  else
    snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, cycle);

  int driver = DB_HDF5;
  file->dbfile = DBOpen(file->filename, driver, file->mode);
  DBSetDir(file->dbfile, "/");

  show_provenance_on_debug_log(file);
#endif

  // Get cycle/time information.
  if (DBInqVarExists(file->dbfile, "dtime"))
  {
    double dtime;
    DBReadVar(file->dbfile, "dtime", &dtime);
    file->time = (real_t)dtime;
  }
  else
    file->time = 0.0;
  if (DBInqVarExists(file->dbfile, "cycle"))
    DBReadVar(file->dbfile, "cycle", &file->cycle);
  else
    file->cycle = -1;

  if (time != NULL)
    *time = file->time;

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

    PMPIO_HandOffBaton(file->baton, (void*)file->dbfile);
    PMPIO_Finish(file->baton);

    if (file->mode == DB_CLOBBER)
    {
      // Write the uber-master file containing any multiobjects if need be.
      write_master_file(file);
    }
    MPI_Barrier(file->comm);

    ptr_array_free(file->subdomain_meshes);
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

  // Clean up.
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

  int pos = 0, *tag, tag_size;
  char* tag_name;
  while (mesh_next_tag(tagger, &pos, &tag_name, &tag, &tag_size))
  {
    int_array_append(elem_lengths, tag_size);
    string_array_append(elem_names, tag_name);
    for (int i = 0; i < tag_size; ++i)
      int_array_append(tag_data, tag[i]);
  }

  // Write the compound array.
  if (elem_names->size > 0)
  {
    DBPutCompoundarray(file->dbfile, tag_list_name, 
                       (char const* const*)elem_names->data, elem_lengths->data,
                       elem_names->size, tag_data->data, tag_data->size, DB_INT, 0);
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
        int size;
        tag = tagger_tag(tagger, (const char*)tag_names[i], &size);
        ASSERT(size == tag_sizes[i]);
      }
      memcpy(tag, &array[j], sizeof(int) * tag_sizes[i]);
      j += tag_sizes[i];
    }
    ASSERT(j == var->nvalues);
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
#if POLYMEC_HAVE_MPI
  int_array_append(array, file->nproc);
  int_array_append(array, file->rank);
#endif

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
  int size;
  int* array = silo_file_read_int_array(file, exchanger_name, &size);

  // Create the exchanger.
  exchanger_t* ex = exchanger_new(comm);
  int i = 0;
#if POLYMEC_HAVE_MPI
  ASSERT(array[i++] == file->nproc);
  ASSERT(array[i++] == file->rank);
#endif
  int num_sends = array[i++];
  for (int j = 0; j < num_sends; ++j)
  {
    int proc = array[i++];
    int num_indices = array[i++];
    exchanger_set_send(ex, proc, &array[i++], num_indices, true);
  }
  int num_receives = array[i++];
  for (int j = 0; j < num_receives; ++j)
  {
    int proc = array[i++];
    int num_indices = array[i++];
    exchanger_set_receive(ex, proc, &array[i++], num_indices, true);
  }
  ASSERT(i == size);

  // Clean up.
  polymec_free(array);
  STOP_FUNCTION_TIMER();

  return ex;
}

void silo_file_write_mesh(silo_file_t* file,
                          const char* mesh_name,
                          mesh_t* mesh)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_CLOBBER);

  // This is optional for now, but we'll give it anyway.
  char *coordnames[3];
  coordnames[0] = (char*)"xcoords";
  coordnames[1] = (char*)"ycoords";
  coordnames[2] = (char*)"zcoords";

  // Node coordinates.
  int num_nodes = mesh->num_nodes;
  double* x = polymec_malloc(sizeof(double) * num_nodes);
  double* y = polymec_malloc(sizeof(double) * num_nodes);
  double* z = polymec_malloc(sizeof(double) * num_nodes);
  for (int i = 0; i < num_nodes; ++i)
  {
    x[i] = (double)mesh->nodes[i].x;
    y[i] = (double)mesh->nodes[i].y;
    z[i] = (double)mesh->nodes[i].z;
  }
  double* coords[3];
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;

  // The polyhedral zone list is referred to in the options list.
  DBoptlist* optlist = DBMakeOptlist(10);
  char zonelist_name[FILENAME_MAX];
  snprintf(zonelist_name, FILENAME_MAX, "%s_zonelist", mesh_name);
  DBAddOption(optlist, DBOPT_PHZONELIST, zonelist_name);

  // Stick in cycle/time information if needed.
  if (file->cycle >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &file->cycle);
  if (file->time != -FLT_MAX)
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
    polymec_error("silo_file_write_mesh: Could not write mesh '%s'.", mesh_name);

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
  int* cell_face_counts = polymec_malloc(sizeof(int) * (num_cells + mesh->num_ghost_cells));
  memset(cell_face_counts, 0, sizeof(int) * (num_cells + mesh->num_ghost_cells));
  for (int i = 0; i < num_cells; ++i)
    cell_face_counts[i] = mesh->cell_face_offsets[i+1] - mesh->cell_face_offsets[i];

  // Write the connectivity information.
  result = DBPutPHZonelist(file->dbfile, zonelist_name, num_faces, face_node_counts,
                           mesh->face_node_offsets[num_faces], mesh->face_nodes,
                           ext_faces, num_cells + mesh->num_ghost_cells, cell_face_counts,
                           mesh->cell_face_offsets[num_cells], mesh->cell_faces,
                           0, 0, num_cells-1, optlist);
  if (result == -1)
    polymec_error("silo_file_write_mesh: Could not write connectivity data for mesh '%s'.", mesh_name);

  // Partial cleanup.
  polymec_free(face_node_counts);
  polymec_free(ext_faces);
  polymec_free(cell_face_counts);

  // Finally, write out the face_cells array.
  {
    char name[FILENAME_MAX];
    snprintf(name, FILENAME_MAX, "%s_face_cells", mesh_name);
    silo_file_write_int_array(file, name, mesh->face_cells, 2*mesh->num_faces);
  }

  // Write out tag information.
  {
    char tag_name[FILENAME_MAX];
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
    char ex_name[FILENAME_MAX];
    snprintf(ex_name, FILENAME_MAX, "%s_exchanger", mesh_name);
    silo_file_write_exchanger(file, ex_name, mesh_exchanger(mesh));
  }

  // Write out the number of mesh cells/faces/nodes to special variables.
  char num_cells_var[FILENAME_MAX], num_faces_var[FILENAME_MAX],
       num_nodes_var[FILENAME_MAX];
  snprintf(num_cells_var, FILENAME_MAX, "%s_mesh_num_cells", mesh_name);
  snprintf(num_faces_var, FILENAME_MAX, "%s_mesh_num_faces", mesh_name);
  snprintf(num_nodes_var, FILENAME_MAX, "%s_mesh_num_nodes", mesh_name);
  int one = 1;
  DBWrite(file->dbfile, num_cells_var, &mesh->num_cells, &one, 1, DB_INT);
  DBWrite(file->dbfile, num_faces_var, &mesh->num_faces, &one, 1, DB_INT);
  DBWrite(file->dbfile, num_nodes_var, &mesh->num_nodes, &one, 1, DB_INT);
  
  // Clean up.
  DBFreeOptlist(optlist);

#if POLYMEC_HAVE_MPI
  // For parallel environments, add a subdomain entry.
  if (file->nproc > 1)
    silo_file_add_subdomain_mesh(file, mesh_name, DB_UCDMESH, NULL);
#endif

  STOP_FUNCTION_TIMER();
}

mesh_t* silo_file_read_mesh(silo_file_t* file,
                            const char* mesh_name)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_READ);

  DBucdmesh* ucd_mesh = DBGetUcdmesh(file->dbfile, mesh_name);
  if (ucd_mesh == NULL)
    polymec_error("No mesh named '%s' was found within the Silo file.", mesh_name);
  ASSERT(ucd_mesh->ndims == 3);

  // Also get the polyhedral zone list.
  char phzl_name[FILENAME_MAX];
  snprintf(phzl_name, FILENAME_MAX, "%s_zonelist", mesh_name);
  DBphzonelist* ph_zonelist = DBGetPHZonelist(file->dbfile, phzl_name);
  if (ph_zonelist == NULL)
    polymec_error("Mesh '%s' is not a polymec polyhedral mesh.", mesh_name);

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
  mesh_t* mesh = mesh_new(comm, num_cells, num_ghost_cells, 
                          num_faces, num_nodes);

  // Set node positions.
  double* x = ucd_mesh->coords[0];
  double* y = ucd_mesh->coords[1];
  double* z = ucd_mesh->coords[2];
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
  mesh_reserve_connectivity_storage(mesh);

  // Read in the face_cells array.
  {
    char name[FILENAME_MAX];
    snprintf(name, FILENAME_MAX, "%s_face_cells", mesh_name);
    int num_face_cells;
    int* face_cells = silo_file_read_int_array(file, name, &num_face_cells);
    ASSERT(num_face_cells == 2*mesh->num_faces);
    memcpy(mesh->face_cells, face_cells, sizeof(int) * 2 * mesh->num_faces);
    polymec_free(face_cells);
  }

  // Fill in cell faces and face nodes.
  memcpy(mesh->cell_faces, ph_zonelist->facelist, sizeof(int) * mesh->cell_face_offsets[mesh->num_cells]);
  memcpy(mesh->face_nodes, ph_zonelist->nodelist, sizeof(int) * mesh->face_node_offsets[mesh->num_faces]);

  // Finish constructing the mesh.
  mesh_construct_edges(mesh);
  mesh_compute_geometry(mesh);

  // Read in tag information.
  {
    char tag_name[FILENAME_MAX];
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
    char ex_name[FILENAME_MAX];
    snprintf(ex_name, FILENAME_MAX, "%s_exchanger", mesh_name);
    mesh_set_exchanger(mesh, silo_file_read_exchanger(file, ex_name, mesh->comm));
  }

  // Clean up.
  DBFreeUcdmesh(ucd_mesh);
  DBFreePHZonelist(ph_zonelist);

  STOP_FUNCTION_TIMER();
  return mesh;
}

bool silo_file_contains_mesh(silo_file_t* file, const char* mesh_name)
{
  return (DBInqVarExists(file->dbfile, mesh_name) && 
          (DBInqVarType(file->dbfile, mesh_name) == DB_UCDMESH));
}

static void silo_file_write_mesh_field(silo_file_t* file,
                                       mesh_centering_t centering,
                                       const char** field_component_names,
                                       const char* mesh_name,
                                       real_t* field_data,
                                       int num_components,
                                       silo_field_metadata_t** field_metadata)
{
  ASSERT(file->mode == DB_CLOBBER);

  // How many elements does our mesh have?
  char num_elems_var[FILENAME_MAX];
  int cent;
  switch(centering)
  {
    case MESH_CELL: cent = DB_ZONECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_cells", mesh_name); break;
    case MESH_FACE: cent = DB_FACECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_faces", mesh_name); break;
    case MESH_EDGE: cent = DB_EDGECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_edges", mesh_name); break;
    case MESH_NODE: cent = DB_NODECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_nodes", mesh_name);
  }
  ASSERT(DBInqVarExists(file->dbfile, num_elems_var));
  int num_elems;
  DBReadVar(file->dbfile, num_elems_var, &num_elems);
  if (num_components == 1)
  {
    silo_field_metadata_t* metadata = (field_metadata != NULL) ? field_metadata[0] : NULL;
    DBoptlist* optlist = optlist_from_metadata(metadata);
    DBPutUcdvar1(file->dbfile, field_component_names[0], mesh_name, field_data, num_elems, NULL, 0, SILO_FLOAT_TYPE, cent, optlist);

#if POLYMEC_HAVE_MPI
  // Add a subdomain entry for parallel environments.
  if (file->nproc > 1)
    silo_file_add_subdomain_field(file, mesh_name, field_component_names[0], DB_UCDVAR, optlist);
#endif
  }
  else
  {
    real_t* comp_data = polymec_malloc(sizeof(real_t) * num_elems); 
    for (int c = 0; c < num_components; ++c)
    {
      for (int i = 0; i < num_elems; ++i)
        comp_data[i] = field_data[num_components*i+c];
      silo_field_metadata_t* metadata = (field_metadata != NULL) ? field_metadata[c] : NULL;
      DBoptlist* optlist = optlist_from_metadata(metadata);
      DBPutUcdvar1(file->dbfile, field_component_names[c], mesh_name, comp_data, num_elems, NULL, 0, SILO_FLOAT_TYPE, cent, optlist);

#if POLYMEC_HAVE_MPI
  // Add a subdomain entry for parallel environments.
  if (file->nproc > 1)
    silo_file_add_subdomain_field(file, mesh_name, field_component_names[c], DB_UCDVAR, optlist);
#endif
    }
    polymec_free(comp_data);
  }
}

static void silo_file_write_scalar_mesh_field(silo_file_t* file,
                                              mesh_centering_t centering,
                                              const char* field_name,
                                              const char* mesh_name,
                                              real_t* field_data,
                                              silo_field_metadata_t* field_metadata)
{
  silo_file_write_mesh_field(file, centering, &field_name, mesh_name, field_data, 1, (field_metadata != NULL) ? &field_metadata : NULL);
}

static real_t* silo_file_read_mesh_field(silo_file_t* file,
                                         mesh_centering_t centering,
                                         const char** field_component_names,
                                         const char* mesh_name,
                                         int num_components,
                                         silo_field_metadata_t** field_metadata)
{
  ASSERT(file->mode == DB_READ);

  // How many elements does our mesh have?
  char num_elems_var[FILENAME_MAX];
  int cent = DB_ZONECENT;
  switch(centering)
  {
    case MESH_CELL: cent = DB_ZONECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_cells", mesh_name); break;
    case MESH_FACE: cent = DB_FACECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_faces", mesh_name); break;
    case MESH_EDGE: cent = DB_EDGECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_edges", mesh_name); break;
    case MESH_NODE: cent = DB_NODECENT; snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_nodes", mesh_name);
  }
  ASSERT(DBInqVarExists(file->dbfile, num_elems_var));
  int num_elems;
  DBReadVar(file->dbfile, num_elems_var, &num_elems);
  real_t* field = polymec_malloc(sizeof(real_t) * num_components * num_elems); 
  if (num_components == 1)
  {
    DBucdvar* var = DBGetUcdvar(file->dbfile, (char*)field_component_names[0]);
    if (var == NULL)
      polymec_error("Field '%s' was not found in the Silo file.", field_component_names[0]);
    if (var->centering != cent)
      polymec_error("Field '%s' has the incorrect centering.", field_component_names[0]);
    memcpy(field, var->vals[0], sizeof(real_t) * num_elems);
    silo_field_metadata_t* metadata = (field_metadata[0] != NULL) ? field_metadata[0] : NULL;
    read_mesh_metadata(var, metadata);
    DBFreeUcdvar(var);
  }
  else
  {
    for (int c = 0; c < num_components; ++c)
    {
      DBucdvar* var = DBGetUcdvar(file->dbfile, (char*)field_component_names[c]);
      if (var == NULL)
        polymec_error("Field '%s' was not found in the Silo file.", field_component_names[c]);
      if (var->centering != cent)
        polymec_error("Field '%s' has the incorrect centering.", field_component_names[c]);
      real_t* data = var->vals[c];
      for (int i = 0; i < num_elems; ++i)
        field[num_components*i+c] = data[i];
      read_mesh_metadata(var, (field_metadata != NULL) ? field_metadata[c] : NULL);
      DBFreeUcdvar(var);
    }
  }
  return field;
}

static real_t* silo_file_read_scalar_mesh_field(silo_file_t* file,
                                                mesh_centering_t centering,
                                                const char* field_name,
                                                const char* mesh_name,
                                                silo_field_metadata_t* field_metadata)
{
  return silo_file_read_mesh_field(file, centering, &field_name, mesh_name, 1, (field_metadata != NULL) ? &field_metadata : NULL);
}

static bool silo_file_contains_mesh_field(silo_file_t* file, 
                                          mesh_centering_t centering, 
                                          const char* field_name, 
                                          const char* mesh_name)
{
  bool result = (silo_file_contains_mesh(file, mesh_name) &&  // mesh exists...
                 (DBInqVarType(file->dbfile, mesh_name) == DB_UCDMESH) &&  // mesh is a point cloud...
                 (DBInqVarExists(file->dbfile, field_name) &&  // field exists...
                  (DBInqVarType(file->dbfile, field_name) == DB_UCDVAR))); // field is actually a variable.

  if (result)
  {
    // Make sure that our field is associated with our mesh.
    char field_mesh_name[FILENAME_MAX];
    int stat = DBInqMeshname(file->dbfile, field_name, field_mesh_name);
    result = ((stat == 0) && (strcmp(mesh_name, field_mesh_name) == 0));
    if (result)
    {
      // Make sure the field has the right size.
      char num_elems_var[FILENAME_MAX];
      switch(centering)
      {
        case MESH_CELL: snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_cells", mesh_name); break;
        case MESH_FACE: snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_faces", mesh_name); break;
        case MESH_EDGE: snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_edges", mesh_name); break;
        case MESH_NODE: snprintf(num_elems_var, FILENAME_MAX, "%s_mesh_num_nodes", mesh_name);
      }
      int num_elems;
      int stat = DBReadVar(file->dbfile, num_elems_var, &num_elems);
      result = ((stat == 0) && (num_elems == DBGetVarLength(file->dbfile, field_name)));
    }
  }

  return result;
}

void silo_file_write_scalar_cell_field(silo_file_t* file,
                                       const char* field_name,
                                       const char* mesh_name,
                                       real_t* field_data,
                                       silo_field_metadata_t* field_metadata)
{
  START_FUNCTION_TIMER();
  silo_file_write_scalar_mesh_field(file, MESH_CELL, field_name, mesh_name, field_data, field_metadata);
  STOP_FUNCTION_TIMER();
}

real_t* silo_file_read_scalar_cell_field(silo_file_t* file,
                                         const char* field_name,
                                         const char* mesh_name,
                                         silo_field_metadata_t* field_metadata)
{
  START_FUNCTION_TIMER();
  real_t* field = silo_file_read_scalar_mesh_field(file, MESH_CELL, field_name, mesh_name, field_metadata);
  STOP_FUNCTION_TIMER();
  return field;
}


void silo_file_write_cell_field(silo_file_t* file,
                                const char** field_component_names,
                                const char* mesh_name,
                                real_t* field_data,
                                int num_components,
                                silo_field_metadata_t** field_metadata)
{
  START_FUNCTION_TIMER();
  silo_file_write_mesh_field(file, MESH_CELL, field_component_names, mesh_name, field_data,
                             num_components, field_metadata);
  STOP_FUNCTION_TIMER();
}

real_t* silo_file_read_cell_field(silo_file_t* file,
                                  const char** field_component_names,
                                  const char* mesh_name,
                                  int num_components,
                                  silo_field_metadata_t** field_metadata)
{
  START_FUNCTION_TIMER();
  real_t* field = silo_file_read_mesh_field(file, MESH_CELL, field_component_names, mesh_name, 
                                            num_components, field_metadata);
  STOP_FUNCTION_TIMER();
  return field;
}

bool silo_file_contains_cell_field(silo_file_t* file, const char* field_name, const char* mesh_name)
{
  START_FUNCTION_TIMER();
  bool result = silo_file_contains_mesh_field(file, MESH_CELL, field_name, mesh_name);
  STOP_FUNCTION_TIMER();
  return result;
}

void silo_file_write_scalar_face_field(silo_file_t* file,
                                       const char* field_name,
                                       const char* mesh_name,
                                       real_t* field_data,
                                       silo_field_metadata_t* field_metadata)
{
  START_FUNCTION_TIMER();
  silo_file_write_scalar_mesh_field(file, MESH_FACE, field_name, mesh_name, field_data, field_metadata);
  STOP_FUNCTION_TIMER();
}

real_t* silo_file_read_scalar_face_field(silo_file_t* file,
                                         const char* field_name,
                                         const char* mesh_name,
                                         silo_field_metadata_t* field_metadata)
{
  START_FUNCTION_TIMER();
  real_t* field = silo_file_read_scalar_mesh_field(file, MESH_FACE, field_name, mesh_name, field_metadata);
  STOP_FUNCTION_TIMER();
  return field;
}

void silo_file_write_face_field(silo_file_t* file,
                                const char** field_component_names,
                                const char* mesh_name,
                                real_t* field_data,
                                int num_components,
                                silo_field_metadata_t** field_metadata)
{
  START_FUNCTION_TIMER();
  silo_file_write_mesh_field(file, MESH_FACE, field_component_names, mesh_name, field_data,
                             num_components, field_metadata);
  STOP_FUNCTION_TIMER();
}

real_t* silo_file_read_face_field(silo_file_t* file,
                                  const char** field_component_names,
                                  const char* mesh_name,
                                  int num_components,
                                  silo_field_metadata_t** field_metadata)
{
  START_FUNCTION_TIMER();
  real_t* field = silo_file_read_mesh_field(file, MESH_FACE, field_component_names, mesh_name, 
                                            num_components, field_metadata);
  ASSERT(file->mode == DB_READ);
  STOP_FUNCTION_TIMER();
  return field;
}

bool silo_file_contains_face_field(silo_file_t* file, const char* field_name, const char* mesh_name)
{
  START_FUNCTION_TIMER();
  bool result = silo_file_contains_mesh_field(file, MESH_FACE, field_name, mesh_name);
  STOP_FUNCTION_TIMER();
  return result;
}

void silo_file_write_scalar_node_field(silo_file_t* file,
                                       const char* field_name,
                                       const char* mesh_name,
                                       real_t* field_data,
                                       silo_field_metadata_t* field_metadata)
{
  START_FUNCTION_TIMER();
  silo_file_write_scalar_mesh_field(file, MESH_NODE, field_name, mesh_name, field_data, field_metadata);
  STOP_FUNCTION_TIMER();
}

real_t* silo_file_read_scalar_node_field(silo_file_t* file,
                                         const char* field_name,
                                         const char* mesh_name,
                                         silo_field_metadata_t* field_metadata)
{
  START_FUNCTION_TIMER();
  real_t* field = silo_file_read_scalar_mesh_field(file, MESH_NODE, field_name, mesh_name, field_metadata);
  STOP_FUNCTION_TIMER();
  return field;
}

void silo_file_write_node_field(silo_file_t* file,
                                const char** field_component_names,
                                const char* mesh_name,
                                real_t* field_data,
                                int num_components,
                                silo_field_metadata_t** field_metadata)
{
  START_FUNCTION_TIMER();
  silo_file_write_mesh_field(file, MESH_NODE, field_component_names, mesh_name, field_data,
                             num_components, field_metadata);
  STOP_FUNCTION_TIMER();
}

real_t* silo_file_read_node_field(silo_file_t* file,
                                  const char** field_component_names,
                                  const char* mesh_name,
                                  int num_components,
                                  silo_field_metadata_t** field_metadata)
{
  START_FUNCTION_TIMER();
  real_t* field = silo_file_read_mesh_field(file, MESH_NODE, field_component_names, mesh_name, 
                                            num_components, field_metadata);
  STOP_FUNCTION_TIMER();
  return field;
}

bool silo_file_contains_node_field(silo_file_t* file, const char* field_name, const char* mesh_name)
{
  START_FUNCTION_TIMER();
  bool result = silo_file_contains_mesh_field(file, MESH_NODE, field_name, mesh_name);
  STOP_FUNCTION_TIMER();
  return result;
}

void silo_file_write_scalar_edge_field(silo_file_t* file,
                                       const char* field_name,
                                       const char* mesh_name,
                                       real_t* field_data,
                                       silo_field_metadata_t* field_metadata)
{
  START_FUNCTION_TIMER();
  silo_file_write_scalar_mesh_field(file, MESH_EDGE, field_name, mesh_name, field_data, field_metadata);
  STOP_FUNCTION_TIMER();
}

real_t* silo_file_read_scalar_edge_field(silo_file_t* file,
                                         const char* field_name,
                                         const char* mesh_name,
                                         silo_field_metadata_t* field_metadata)
{
  START_FUNCTION_TIMER();
  real_t* field = silo_file_read_scalar_mesh_field(file, MESH_EDGE, field_name, mesh_name, field_metadata);
  STOP_FUNCTION_TIMER();
  return field;
}

void silo_file_write_edge_field(silo_file_t* file,
                                const char** field_component_names,
                                const char* mesh_name,
                                real_t* field_data,
                                int num_components,
                                silo_field_metadata_t** field_metadata)
{
  START_FUNCTION_TIMER();
  silo_file_write_mesh_field(file, MESH_EDGE, field_component_names, mesh_name, field_data,
                             num_components, field_metadata);
  STOP_FUNCTION_TIMER();
}

real_t* silo_file_read_edge_field(silo_file_t* file,
                                  const char** field_component_names,
                                  const char* mesh_name,
                                  int num_components,
                                  silo_field_metadata_t** field_metadata)
{
  START_FUNCTION_TIMER();
  real_t* field = silo_file_read_mesh_field(file, MESH_EDGE, field_component_names, mesh_name, 
                                            num_components, field_metadata);
  STOP_FUNCTION_TIMER();
  return field;
}

bool silo_file_contains_edge_field(silo_file_t* file, const char* field_name, const char* mesh_name)
{
  START_FUNCTION_TIMER();
  bool result = silo_file_contains_mesh_field(file, MESH_EDGE, field_name, mesh_name);
  STOP_FUNCTION_TIMER();
  return result;
}

void silo_file_write_point_cloud(silo_file_t* file,
                                 const char* cloud_name,
                                 point_cloud_t* cloud)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_CLOBBER);

  // Point coordinates.
  int num_points = cloud->num_points;
  point_t* points = cloud->points;
  real_t* x = polymec_malloc(sizeof(real_t) * num_points);
  real_t* y = polymec_malloc(sizeof(real_t) * num_points);
  real_t* z = polymec_malloc(sizeof(real_t) * num_points);
  for (int i = 0; i < num_points; ++i)
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
  DBPutPointmesh(file->dbfile, (char*)cloud_name, 3, coords, num_points, SILO_FLOAT_TYPE, NULL); 
  polymec_free(x);
  polymec_free(y);
  polymec_free(z);

  // Write out the number of points to a special variable.
  char num_points_var[FILENAME_MAX];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", cloud_name);
  int one = 1;
  DBWrite(file->dbfile, num_points_var, &num_points, &one, 1, DB_INT);
  
  // Write out tag information.
  {
    char tag_name[FILENAME_MAX];
    snprintf(tag_name, FILENAME_MAX, "%s_node_tags", cloud_name);
    silo_file_write_tags(file, cloud->tags, tag_name);
  }

#if POLYMEC_HAVE_MPI
  if (file->nproc > 1)
    silo_file_add_subdomain_mesh(file, cloud_name, DB_POINTMESH, NULL);
#endif

  STOP_FUNCTION_TIMER();
}

point_cloud_t* silo_file_read_point_cloud(silo_file_t* file,
                                          const char* cloud_name)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_READ);

  // How many points does our cloud have?
  int num_points;
  char num_points_var[FILENAME_MAX];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", cloud_name);
  ASSERT(DBInqVarExists(file->dbfile, num_points_var));
  DBReadVar(file->dbfile, num_points_var, &num_points);

  DBpointmesh* pm = DBGetPointmesh(file->dbfile, (char*)cloud_name);
  if (pm == NULL)
    polymec_error("Point mesh '%s' was not found in the Silo file.", cloud_name);
  point_t* points = polymec_malloc(sizeof(point_t) * num_points);
  double* x = pm->coords[0];
  double* y = pm->coords[1];
  double* z = pm->coords[2];
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
    char tag_name[FILENAME_MAX];
    snprintf(tag_name, FILENAME_MAX, "%s_node_tags", cloud_name);
    silo_file_read_tags(file, tag_name, cloud->tags);
  }

  STOP_FUNCTION_TIMER();
  return cloud;
}

bool silo_file_contains_point_cloud(silo_file_t* file, const char* cloud_name)
{
  return (DBInqVarExists(file->dbfile, cloud_name) && 
          (DBInqVarType(file->dbfile, cloud_name) == DB_POINTMESH));
}

void silo_file_write_scalar_point_field(silo_file_t* file,
                                        const char* field_name,
                                        const char* cloud_name,
                                        real_t* field_data,
                                        silo_field_metadata_t* field_metadata)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_CLOBBER);

  // How many points does our mesh have?
  char num_points_var[FILENAME_MAX];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", cloud_name);
  ASSERT(DBInqVarExists(file->dbfile, num_points_var));
  int num_points;
  DBReadVar(file->dbfile, num_points_var, &num_points);

  // Write the field.
  DBoptlist* optlist = optlist_from_metadata(field_metadata);
  DBPutPointvar1(file->dbfile, field_name, cloud_name, field_data, num_points, SILO_FLOAT_TYPE, NULL);
  optlist_free(optlist);

#if POLYMEC_HAVE_MPI
  if (file->nproc > 1)
    silo_file_add_subdomain_field(file, cloud_name, field_name, DB_POINTVAR, optlist);
#endif
  STOP_FUNCTION_TIMER();
}

real_t* silo_file_read_scalar_point_field(silo_file_t* file,
                                          const char* field_name,
                                          const char* cloud_name,
                                          silo_field_metadata_t* field_metadata)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_READ);

  DBmeshvar* var = DBGetPointvar(file->dbfile, (char*)field_name);
  if (var == NULL)
    polymec_error("Field '%s' was not found in the Silo file.", field_name);
  real_t* field = polymec_malloc(sizeof(real_t) * var->nels);
  memcpy(field, var->vals[0], sizeof(real_t) * var->nels);
  read_point_metadata(var, field_metadata);
  DBFreePointvar(var);
  STOP_FUNCTION_TIMER();
  return field;
}

void silo_file_write_point_field(silo_file_t* file,
                                 const char** field_component_names,
                                 const char* cloud_name,
                                 real_t* field_data,
                                 int num_components,
                                 silo_field_metadata_t** field_metadata)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_CLOBBER);

  for (int c = 0; c < num_components; ++c)
  {
    silo_field_metadata_t* metadata = (field_metadata != NULL) ? field_metadata[c] : NULL;
    silo_file_write_scalar_point_field(file, field_component_names[c], 
                                       cloud_name, field_data, metadata);
  }
  STOP_FUNCTION_TIMER();
}

real_t* silo_file_read_point_field(silo_file_t* file,
                                   const char** field_component_names,
                                   const char* cloud_name,
                                   int num_components,
                                   silo_field_metadata_t** field_metadata)
{
  START_FUNCTION_TIMER();
  ASSERT(file->mode == DB_READ);

  // How many points does our mesh have?
  char num_points_var[FILENAME_MAX];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", cloud_name);
  ASSERT(DBInqVarExists(file->dbfile, num_points_var));
  int num_points;
  DBReadVar(file->dbfile, num_points_var, &num_points);
  real_t* field = polymec_malloc(sizeof(real_t) * num_components * num_points); 
  for (int c = 0; c < num_components; ++c)
  {
    silo_field_metadata_t* metadata = (field_metadata != NULL) ? field_metadata[c] : NULL;
    real_t* comp_data = silo_file_read_scalar_cell_field(file, field_component_names[c], cloud_name, metadata);
    for (int i = 0; i < num_points; ++i)
      field[num_components*i+c] = comp_data[i];
    polymec_free(comp_data);
  }
  STOP_FUNCTION_TIMER();
  return field;
}

bool silo_file_contains_point_field(silo_file_t* file, 
                                    const char* field_name, 
                                    const char* cloud_name)
{
  bool result = (silo_file_contains_point_cloud(file, cloud_name) &&  // point cloud exists...
                 (DBInqVarType(file->dbfile, cloud_name) == DB_POINTMESH) &&  // cloud is a point cloud...
                 (DBInqVarExists(file->dbfile, field_name) &&  // field exists...
                  (DBInqVarType(file->dbfile, field_name) == DB_POINTVAR))); // field is actually a variable.

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

void silo_file_write_string(silo_file_t* file,
                            const char* string_name,
                            char* string_data)
{
  ASSERT(file->mode == DB_CLOBBER);
  ASSERT(string_data != NULL);

  int string_size = strlen(string_data);
  char size_name[FILENAME_MAX];
  snprintf(size_name, FILENAME_MAX, "%s_string_size", string_name);
  int one = 1;
  int result = DBWrite(file->dbfile, size_name, &string_size, &one, 1, DB_INT);
  if (result != 0)
    polymec_error("silo_file_write_string: write of string '%s' (length) failed.", string_name);
  if (string_size > 0)
  {
    char string_data_name[FILENAME_MAX];
    snprintf(string_data_name, FILENAME_MAX, "%s_string", string_name);
    result = DBWrite(file->dbfile, string_data_name, string_data, &string_size, 1, DB_CHAR);
    if (result != 0)
      polymec_error("silo_file_write_string: write of string '%s' failed.", string_name);
  }
}

char* silo_file_read_string(silo_file_t* file,
                            const char* string_name)
{
  ASSERT(file->mode == DB_READ);

  char string_data_name[FILENAME_MAX];
  snprintf(string_data_name, FILENAME_MAX, "%s_string", string_name);
  if (!DBInqVarExists(file->dbfile, string_data_name))
  {
    polymec_error("silo_file_read_string: Could not read string '%s'.", string_name);
    return NULL;
  }
  char size_name[FILENAME_MAX];
  snprintf(size_name, FILENAME_MAX, "%s_string_size", string_name);
  ASSERT(DBInqVarExists(file->dbfile, size_name));
  int string_size;
  DBReadVar(file->dbfile, size_name, &string_size);
  char* string = polymec_malloc(sizeof(char) * (string_size + 1));
  if (string_size > 0)
    DBReadVar(file->dbfile, string_data_name, string);
  string[string_size] = '\0';
  return string;
}

void silo_file_write_real_array(silo_file_t* file,
                                const char* array_name,
                                real_t* array_data,
                                int array_size)
{
  ASSERT(file->mode == DB_CLOBBER);
  ASSERT(array_data != NULL);
  ASSERT(array_size >= 0);

  char size_name[FILENAME_MAX];
  snprintf(size_name, FILENAME_MAX, "%s_real_array_size", array_name);
  int one = 1;
  int result = DBWrite(file->dbfile, size_name, &array_size, &one, 1, DB_INT);
  if (result != 0)
    polymec_error("silo_file_write_real_array: write of array '%s' failed.", array_name);
  if (array_size > 0)
  {
    char real_array_name[FILENAME_MAX];
    snprintf(real_array_name, FILENAME_MAX, "%s_real_array", array_name);
    result = DBWrite(file->dbfile, real_array_name, array_data, &array_size, 1, SILO_FLOAT_TYPE);
    if (result != 0)
      polymec_error("silo_file_write_real_array: write of array '%s' failed.", array_name);
  }
}

real_t* silo_file_read_real_array(silo_file_t* file,
                                  const char* array_name,
                                  int* array_size)
{
  ASSERT(file->mode == DB_READ);
  ASSERT(array_size != NULL);

  char real_array_name[FILENAME_MAX];
  snprintf(real_array_name, FILENAME_MAX, "%s_real_array", array_name);
  if (!DBInqVarExists(file->dbfile, real_array_name))
  {
    polymec_error("silo_file_read_real_array: Could not read array '%s'.", array_name);
    return NULL;
  }
  char size_name[FILENAME_MAX];
  snprintf(size_name, FILENAME_MAX, "%s_real_array_size", array_name);
  ASSERT(DBInqVarExists(file->dbfile, size_name));
  DBReadVar(file->dbfile, size_name, array_size);
  if (*array_size > 0)
  {
    real_t* array = polymec_malloc(sizeof(real_t) * *array_size);
    DBReadVar(file->dbfile, real_array_name, array);
    return array;
  }
  else
    return NULL;
}

void silo_file_write_int_array(silo_file_t* file,
                               const char* array_name,
                               int* array_data,
                               int array_size)
{
  ASSERT(file->mode == DB_CLOBBER);
  ASSERT(array_data != NULL);
  ASSERT(array_size >= 0);

  char size_name[FILENAME_MAX];
  snprintf(size_name, FILENAME_MAX, "%s_int_array_size", array_name);
  int one = 1;
  int result = DBWrite(file->dbfile, size_name, &array_size, &one, 1, DB_INT);
  if (result != 0)
    polymec_error("silo_file_write_int_array: write of array '%s' failed.", array_name);
  if (array_size > 0)
  {
    char int_array_name[FILENAME_MAX];
    snprintf(int_array_name, FILENAME_MAX, "%s_int_array", array_name);
    result = DBWrite(file->dbfile, int_array_name, array_data, &array_size, 1, DB_INT);
    if (result != 0)
      polymec_error("silo_file_write_int_array: write of array '%s' failed.", array_name);
  }
}

int* silo_file_read_int_array(silo_file_t* file,
                              const char* array_name,
                              int* array_size)
{
  ASSERT(file->mode == DB_READ);
  ASSERT(array_size != NULL);

  char int_array_name[FILENAME_MAX];
  snprintf(int_array_name, FILENAME_MAX, "%s_int_array", array_name);
  if (!DBInqVarExists(file->dbfile, int_array_name))
  {
    polymec_error("silo_file_read_int_array: Could not read array '%s'.", array_name);
    return NULL;
  }
  char size_name[FILENAME_MAX];
  snprintf(size_name, FILENAME_MAX, "%s_int_array_size", array_name);
  ASSERT(DBInqVarExists(file->dbfile, size_name));
  DBReadVar(file->dbfile, size_name, array_size);
  if (*array_size > 0)
  {
    int* array = polymec_malloc(sizeof(int) * *array_size);
    DBReadVar(file->dbfile, int_array_name, array);
    return array;
  }
  else
    return NULL;
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

