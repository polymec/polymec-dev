// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <unistd.h>
#include <time.h>
#include <stdarg.h>
#include <signal.h>
#include <gc/gc.h>
#include <dlfcn.h>
#include "arena/proto.h"
#include "arena/pool.h"
#include "silo.h"
#include "core/polymec.h"
#include "core/polymec_version.h"
#include "core/options.h"
#include "core/array.h"
#include "core/timer.h"
#include "core/memory_info.h"
#include "core/file_utils.h"
#include "core/string_utils.h"

#if POLYMEC_HAVE_OPENMP
#include <omp.h>
#endif

#if POLYMEC_HAVE_VALGRIND
#include "valgrind.h"
#endif

// Standard C support for floating point environment.
#ifdef LINUX
#define _GNU_SOURCE // for older GNU compilers
#define __USE_GNU   // for newer GNU compilers
#endif
#include <fenv.h>

#ifdef LINUX
// FE_INEXACT    inexact result
// FE_DIVBYZERO  division by zero
// FE_UNDERFLOW  result not representable due to underflow
// FE_OVERFLOW   result not representable due to overflow
// FE_INVALID    invalid operation
#endif

#ifdef APPLE
#include <xmmintrin.h>
// _MM_MASK_DIV_INEXACT     inexact result
// _MM_MASK_DIV_ZERO        division by zero
// _MM_MASK_DIV_OVERFLOW    result not reproducible due to overflow
// _MM_MASK_DIV_UNDERFLOW   result not reproducible due to underflow
// _MM_MASK_INVALID         invalid operation
// _MM_MASK_DENORM          denormalized number
#endif

// FIXME: This isn't a standard C function.
extern int gethostname(char *name, size_t len);

// This flag tells us whether polymec is already initialized.
static bool polymec_initialized = false;

// MPI_COMM_WORLD rank, nprocs.
static int world_rank = 0;
static int world_nprocs = 1;

// This function must be called to take advantage of Jonathan Shewchuk's 
// exact geometric predicates.
extern void exactinit();

// Command line arguments (used for provenance information).
static int polymec_argc = 0;
static char** polymec_argv = NULL;
static char polymec_exe_name[FILENAME_MAX+1];

// Invocation string, time, directory.
static char* polymec_invoc_str = NULL;
static time_t polymec_invoc_time = 0;
static char* polymec_invoc_dir = NULL;

// Extra provenance information.
static string_array_t* polymec_extra_provenance = NULL;

// Error handler.
static polymec_error_handler_function error_handler = NULL;

#if POLYMEC_HAVE_MPI
// MPI error handler.
static MPI_Errhandler mpi_error_handler;
#endif

// Are we responsible for finalizing MPI?
static bool polymec_initialized_mpi = false;

// Directories to search for dynamically loadable libs.
static string_array_t* _dl_paths = NULL;

// Functions to call after initialization.
typedef void (*atinit_func)(int argc, char** argv);
DEFINE_ARRAY(atinit_array, atinit_func)
static atinit_array_t* _atinit_funcs = NULL;

// Functions to call on exit.
typedef void (*atexit_func)();
DEFINE_ARRAY(atexit_array, atexit_func)
static atexit_array_t* _atexit_funcs = NULL;

static void shutdown()
{
  ASSERT(polymec_initialized);

  // Print memory information if we're debugging.
  if (log_level() == LOG_DEBUG)
  {
    memory_info_t meminfo;
    get_memory_info(&meminfo);
    log_debug("polymec: Shutting down...");
    log_debug("polymec: Final memory usage: %zu kB resident, %zu kB virtual.", 
              meminfo.process_resident_size, meminfo.process_virtual_size);
    log_debug("polymec: Peak memory usage: %zu kB.", 
              meminfo.process_peak_resident_size);
  }

  // Kill extra provenance data.
  if (polymec_extra_provenance != NULL)
    string_array_free(polymec_extra_provenance);

  // Stop timing and generate a report.
  polymec_timer_stop_all();
  polymec_timer_report();

  // Kill command line arguments.
  string_free(polymec_invoc_str);
  string_free(polymec_invoc_dir);
  for (int i = 0; i < polymec_argc; ++i)
    string_free(polymec_argv[i]);
  polymec_free(polymec_argv);

  // Finalize any remaining garbage-collected objects.
  if (GC_should_invoke_finalizers())
  {
    int num_finalized = GC_invoke_finalizers();
    log_debug("polymec: finalized %d garbage-collected objects.", num_finalized);
  }

  // Call shutdown functions.
  log_debug("polymec: Calling %d shutdown functions.", _atexit_funcs);
  for (size_t i = 0; i < _atexit_funcs->size; ++i)
    _atexit_funcs->data[i]();

  // Final clean up.
  string_array_free(_dl_paths);
  _dl_paths = NULL;
  atinit_array_free(_atinit_funcs);
  _atinit_funcs = NULL;
  atexit_array_free(_atexit_funcs);
  _atexit_funcs = NULL;

  // Shut down MPI if we're supposed to.
  if (polymec_initialized_mpi)
    MPI_Finalize();

#ifndef NDEBUG
  polymec_disable_fpe();
#endif
}

static noreturn void handle_sigint(int signal)
{
  log_detail("polymec: intercepted Ctrl-C: exiting.");
  polymec_disable_fpe();
  exit(0);
  polymec_unreachable();
}

#if POLYMEC_HAVE_MPI
// This MPI error handler can be used to intercept MPI errors.
static noreturn void mpi_fatal_error_handler(MPI_Comm* comm, int* error_code, ...)
{
  int len;
  char error_string[1024];
  MPI_Error_string(*error_code, error_string, &len);
  ASSERT(len < 1024);
  polymec_error("%s on rank %d\n", error_string, world_rank);
}
#endif

// This Silo error handler intercepts I/O-related errors.
static noreturn void handle_silo_error(char* message)
{
  polymec_error("%s: %s", DBErrFuncname(), message);
}

// This sets the logging level/mode if given as a command-line argument.
static void set_up_logging()
{
  options_t* opts = options_argv();
  bool free_logging = false, free_logging_mode = false;
  char* logging = options_value(opts, "logging");
  char* logging_mode = options_value(opts, "logging_mode");
  if (logging == NULL)
  {
    // Are we maybe in a test environment, in which the logging=xxx 
    // argument is the first one passed?
    char* arg = options_argument(opts, 1);
    if ((arg != NULL) && (strstr(arg, "logging") != NULL) && (string_casecmp(arg, "logging=") != 0))
    {
      int num_words;
      char** words = string_split(arg, "=", &num_words);
      if (num_words == 2)
      {
        if (string_casecmp(words[0], "logging_mode") == 0)
        {
          free_logging = true;
          logging = string_dup(words[1]);
        }
        else if (string_casecmp(words[0], "logging_mode") == 0)
        {
          free_logging_mode = true;
          logging_mode = string_dup(words[1]);
        }
      }
      for (int i = 0; i < num_words; ++i)
        string_free(words[i]);
      polymec_free(words);
    }
  }
  opts = NULL;
  if (logging != NULL)
  {
    if (!string_casecmp(logging, "debug"))
      set_log_level(LOG_DEBUG);
    else if (!string_casecmp(logging, "detail"))
      set_log_level(LOG_DETAIL);
    else if (!string_casecmp(logging, "info"))
      set_log_level(LOG_INFO);
    else if (!string_casecmp(logging, "urgent"))
      set_log_level(LOG_URGENT);
    else if (!string_casecmp(logging, "off"))
      set_log_level(LOG_NONE);
    if (free_logging)
      string_free(logging);
  }
  if (logging_mode != NULL)
  {
    if (!string_casecmp(logging_mode, "all"))
    {
      log_debug("polymec_init: logging all MPI ranks.");
      set_log_mode(LOG_TO_ALL_RANKS);
    }
    else if (!string_casecmp(logging_mode, "single"))
    {
      log_debug("polymec_init: logging MPI rank 0.");
      set_log_mode(LOG_TO_SINGLE_RANK);
    }
#if POLYMEC_HAVE_MPI
    else if (string_is_integer(logging_mode))
    {
      char* ptr;
      int p = (int)strtol(logging_mode, &ptr, 10);
      if ((p < 0) || (p >= world_nprocs))
        polymec_error("polymec_init: invalid logging rank: %d", p);
      log_debug("polymec_init: logging MPI rank %d.", p);
      set_log_mode(LOG_TO_SINGLE_RANK);
      set_log_mpi_rank(log_level(), p);
    }
#endif
    if (free_logging_mode)
      string_free(logging_mode);
  }
}

// This sets up numbers of threads, etc.
static void set_up_threads()
{
#if POLYMEC_HAVE_OPENMP
  log_debug("polymec: Max number of OpenMP threads available: %d", omp_get_max_threads());

  options_t* opts = options_argv();
  char* num_threads_str = options_value(opts, "num_threads");
  if ((num_threads_str != NULL) && string_is_number(num_threads_str))
  {
    int num_threads = atoi(num_threads_str);
    log_debug("polymec: Setting number of OpenMP threads to %d.", num_threads);
    omp_set_num_threads(num_threads);
  }
  opts = NULL;
#endif
}

// This somewhat delicate procedure implements a simple mechanism to pause 
// and allow a developer to attach a debugger.
static void pause_if_requested()
{
  options_t* opts = options_argv();
  char* delay = options_value(opts, "pause");
  bool free_delay = false;
  if (delay == NULL)
  {
    // Are we maybe in a test environment, in which the pause=xxx 
    // argument is the first one passed?
    char* arg = options_argument(opts, 1);
    if ((arg != NULL) && (strstr(arg, "pause") != NULL) && (string_casecmp(arg, "pause=") != 0))
    {
      int num_words;
      char** words = string_split(arg, "=", &num_words);
      if (num_words == 2)
      {
        ASSERT(string_casecmp(words[0], "pause") == 0);
        delay = words[1];
        free_delay = true;
      }
      for (int i = 0; i < num_words; ++i)
      {
        if (i != 1)
          string_free(words[i]);
      }
      polymec_free(words);
    }
  }
  opts = NULL;

  if (delay != NULL)
  {
    int secs = atoi((const char*)delay);
    if (secs <= 0)
      polymec_error("Cannot pause for a non-positive interval.");
    if (world_nprocs > 1)
    {
      log_urgent("Pausing for %d seconds. PIDS: ", secs);
      int pid = (int)getpid();
#if POLYMEC_HAVE_MPI
      char hostname[MPI_MAX_PROCESSOR_NAME];
      int name_len;
      MPI_Get_processor_name(hostname, &name_len);
#else
      char hostname[256];
      gethostname(hostname, 256);
#endif
      int pids[world_nprocs];
      char hostnames[32*world_nprocs];
      MPI_Gather(&pid, 1, MPI_INT, pids, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Gather(hostname, 32, MPI_CHAR, hostnames, 32, MPI_CHAR, 0, MPI_COMM_WORLD);
      if (world_rank == 0)
      {
        for (int p = 0; p < world_nprocs; ++p)
          log_urgent("Rank %d (%s): %d", p, &hostnames[32*p], pids[p]);
      }
    }
    else
    {
      int pid = (int)getpid();
      log_urgent("Pausing for %d seconds (PID = %d).", secs, pid);
    }
    sleep((unsigned)secs);
    if (free_delay)
      string_free(delay);
  }
}

void polymec_atinit(void (*func)(int argc, char** argv))
{
  ASSERT(!polymec_initialized);
  atinit_array_append(_atinit_funcs, func);
}

void polymec_init(int argc, char** argv)
{
  if (!polymec_initialized)
  {
    // Jot down the invocation time.
    polymec_invoc_time = time(NULL);

    // Set up initialization and exit function arrays.
    _atinit_funcs = atinit_array_new();
    _atexit_funcs = atexit_array_new();

    // Paths to search for dynamically loadable libraries.
    _dl_paths = string_array_new();

    // Jot down command line args.
    polymec_argc = argc;
    polymec_argv = polymec_malloc(sizeof(char*) * argc);
    for (int i = 0; i < argc; ++i)
      polymec_argv[i] = string_dup(argv[i]);

    // Figure out the name (and only the name) of the executable.
    char dir[FILENAME_MAX+1];
    parse_path(polymec_argv[0], dir, polymec_exe_name);

    // Construct the invocation string.
    int invoc_len = 2;
    for (int i = 0; i < polymec_argc; ++i)
      invoc_len += 1 + strlen(polymec_argv[i]);
    polymec_invoc_str = polymec_malloc(sizeof(char) * invoc_len);
    polymec_invoc_str[0] = '\0';
    for (int i = 0; i < polymec_argc; ++i)
    {
      char arg[strlen(polymec_argv[i])+2];
      if (i == polymec_argc-1)
        sprintf(arg, "%s", polymec_argv[i]); 
      else
        sprintf(arg, "%s ", polymec_argv[i]);
      strcat(polymec_invoc_str, arg);
    }
    char* pwd = getenv("PWD");
    if (pwd != NULL)
      polymec_invoc_dir = string_dup(pwd);
    else
      polymec_invoc_dir = string_dup("(unknown)");

    // Start up MPI if needed.
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
      MPI_Init(&argc, &argv);
      polymec_initialized_mpi = true;
    }

#if POLYMEC_HAVE_MPI
    // Cache our rank and number of processes.
    MPI_Comm_size(MPI_COMM_WORLD, &world_nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Set up the MPI error handler.
    MPI_Comm_create_errhandler(mpi_fatal_error_handler, &mpi_error_handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, mpi_error_handler);
#endif

    // Set up the Silo I/O error handler.
    DBShowErrors(DB_ALL, handle_silo_error);

    // Start up the garbage collector.
    GC_INIT();

    // Initialize variables for exact arithmetic.
    exactinit();

    // Register a shutdown function.
    atexit(&shutdown);

    // Parse command line options.
    options_parse(argc, argv);

    // If we are asked to set a specific logging level/mode, do so.
    set_up_logging();

    // If we are asked to set up threads specifically, do so.
    set_up_threads();

    // Start timing the main program.
    polymec_timer_t* polymec_timer = polymec_timer_get("polymec");
    polymec_timer_start(polymec_timer);

    // Set up any required signal handlers (if they haven't already been 
    // claimed).
    void (*prev_handler)(int) = signal(SIGINT, handle_sigint);
    if (prev_handler != NULL)
      signal(SIGINT, prev_handler);

    // Okay! We're initialized.
    polymec_initialized = true;

    // If we are asked to pause, do so.
    pause_if_requested();

#ifndef NDEBUG
    // By default, we enable floating point exceptions for debug builds.
    polymec_enable_fpe();
#endif

    // Call initialization functions.
    log_debug("polymec: Calling %d initialization functions.", _atinit_funcs->size);
    for (size_t i = 0; i < _atinit_funcs->size; ++i)
      _atinit_funcs->data[i](argc, argv);

    // Print memory information if we're debugging.
    if (log_level() == LOG_DEBUG)
    {
      memory_info_t meminfo;
      get_memory_info(&meminfo);
      log_debug("polymec: initialization complete (%zu kB resident, %zu kB virtual).", 
                meminfo.process_resident_size, meminfo.process_virtual_size);
    }
  }
}

// Default error handler.
static noreturn void default_error_handler(const char* message)
{
  printf("%d: Fatal error: %s\n", world_rank, message);
#if POLYMEC_HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, -1);
#else
  exit(-1);
#endif
  polymec_unreachable();
}

#if POLYMEC_HAVE_MPI
// Here are the error handlers for the Scotch partitioning library.
// These are not part of polymec's public API, but are used by other 
// parts of the library, so must be available externally.
noreturn void SCOTCH_errorPrint(const char* const errstr, ...);
noreturn void SCOTCH_errorPrint(const char* const errstr, ...)
{
  va_list argp;
  va_start(argp, errstr);
  polymec_error(errstr, argp);
  va_end(argp);
}

void SCOTCH_errorPrintW(const char* const errstr, ...);
void SCOTCH_errorPrintW(const char* const errstr, ...)
{
  va_list argp;
  va_start(argp, errstr);
  polymec_warn(errstr, argp);
  va_end(argp);
}
#endif

void polymec_abort(const char* message, ...)
{
  // Extract the variadic arguments and splat them into a string.
  char err[1024];
  va_list argp;
  va_start(argp, message);
  vsnprintf(err, 1024, message, argp);
  va_end(argp);

  // Abort.
  fprintf(stderr, "%d: %s\n", world_rank, err);
#if POLYMEC_HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, -1); 
#else
  abort(); 
#endif
  polymec_unreachable();
}

void polymec_error(const char* message, ...)
{
  // Set the default error handler if no handler is set.
  if (error_handler == NULL)
    error_handler = &default_error_handler;

  // Extract the variadic arguments and splat them into a string.
  char err[1024];
  va_list argp;
  va_start(argp, message);
  vsnprintf(err, 1024, message, argp);
  va_end(argp);

  // Call the handler.
  error_handler(err);

  // Make sure we don't return.
  exit(-1);
  polymec_unreachable();
}

void polymec_set_error_handler(polymec_error_handler_function handler)
{
  error_handler = handler;
}

void polymec_warn(const char* message, ...)
{
  // Extract the variadic arguments and splat them into a string.
  va_list argp;
  va_start(argp, message);
  vfprintf(stderr, message, argp);
  va_end(argp);
}

static noreturn void handle_fpe_signal(int signal)
{
  polymec_error("%d: Detected a floating point exception signal.", world_rank);
}

void polymec_enable_fpe()
{
  feclearexcept(FE_ALL_EXCEPT);
#ifdef LINUX
  int flags = FE_DIVBYZERO |
              FE_INVALID   |
//              FE_UNDERFLOW |
              FE_OVERFLOW;

  feenableexcept(flags);
#endif

#ifdef APPLE
  // Catch all the interesting ones.
  _MM_SET_EXCEPTION_MASK(_MM_MASK_INEXACT | _MM_MASK_UNDERFLOW);
#endif
  signal(SIGFPE, handle_fpe_signal);
  log_debug("polymec: Enabled floating point exception support.");
}

void polymec_disable_fpe()
{
  fesetenv(FE_DFL_ENV);
  signal(SIGFPE, SIG_DFL);
  log_debug("polymec: Disabled floating point exception support.");
}


// The following suspend/restore mechanism uses standard C99 floating point 
// kung fu.
static fenv_t fpe_env;

void polymec_suspend_fpe()
{
  // Hold exceptions till further notice.
  feholdexcept(&fpe_env);
}

void polymec_restore_fpe()
{
  // Clear all exception flags.
  feclearexcept(FE_ALL_EXCEPT);

  // Now restore the previous floating point environment.
  fesetenv(&fpe_env);
}

void polymec_not_implemented(const char* component)
{
  if (world_rank == 0)
    fprintf(stderr, "polymec: not implemented: %s\n", component);
  exit(-1);
  polymec_unreachable();
}

void polymec_atexit(void (*func)()) 
{
  atexit_array_append(_atexit_funcs, func);
}


void polymec_add_dl_path(const char* path)
{
  ASSERT(polymec_initialized);
  ASSERT(path != NULL);
  if (directory_exists(path))
    string_array_append(_dl_paths, (char*)path);
}

void* polymec_dlopen(const char* lib_name)
{
  void* lib = NULL;
  for (size_t i = 0; i < _dl_paths->size; ++i)
  {
    char full_path[FILENAME_MAX+1];
    snprintf(full_path, FILENAME_MAX, "%s/lib%s%s", _dl_paths->data[i], 
             lib_name, SHARED_LIBRARY_SUFFIX);
    if (file_exists(full_path))
      lib = dlopen(full_path, RTLD_NOW);
    if (lib != NULL)
      break;
  }
  return lib;
}

void polymec_version_fprintf(const char* exe_name, FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "%s v%s\n", exe_name, POLYMEC_VERSION);
}

void polymec_provenance_fprintf(FILE* stream)
{
  ASSERT(polymec_initialized);
  if (stream == NULL) return;

  fprintf(stream, "=======================================================================\n");
  fprintf(stream, "                                Provenance:\n");
  fprintf(stream, "=======================================================================\n");
  fprintf(stream, "Version: %s\n", POLYMEC_VERSION);
  fprintf(stream, "Invoked with: %s\n", polymec_invoc_str);
  fprintf(stream, "Invoked on: %s", ctime(&polymec_invoc_time)); // No \n because of ctime
  fprintf(stream, "Invocation dir: %s\n\n", polymec_invoc_dir);

  if (polymec_extra_provenance != NULL)
  {
    fprintf(stream, "Additional provenance information:\n");
    for (int i = 0; i < polymec_extra_provenance->size; ++i)
    {
      char* str = polymec_extra_provenance->data[i];
      size_t len = strlen(str);
      if (str[len-1] == '\n')
        fprintf(stream, "%s", str);
      else
        fprintf(stream, "%s\n", str);
    }
    fprintf(stream, "\n");
  }

  if (POLYMEC_NUM_GIT_DIFFS > 0)
  {
    fprintf(stream, "=======================================================================\n");
    fprintf(stream, "Modifications to revision:\n");
    for (int i = 0; i < POLYMEC_NUM_GIT_DIFFS; ++i)
      fprintf(stream, "%s", POLYMEC_GIT_DIFFS[i]);
    fprintf(stream, "\n\n");
  }

  // If we received an input script, write out its contents.
  options_t* options = options_argv();

  // Start going backwards through the arguments, throwing out all the ones
  // with equals signs.
  char* input = NULL;
  int arg = options_num_arguments(options) - 1;
  while ((input == NULL) && (arg > 0))
  {
    char* candidate = options_argument(options, arg);
    if (string_contains(candidate, "="))
      --arg;
    else 
    {
      char path[FILENAME_MAX];
      if (strcmp(polymec_invoc_dir, "(unknown)") != 0)
        snprintf(path, FILENAME_MAX, "%s/%s", polymec_invoc_dir, candidate);
      else
        strcpy(path, candidate);
      FILE* fp = fopen(candidate, "r");
      if (fp != NULL)
      {
        input = candidate;
        fclose(fp);
      }
      else
        --arg;
    }
  }
  if (input != NULL)
  {
    FILE* fp = fopen(input, "r");
    if (fp == NULL)
      fprintf(stream, "Invalid input specified.");
    else
    {
      char buff[1024];
      fseek(fp, 0L, SEEK_END);
      long end = ftell(fp);
      rewind(fp);
      static const long input_len_limit = 10 * 1024 * 1024; // 10kb input limit.
      bool truncated = false;
      if (input_len_limit < end)
      {
        truncated = true;
        end = input_len_limit; 
      }
      fprintf(stream, "=======================================================================\n");
      fprintf(stream, "Contents of input script:\n");
      int offset = 0;
      while (offset < end)
      {
        size_t bytes_read = fread(buff, 1, MIN(1000, end-offset), fp);
        if (bytes_read < 1000)
          buff[bytes_read] = '\0';
        fprintf(stream, "%s", buff);
        offset += MIN(1000, bytes_read);
      }
      if (truncated)
        fprintf(stream, "\n<<< truncated >>>\n");
      fclose(fp);
    }
    fprintf(stream, "\n");
  }
  options = NULL;

  fprintf(stream, "=======================================================================\n\n");
}

const char* polymec_executable_name()
{
  ASSERT(polymec_initialized);
  return (const char*)polymec_exe_name;
}

const char* polymec_invocation()
{
  ASSERT(polymec_initialized);
  return (const char*)polymec_invoc_str;
}

time_t polymec_invocation_time()
{
  return polymec_invoc_time;
}

void polymec_append_provenance_data(const char* provenance_data)
{
  if (polymec_extra_provenance == NULL)
    polymec_extra_provenance = string_array_new();
  string_array_append_with_dtor(polymec_extra_provenance, 
                                string_dup(provenance_data),
                                string_free);
}

int polymec_num_cores()
{
  return (int)sysconf(_SC_NPROCESSORS_ONLN);
}

bool polymec_running_in_valgrind()
{
#if POLYMEC_HAVE_VALGRIND
  return (RUNNING_ON_VALGRIND);
#else
  return false;
#endif
}
