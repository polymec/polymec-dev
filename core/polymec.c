// Copyright (c) 2012-2015, Jeffrey N. Johnson
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
#include "arena/proto.h"
#include "arena/pool.h"
#include "core/polymec.h"
#include "core/polymec_version.h"
#include "core/options.h"
#include "core/array.h"
#include "core/timer.h"

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

// This flag tells us whether polymec is already initialized.
static bool polymec_initialized = false;

// This function must be called to take advantage of Jonathan Shewchuk's 
// exact geometric predicates.
extern void exactinit();

// Command line arguments (used for provenance information).
int polymec_argc = 0;
char** polymec_argv = NULL;

// Invocation string, time, directory.
char* polymec_invoc_str = NULL;
time_t polymec_invoc_time = 0;
char* polymec_invoc_dir = NULL;

// Extra provenance information.
string_array_t* polymec_extra_provenance = NULL;

// Error handler.
static polymec_error_handler_function error_handler = NULL;

#if POLYMEC_HAVE_MPI
// MPI error handler.
static MPI_Errhandler mpi_error_handler;
#endif

// Functions to call on exit.
typedef void (*at_exit_func)();
static at_exit_func _atexit_funcs[32];
static int _num_atexit_funcs = 0;

static void shutdown()
{
  log_debug("polymec: Shutting down...");
  ASSERT(polymec_initialized);

  // Kill command line arguments.
  free(polymec_invoc_str);
  free(polymec_invoc_dir);
  for (int i = 0; i < polymec_argc; ++i)
    free(polymec_argv[i]);
  free(polymec_argv);

  // Kill extra provenance data.
  if (polymec_extra_provenance != NULL)
    string_array_free(polymec_extra_provenance);

  // Stop timing and generate a report.
  polymec_timer_t* polymec_timer = polymec_timer_get("polymec");
  polymec_timer_stop(polymec_timer);
  polymec_timer_report();

  // Call shutdown functions.
  log_debug("polymec: Calling %d shutdown functions.", _num_atexit_funcs);
  for (int i = 0; i < _num_atexit_funcs; ++i)
    _atexit_funcs[i]();

  MPI_Finalize();
#ifndef NDEBUG
  polymec_disable_fpe();
#endif
}

#if POLYMEC_HAVE_MPI
// This MPI error handler can be used to intercept MPI errors.
static void mpi_fatal_error_handler(MPI_Comm* comm, int* error_code, ...)
{
  int rank;
  MPI_Comm_rank(*comm, &rank);
  int len;
  char error_string[1024];
  MPI_Error_string(*error_code, error_string, &len);
  ASSERT(len < 1024);
  polymec_error("%s on rank %d\n", error_string, rank);
}
#endif

// This sets the logging level if it is given as a command-line argument.
static void set_up_logging()
{
  options_t* opts = options_argv();
  bool free_logging = false;
  char* logging = options_value(opts, "logging");
  if (logging == NULL)
  {
    // Are we maybe in a test environment, in which the logging=xxx 
    // argument is the first one passed?
    char* arg = options_argument(opts, 1);
    if ((arg != NULL) && (strstr(arg, "logging") != NULL) && (strcasecmp(arg, "logging=") != 0))
    {
      int num_words;
      char** words = string_split(arg, "=", &num_words);
      if (num_words == 2)
      {
        ASSERT(strcasecmp(words[0], "logging") == 0);
        free_logging = true;
        logging = string_dup(words[1]);
      }
      for (int i = 0; i < num_words; ++i)
        string_free(words[i]);
      polymec_free(words);
    }
  }
  if (logging != NULL)
  {
    if (!strcasecmp(logging, "debug"))
      set_log_level(LOG_DEBUG);
    else if (!strcasecmp(logging, "detail"))
      set_log_level(LOG_DETAIL);
    else if (!strcasecmp(logging, "info"))
      set_log_level(LOG_INFO);
    else if (!strcasecmp(logging, "urgent"))
      set_log_level(LOG_URGENT);
    else if (!strcasecmp(logging, "off"))
      set_log_level(LOG_NONE);
    if (free_logging)
      string_free(logging);
  }
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
    if ((arg != NULL) && (strstr(arg, "pause") != NULL) && (strcasecmp(arg, "pause=") != 0))
    {
      int num_words;
      char** words = string_split(arg, "=", &num_words);
      if (num_words == 2)
      {
        ASSERT(strcasecmp(words[0], "pause") == 0);
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

  if (delay != NULL)
  {
    int secs = atoi((const char*)delay);
    if (secs <= 0)
      polymec_error("Cannot pause for a non-positive interval.");
    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (nprocs > 1)
    {
      log_urgent("Pausing for %d seconds. PIDS: ", secs);
      int pid = (int)getpid();
      char hostname[32];
      gethostname(hostname, 32);
      int pids[nprocs];
      char hostnames[32*nprocs];
      MPI_Gather(&pid, 1, MPI_INT, pids, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Gather(hostname, 32, MPI_CHAR, hostnames, 32, MPI_CHAR, 0, MPI_COMM_WORLD);
      if (rank == 0)
      {
        for (int p = 0; p < nprocs; ++p)
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

void polymec_init(int argc, char** argv)
{
  if (!polymec_initialized)
  {
    // Jot down the invocation time.
    polymec_invoc_time = time(NULL);

    // Jot down command line args (use regular malloc).
    polymec_argc = argc;
    polymec_argv = malloc(sizeof(char*) * argc);
    for (int i = 0; i < argc; ++i)
      polymec_argv[i] = string_dup(argv[i]);

    // Construct the invocation string.
    int invoc_len = 2;
    for (int i = 0; i < polymec_argc; ++i)
      invoc_len += 1 + strlen(polymec_argv[i]);
    polymec_invoc_str = malloc(sizeof(char) * invoc_len);
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
      MPI_Init(&argc, &argv);

#if POLYMEC_HAVE_MPI
    // Set up the MPI error handler.
    MPI_Comm_create_errhandler(mpi_fatal_error_handler, &mpi_error_handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, mpi_error_handler);
#endif

    // Start up the garbage collector.
    GC_INIT();

    // Initialize variables for exact arithmetic.
    exactinit();

    // Register a shutdown function.
    atexit(&shutdown);

    // Parse command line options.
    options_parse(argc, argv);

    // If we are asked to set a specific logging level, do so.
    set_up_logging();

    // Start timing the main program.
    polymec_timer_t* polymec_timer = polymec_timer_get("polymec");
    polymec_timer_start(polymec_timer);

    // Okay! We're initialized.
    polymec_initialized = true;

    // If we are asked to pause, do so.
    pause_if_requested();

#ifndef NDEBUG
    // By default, we enable floating point exceptions for debug builds.
    polymec_enable_fpe();
#endif
  }
}

// Default error handler.
static noreturn void default_error_handler(const char* message)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("%d: Fatal error: %s\n", rank, message);
#if POLYMEC_HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, -1);
#else
  exit(-1);
#endif
  polymec_unreachable();
}

#if POLYMEC_HAVE_MPI
// Here are the error handlers for the Scotch partitioning library.
void SCOTCH_errorPrint(const char* const errstr, ...)
{
  va_list argp;
  va_start(argp, errstr);
  polymec_error(errstr, argp);
  va_end(argp);
}

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
  fprintf(stderr, "%s\n", err);
#if USE_MPI
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

static void handle_fpe_signal(int signal)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  polymec_error("%d: Detected a floating point exception signal.", rank);
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
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    fprintf(stderr, "polymec: not implemented: %s\n", component);
  exit(-1);
  polymec_unreachable();
}

void polymec_atexit(void (*func)()) 
{
  ASSERT(_num_atexit_funcs <= 32);
  _atexit_funcs[_num_atexit_funcs++] = func;
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
      int len = strlen(str);
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
      int end = ftell(fp);
      rewind(fp);
      static const int input_len_limit = 10 * 1024 * 1024; // 10kb input limit.
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
        fread(buff, 1, MIN(1000, end-offset), fp);
        if (end-offset < 1000)
          buff[end-offset] = '\0';
        fprintf(stream, "%s", buff);
        offset += MIN(1000, end-offset);
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
#ifdef LINUX
  return sysconf(_SC_NPROCESSORS_ONLN);
#elif defined APPLE
  return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

