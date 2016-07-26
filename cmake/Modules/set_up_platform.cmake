# This macro identifies compilers and third-party library needs 
# for particular hosts.
macro(set_up_platform)

  # Set library suffix based on whether we're building shared/static.
  # FIXME: We have to hack this together here, since CMAKE_SHARED_LIBRARY_SUFFIX
  # FIXME: isn't available before project() is called, which is when set_up_platform()
  # FIXME: is invoked. Gross.
  if (BUILD_SHARED_LIBS)
    if (APPLE)
      set(LIB_SUFFIX .dylib)
    elseif (WIN32)
      set(LIB_SUFFIX .dll)
    else()
      set(LIB_SUFFIX .so)
    endif()
  else()
    set(LIB_SUFFIX .a)
  endif()

  # Set defaults for the various third-party libraries. These defaults
  # are hardwired because the project can't have been defined before 
  # this macro is executed, and so PROJECT_BINARY_DIR is unavailable.
  set(Z_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/libz.a")
  set(Z_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
  get_filename_component(Z_LIBRARY_DIR ${Z_LIBRARY} DIRECTORY)
  set(HDF5_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/libhdf5${LIB_SUFFIX}")
  set(HDF5_HL_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/libhdf5_hl${LIB_SUFFIX}")
  set(HDF5_LIBRARIES hdf5_hl;hdf5)
  set(HDF5_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
  get_filename_component(HDF5_LIBRARY_DIR ${Z_LIBRARY} DIRECTORY)
  set(SILO_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/libsiloh5.a")
  set(SILO_LIBRARIES siloh5)
  if (APPLE)
    set(NEED_LAPACK FALSE)
  else()
    set(NEED_LAPACK TRUE)
  endif()

  # Certain tools (e.g. patch) require TMPDIR to be defined. If it is not, 
  # we do so here.
  set(TMPDIR_VAR $ENV{TMPDIR})
  if (NOT TMPDIR_VAR)
    # FIXME: Does this exist everywhere?
    set(ENV{TMPDIR} "/tmp")
  endif()

  # Now we check for a specific machine name.
  if (DEFINED POLYMEC_MACHINE)

    if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/machines)
      # Clone the polymec-machines repo into machines/.
      message(STATUS "Fetching machine configurations to ${CMAKE_CURRENT_SOURCE_DIR}/machines...")
      execute_process(COMMAND git clone https://github.com/polymec/polymec-machines ${CMAKE_CURRENT_SOURCE_DIR}/machines
                      OUTPUT_VARIABLE shhh ERROR_VARIABLE shhh RESULT_VARIABLE stat)
      if (NOT ${stat} EQUAL 0)
        message(FATAL_ERROR "Failed to retrieve Polymec machine files.")
      endif()
    endif()

    # See whether we have a file for the given machine.
    set(MACHINE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/machines/${POLYMEC_MACHINE}.cmake")
    if (NOT EXISTS ${MACHINE_FILE})
      message(FATAL_ERROR "Invalid machine name: ${POLYMEC_MACHINE}. See machines/ for available options.")
    else()
      message(STATUS "Selected machine '${POLYMEC_MACHINE}'.")
      message(STATUS "See ${MACHINE_FILE} for special build/install/test instructions.")
    endif()

    # Now read the machine file.
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/machines/")
    include(${POLYMEC_MACHINE})
  endif()

  # Based on the batch system we're using, we set up our parallel 
  # execution environment.
  if (BATCH_SYSTEM STREQUAL "slurm")
    find_package(SLURM)
    set(MPIEXEC ${SLURM_SRUN_COMMAND})
    set(MPIEXEC_NUMPROC_FLAG -n)
    if (DEFINED PROCS_PER_NODE)
      set(MPIEXEC_PREFLAGS --ntasks-per-node=${PROCS_PER_NODE})
    endif()
    set(MPIEXEC_POSTFLAGS )
  else()
    # Regular old MPI execution environment
    find_package(MPI)
  endif()

endmacro()
