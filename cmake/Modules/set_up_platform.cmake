# This macro identifies compilers and third-party library needs 
# for particular hosts.
macro(set_up_platform)

  # Set defaults for the various third-party libraries. These defaults
  # are hardwired because the project can't have been defined before 
  # this macro is executed, and so PROJECT_BINARY_DIR is unavailable.
  set(Z_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/libz.a")
  set(HDF5_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/libhdf5.a")
  set(SILO_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/libsiloh5.a")

  # Certain tools (e.g. patch) require TMPDIR to be defined. If it is not, 
  # we do so here.
  set(TMPDIR_VAR $ENV{TMPDIR})
  if (NOT TMPDIR_VAR)
    # FIXME: Does this exist everywhere?
    set(ENV{TMPDIR} "/tmp")
  endif()

  # Get the hostname for this machine. 
  site_name(HOSTNAME)

  if (HOSTNAME MATCHES "edison") # NERSC Edison
    # Edison likes Intel's compilers...
    # ...but Intel's compilers don't do C11.
    #set(CMAKE_C_COMPILER icc)
    set(CMAKE_C_COMPILER gcc)
    set(CMAKE_CXX_COMPILER icpc)
    set(CMAKE_Fortran_COMPILER ifort)

    # We expect the following libraries to be available.
    set(Z_LIBRARY /usr/lib64/libz.a)

    set(HDF5_LOC $ENV{HDF5_DIR})
    if (NOT HDF5_LOC)
      message(FATAL_ERROR "HDF5_DIR not found. Please load the cray-hdf5 module.")
    endif()
    include_directories(${HDF5_LOC}/include)
    set(HDF5_LIBRARY ${HDF5_LOC}/lib/libhdf5.a)

    set(SILO_LOC $ENV{SILO_DIR})
    if (NOT SILO_LOC)
      message(FATAL_ERROR "SILO_DIR not found. Please load the silo module.")
    endif()
    include_directories(${SILO_LOC}/include)
    set(SILO_LIBRARY ${SILO_LOC}/lib/libsiloh5.a)
  endif()

endmacro()
