include_directories(${PROJECT_BINARY_DIR}/include)

# Build libgc, the C garbage collector library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libgc.a)
  set(GC_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared --enable-threads=pthreads --disable-gcj-support)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(GC_CONFIG_OPTS ${GC_CONFIG_OPTS} --enable-gc-debug)
  endif()

  message("Preparing gc...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/gc ${CMAKE_CURRENT_BINARY_DIR}/gc
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message("Configuring gc...")
  execute_process(COMMAND env CC=${CC} ./configure ${GC_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gc
                  OUTPUT_VARIABLE gc_config_log ERROR_VARIABLE gc_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/gc_config.log ${gc_config_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/gc_config_errors.log ${gc_config_err})
    message(FATAL_ERROR "Configuration of gc library failed.")
  endif()

  message("Building and installing gc...")
  execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gc
                  OUTPUT_VARIABLE gc_build_log ERROR_VARIABLE gc_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/gc_build.log ${gc_build_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/gc_build_errors.log ${gc_build_err})
    message(FATAL_ERROR "Build of gc library failed.")
  endif()

endif()
file(GLOB_RECURSE headers "${PROJECT_BINARY_DIR}/include/gc/*.h")
install(FILES ${headers} DESTINATION include/polymec/gc)
install(FILES ${PROJECT_BINARY_DIR}/lib/libgc.a DESTINATION lib)
set(POLYMEC_TP_LIBS gc;${POLYMEC_TP_LIBS})

# Install libarena -- A fast C arena/memory pool implementation.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/include/arena)

  message("Preparing libarena...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/libarena ${CMAKE_CURRENT_BINARY_DIR}/libarena
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message("Copying libarena headers into place...")
  file(GLOB arena_includes "${CMAKE_CURRENT_BINARY_DIR}/libarena/src/*.h")
  execute_process(COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/include/arena)
  foreach(inc ${arena_includes})
    execute_process(COMMAND cmake -E copy ${inc} ${PROJECT_BINARY_DIR}/include/arena)
  endforeach()

endif()
file(GLOB headers "${PROJECT_BINARY_DIR}/include/arena/*.h")
foreach(inc ${headers})
  install(FILES ${inc} DESTINATION include/polymec/arena)
endforeach()
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libarena.a DESTINATION lib)
file(GLOB arena_sources "${CMAKE_CURRENT_BINARY_DIR}/libarena/src/*.c")
add_library(arena ${arena_sources})
set(POLYMEC_TP_LIBS arena;${POLYMEC_TP_LIBS})

# Build zlib for compression.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libz.a)
  set(ZLIB_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --static)

  message("Preparing zlib...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/zlib ${CMAKE_CURRENT_BINARY_DIR}/zlib
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message("Configuring zlib...")
  execute_process(COMMAND env CC=${CMAKE_C_COMPILER} ./configure ${ZLIB_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/zlib
                  OUTPUT_VARIABLE zlib_config_log ERROR_VARIABLE zlib_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_config.log ${zlib_config_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_config_errors.log ${zlib_config_err})
    message(FATAL_ERROR "Configuration of zlib library failed.")
  endif()

  message("Building and installing zlib...")
  execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/zlib
                  OUTPUT_VARIABLE zlib_build_log ERROR_VARIABLE zlib_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_build.log ${zlib_build_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_build_errors.log ${zlib_build_err})
    message(FATAL_ERROR "Build of zlib library failed.")
  endif()

endif()
foreach(header zconf.h;zlib.h)
  install(FILES ${PROJECT_BINARY_DIR}/include/${header} DESTINATION include/polymec)
endforeach()
install(FILES ${PROJECT_BINARY_DIR}/lib/libz.a DESTINATION lib)
set(POLYMEC_TP_LIBS z;${POLYMEC_TP_LIBS})

# Build the HDF5 parallel I/O library. 
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libhdf5.a)
  set(HDF5_CONFIG_OPTS --enable-static --disable-shared --with-zlib=${PROJECT_BINARY_DIR})
  set(HDF5_CONFIG_OPTS ${HDF5_CONFIG_OPTS} --prefix=${PROJECT_BINARY_DIR})
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(HDF5_CONFIG_OPTS ${HDF5_CONFIG_OPTS} --enable-debug)
  else()
    set(HDF5_CONFIG_OPTS ${HDF5_CONFIG_OPTS} --disable-debug --enable-production)
  endif()
  if (${HAVE_MPI} EQUAL 1)
    set(HDF5_CONFIG_OPTS ${HDF5_CONFIG_OPTS} --enable-parallel)
  endif()

  message("Preparing hdf5...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/hdf5 ${CMAKE_CURRENT_BINARY_DIR}/hdf5
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message("Configuring hdf5...")
  execute_process(COMMAND env CC=${CMAKE_C_COMPILER} FC=${MAKE_Fortran_COMPILER} ./configure ${HDF5_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hdf5
                  OUTPUT_VARIABLE hdf5_config_log ERROR_VARIABLE hdf5_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_config.log ${hdf5_config_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_config_errors.log ${hdf5_config_err})
    message(FATAL_ERROR "Configuration of hdf5 library failed.")
  endif()

  message("Building and installing hdf5...")
  execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hdf5
                  OUTPUT_VARIABLE hdf5_build_log ERROR_VARIABLE hdf5_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_build.log ${hdf5_build_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_build_errors.log ${hdf5_build_err})
    message(FATAL_ERROR "Build of hdf5 library failed.")
  endif()

endif()
file(GLOB headers "${PROJECT_BINARY_DIR}/include/H5*.h")
install(FILES ${headers} DESTINATION include/polymec)
install(FILES ${PROJECT_BINARY_DIR}/lib/libhdf5.a DESTINATION lib)
set(POLYMEC_TP_LIBS hdf5;${POLYMEC_TP_LIBS})

# Install silo with HDF5 support.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libsiloh5.a)

  message("Preparing silo...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/silo ${CMAKE_CURRENT_BINARY_DIR}/silo
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message("Configuring silo...")
  set(SILO_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR})
  set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --with-hdf5=${PROJECT_BINARY_DIR}/include,${PROJECT_BINARY_DIR}/lib)
  set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --disable-fortran --without-readline --disable-browser --disable-silex)
  set(SILO_ENV CC=${CC} CXX=${CXX})
  if (LINUX EQUAL 1)
   # Note: Linux needs -ldl when building Silo because of HDF5's new 
   # Note: dynamically-loaded filters.
   set(SILO_ENV ${SILO_ENV} LIBS=-ldl)
  endif()
  execute_process(COMMAND env ${SILO_ENV} ./configure ${SILO_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/silo
                  OUTPUT_VARIABLE silo_config_log ERROR_VARIABLE silo_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/silo_config.log ${silo_config_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/silo_config_errors.log ${silo_config_err})
    message(FATAL_ERROR "Configuration of silo library failed.")
  endif()

  message("Building and installing silo...")
  execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/silo
                  OUTPUT_VARIABLE silo_build_log ERROR_VARIABLE silo_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/silo_build.log ${silo_build_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/silo_build_errors.log ${silo_build_err})
    message(FATAL_ERROR "Build of silo library failed.")
  endif()

endif()
install(FILES ${PROJECT_BINARY_DIR}/include/silo.h DESTINATION include/polymec)
install(FILES ${PROJECT_BINARY_DIR}/lib/libsiloh5.a DESTINATION lib)
set(POLYMEC_TP_LIBS siloh5;${POLYMEC_TP_LIBS})

# Build SuperLU, a (serial) direct sparse linear solver library that 
# we will use for preconditioning.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libsuperlu.a)
  message("Preparing superlu...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/superlu ${CMAKE_CURRENT_BINARY_DIR}/superlu
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message("Adjusting superlu make.inc...")
  file(READ ${CMAKE_CURRENT_BINARY_DIR}/superlu/make.inc superlu_makeinc)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/superlu/make.inc.old ${superlu_makeinc})
  string(REPLACE "gcc" ${CMAKE_C_COMPILER} superlu_makeinc ${superlu_makeinc})
  string(REPLACE "$(HOME)/Codes/SuperLU/SuperLU_4.3" ${PROJECT_BINARY_DIR} superlu_makeinc ${superlu_makeinc})
  string(REPLACE "libsuperlu_4.3.a" "libsuperlu.a" superlu_makeinc ${superlu_makeinc})
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
      set(superlu_debug_flags "-O0 -g3")
    else()
      set(superlu_debug_flags "-O0 -g")
    endif()
    string(REPLACE "-O3" "${superlu_debug_flags}" superlu_makeinc ${superlu_makeinc})
  endif()
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/superlu/make.inc ${superlu_makeinc})

  message("Building superlu...")
  # Looks like SuperLU wants to be built serially, so no -j4 here.
  execute_process(COMMAND make lib
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/superlu
                  OUTPUT_VARIABLE superlu_build_log ERROR_VARIABLE superlu_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/superlu_build.log ${superlu_build_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/superlu_build_errors.log ${superlu_build_err})
    message(FATAL_ERROR "Build of superlu library failed.")
  endif()

#  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/SuperLU_4.3/libsuperlu.a DESTINATION ${PROJECT_BINARY_DIR}/lib)
  file(GLOB includes "${CMAKE_CURRENT_BINARY_DIR}/superlu/SRC/*.h")
  foreach(inc ${includes})
    execute_process(COMMAND cmake -E copy ${inc} ${PROJECT_BINARY_DIR}/include)
  endforeach()
endif()
install(FILES ${PROJECT_BINARY_DIR}/lib/libsuperlu.a DESTINATION lib)
foreach(inc ${includes})
  install(FILES ${inc} DESTINATION include/polymec/)
endforeach()
set(POLYMEC_TP_LIBS superlu;${POLYMEC_TP_LIBS})

# Build Sundials, a C library for integrating stiff, nonlinear ODEs and 
# differential-algebraic equations.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libsundials_cvode.a)
  set(SUNDIALS_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --with-precision=${POLYMEC_PRECISION} --enable-static --disable-shared --disable-fcmix)
  set(SUNDIALS_CONFIG_OPTS ${SUNDIALS_CONFIG_OPTS})
  if (HAVE_MPI EQUAL 0)
    set(SUNDIALS_CONFIG_OPTS ${SUNDIALS_CONFIG_OPTS} --disable-mpi)
  endif()
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(SUNDIALS_CONFIG_OPTS ${SUNDIALS_CONFIG_OPTS} --with-cflags=-g)
  elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(SUNDIALS_CONFIG_OPTS ${SUNDIALS_CONFIG_OPTS} --with-cflags=-O3)
  endif()
  message("Preparing sundials...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/sundials ${CMAKE_CURRENT_BINARY_DIR}/sundials
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message("Configuring sundials...")
  execute_process(COMMAND env CC=${CC} ./configure ${SUNDIALS_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/sundials
                  OUTPUT_VARIABLE sundials_config_log ERROR_VARIABLE sundials_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/sundials_config.log ${sundials_config_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/sundials_config_errors.log ${sundials_config_err})
    message(FATAL_ERROR "Configuration of sundials library failed.")
  endif()

  message("Building and installing sundials...")
  execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/sundials
                  OUTPUT_VARIABLE sundials_build_log ERROR_VARIABLE sundials_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/sundials_build.log ${sundials_build_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/sundials_build_errors.log ${sundials_build_err})
    message(FATAL_ERROR "Build of sundials library failed.")
  endif()
endif()
if (HAVE_MPI EQUAL 0)
  set (NVEC_LIBS nvecserial)
else()
  set (NVEC_LIBS nvecparallel;nvecserial)
endif()
foreach(header_dir sundials;nvector;cvode;cvodes;ida;idas;kinsol)
  file(GLOB_RECURSE headers ${PROJECT_BINARY_DIR}/include/${header_dir}/*.h)
  install(FILES ${headers} DESTINATION include/polymec/${header_dir})
endforeach()
foreach(component cvode;cvodes;ida;idas;kinsol;${NVEC_LIBS})
  list(APPEND SUNDIALS_LIBS sundials_${component})
  install(FILES ${PROJECT_BINARY_DIR}/lib/libsundials_${component}.a DESTINATION lib)
endforeach()
set(POLYMEC_TP_LIBS ${SUNDIALS_LIBS};${POLYMEC_TP_LIBS})

# Build lua, a simple interpreted language library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/liblua.a)
  if (${APPLE})
    set(LUA_ARCH macosx)
  else() 
    set(LUA_ARCH ansi)
  endif()
  message("Preparing lua...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/lua ${CMAKE_CURRENT_BINARY_DIR}/lua
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message("Adjusting lua Makefile...")
  file(READ ${CMAKE_CURRENT_BINARY_DIR}/lua/src/Makefile lua_makefile)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lua/src/Makefile.old ${lua_makefile})
  string(REPLACE "gcc" ${CMAKE_C_COMPILER} lua_makefile ${lua_makefile})
  string(REPLACE "MYCFLAGS=" "MYCFLAGS=-fPIC ${CMAKE_C_FLAGS}" lua_makefile ${lua_makefile})
  string(REPLACE "MYLDFLAGS=" "MYLDFLAGS=${CMAKE_EXE_LINKER_FLAGS}" lua_makefile ${lua_makefile})
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    string(REPLACE "-O2 -Wall" "${CMAKE_C_FLAGS} -g3" lua_makefile ${lua_makefile})
  endif()
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lua/src/Makefile ${lua_makefile})

  message("Building lua...")
  execute_process(COMMAND make ${LUA_ARCH}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lua
                  OUTPUT_VARIABLE lua_build_log ERROR_VARIABLE lua_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lua_build.log ${lua_build_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lua_build_errors.log ${lua_build_err})
    message(FATAL_ERROR "Build of lua library failed.")
  endif()

  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/lua/src/liblua.a DESTINATION ${PROJECT_BINARY_DIR}/lib)
  file(GLOB includes "${CMAKE_CURRENT_BINARY_DIR}/lua/src/*.h")
  foreach(inc ${includes})
    execute_process(COMMAND cmake -E copy ${inc} ${PROJECT_BINARY_DIR}/include)
  endforeach()
endif()
install(FILES ${PROJECT_BINARY_DIR}/lib/liblua.a DESTINATION lib)
foreach(inc ${includes})
  install(FILES ${inc} DESTINATION include/polymec)
endforeach()
set(POLYMEC_TP_LIBS lua;${POLYMEC_TP_LIBS})
#set(POLYMEC_TP_INCDIRS ${POLYMEC_TP_INCDIRS};${PROJECT_BINARY_DIR}/include)

# Build jrs_predicates, Jonathan Richard Shewchuk's (public domain) library 
# for arbitrary-precision geometric predicates. 
# Since random() is not part of the C standard, we can't lean on it in general.
add_library(jrs_predicates jrs_predicates/predicates.c)
if (LINUX EQUAL 1)
  # On Linux, rand() and random() are the same.
  set_target_properties(jrs_predicates PROPERTIES COMPILE_FLAGS "-Drandom=rand")
endif()
install(TARGETS jrs_predicates DESTINATION lib)
set(POLYMEC_TP_LIBS jrs_predicates;${POLYMEC_TP_LIBS})

# Add all the libraries to the build system at large.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${POLYMEC_TP_C_FLAGS}" PARENT_SCOPE)
set(POLYMEC_LIBS ${POLYMEC_TP_LIBS};${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_INCDIRS ${POLYMEC_INCDIRS};${POLYMEC_TP_INCDIRS} PARENT_SCOPE)

