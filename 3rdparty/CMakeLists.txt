include_directories(${PROJECT_BINARY_DIR}/include)

# We use the rule of thumb that the optimal number of build threads is 
# one more than the number of logical cores.
math(EXPR NUM_BUILD_THREADS "${NUMBER_OF_CORES} + 1")

# Build cmocka, a C unit testing library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libcmocka.a)
  message(STATUS "Configuring cmocka...")
  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/cmocka
                  RESULT_VARIABLE stat)
  set(CMOCKA_C_FLAGS ${CMAKE_C_FLAGS})
  string(REPLACE "-pedantic-errors" "-Wno-format-pedantic" CMOCKA_C_FLAGS ${CMOCKA_C_FLAGS})
  set(CMOCKA_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_C_FLAGS=${CMOCKA_C_FLAGS} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DWITH_STATIC_LIB=ON -DWITH_EXAMPLES=OFF)
  if (CMAKE_C_COMPILER_ID STREQUAL "Intel")
    set(CMOCKA_CMAKE_OPTS ${CMOCKA_CMAKE_OPTS} -DWORDS_BIGENDIAN=FALSE)
  endif()
  execute_process(COMMAND cmake ${CMOCKA_CMAKE_OPTS} ${CMAKE_CURRENT_SOURCE_DIR}/cmocka
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cmocka
                  OUTPUT_VARIABLE cmocka_config.log ERROR_VARIABLE cmocka_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cmocka_config.log ${cmocka_config_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Configuration of cmocka library failed:")
      message(FATAL_ERROR ${cmocka_config_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cmocka_config_errors.log ${cmocka_config_err})
      message(FATAL_ERROR "Configuration of cmocka library failed. See ${CMAKE_CURRENT_BINARY_DIR}/cmocka_config_*.log.")
    endif()
  endif()
  message(STATUS "Building and installing cmocka (${NUM_BUILD_THREADS} threads)...")
  execute_process(COMMAND make install -j${NUM_BUILD_THREADS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cmocka
                  OUTPUT_VARIABLE cmocka_build_log ERROR_VARIABLE cmocka_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cmocka_build.log ${cmocka_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Build of cmocka library failed:")
      message(FATAL_ERROR ${cmocka_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cmocka_build_errors.log ${cmocka_build_err})
      message(FATAL_ERROR "Build of cmocka library failed. See ${CMAKE_CURRENT_BINARY_DIR}/cmocka_build_*.log.")
    endif()
  endif()

  # On Macs we need to change the rpath of shared libraries so we can install them elsewhere.
  if (APPLE AND BUILD_SHARED_LIBS)
    execute_process(COMMAND install_name_tool -id @rpath/libcmocka${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/libcmocka${LIB_SUFFIX})
  endif()

endif()
add_library(cmocka STATIC IMPORTED GLOBAL)
set_target_properties(cmocka PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libcmocka.a)

# Build libgc, the C garbage collector library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libgc${LIB_SUFFIX})
  if (BUILD_SHARED_LIBS)
    set(SHARED_OPTS --enable-shared --disable-static)
  else()
    set(SHARED_OPTS --enable-static --disable-shared)
  endif()
  set(GC_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} ${SHARED_OPTS} --enable-threads=pthreads --disable-gcj-support --disable-dependency-tracking)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(GC_CONFIG_OPTS ${GC_CONFIG_OPTS} --enable-gc-debug)
  endif()

  message(STATUS "Preparing gc...")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/gc ${CMAKE_CURRENT_BINARY_DIR}/gc
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message(STATUS "Configuring gc...")
  execute_process(COMMAND env CC=${CMAKE_C_COMPILER} ./configure ${GC_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gc
                  OUTPUT_VARIABLE gc_config_log ERROR_VARIABLE gc_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/gc_config.log ${gc_config_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Configuration of gc library failed:")
      message(FATAL_ERROR ${gc_config_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/gc_config_errors.log ${gc_config_err})
      message(FATAL_ERROR "Configuration of gc library failed. See ${CMAKE_CURRENT_BINARY_DIR}/gc_config_*.log.")
    endif()
  endif()

  message(STATUS "Building and installing gc (${NUM_BUILD_THREADS} threads)...")
  execute_process(COMMAND make install -j${NUM_BUILD_THREADS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gc
                  OUTPUT_VARIABLE gc_build_log ERROR_VARIABLE gc_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/gc_build.log ${gc_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Build of gc library failed:")
      message(FATAL_ERROR ${gc_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/gc_build_errors.log ${gc_build_err})
      message(FATAL_ERROR "Build of gc library failed. See ${CMAKE_CURRENT_BINARY_DIR}/gc_build_*.log.")
    endif()
  endif()

  # On Macs we need to change the rpath of shared libraries so we can install them elsewhere.
  if (APPLE AND BUILD_SHARED_LIBS)
    execute_process(COMMAND install_name_tool -id @rpath/libgc${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/libgc${LIB_SUFFIX})
  endif()

endif()
set(POLYMEC_TP_LIBS gc;${POLYMEC_TP_LIBS})

# Install libarena -- A fast C arena/memory pool implementation.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/include/arena)

  message(STATUS "Preparing libarena...")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/libarena ${CMAKE_CURRENT_BINARY_DIR}/libarena
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message(STATUS "Copying libarena headers into place...")
  file(GLOB arena_includes "${CMAKE_CURRENT_BINARY_DIR}/libarena/src/*.h")
  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/include/arena)
  foreach(inc ${arena_includes})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${inc} ${PROJECT_BINARY_DIR}/include/arena)
  endforeach()

endif()
file(GLOB arena_sources "${CMAKE_CURRENT_BINARY_DIR}/libarena/src/*.c")
add_library(arena ${arena_sources})
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libarena${LIB_SUFFIX} DESTINATION lib)
endif()
set(POLYMEC_TP_LIBS arena;${POLYMEC_TP_LIBS})

# Build zlib for compression.
if (NOT EXISTS ${Z_LIBRARY})
  set(ZLIB_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --static)

  message(STATUS "Preparing zlib...")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/zlib ${CMAKE_CURRENT_BINARY_DIR}/zlib
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message(STATUS "Configuring zlib...")
  execute_process(COMMAND env CC=${CMAKE_C_COMPILER} CFLAGS=${CMAKE_C_FLAGS} ./configure ${ZLIB_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/zlib
                  OUTPUT_VARIABLE zlib_config_log ERROR_VARIABLE zlib_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_config.log ${zlib_config_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Configuration of zlib library failed:")
      message(FATAL_ERROR ${zlib_config_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_config_errors.log ${zlib_config_err})
      message(FATAL_ERROR "Configuration of zlib library failed. See ${CMAKE_CURRENT_BINARY_DIR}/zlib_config_*.log.")
    endif()
  endif()

  message(STATUS "Building and installing zlib (${NUM_BUILD_THREADS} threads)...")
  execute_process(COMMAND make install -j${NUM_BUILD_THREADS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/zlib
                  OUTPUT_VARIABLE zlib_build_log ERROR_VARIABLE zlib_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_build.log ${zlib_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Build of zlib library failed:")
      message(FATAL_ERROR ${zlib_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_build_errors.log ${zlib_build_err})
      message(FATAL_ERROR "Build of zlib library failed. See ${CMAKE_CURRENT_BINARY_DIR}/zlib_build_*.log.")
    endif()
  endif()

endif()
set(POLYMEC_TP_LIBS z;${POLYMEC_TP_LIBS})

# Build the HDF5 parallel I/O library. 
if (${SILO_LIBRARY} MATCHES "h5")
  if (NOT EXISTS ${HDF5_LIBRARY})
    message(STATUS "Configuring hdf5...")
    set(HDF5_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER})
    set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}) 
    set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DZLIB_INCLUDE_DIRS=${PROJECT_BINARY_DIR}/include -DZLIB_LIBARIES=z -DZLIB_USE_EXTERNAL=OFF -DHDF5_ENABLE_PARALLEL=${HAVE_MPI})
    set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DHDF5_BUILD_CPP_LIB=OFF)
    set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DHDF5_BUILD_FORTRAN=OFF)
    set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DBUILD_TESTING=OFF -DHDF5_BUILD_EXAMPLES=OFF)
    if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/hdf5)
      execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/hdf5/CMakeCache.txt)
    else()
      execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/hdf5)
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR}/hdf5 ${HDF5_CMAKE_OPTS}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hdf5
                    OUTPUT_VARIABLE hdf5_config_log ERROR_VARIABLE hdf5_config_err
                    RESULT_VARIABLE stat)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_config.log ${hdf5_config_log})
    if (NOT stat EQUAL 0)
      if (TRAVIS_CI)
        message(WARNING "Configuration of hdf5 library failed:")
        message(FATAL_ERROR ${hdf5_config_err})
      else()
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_config_errors.log ${hdf5_config_err})
        message(FATAL_ERROR "Configuration of hdf5 library failed. See ${CMAKE_CURRENT_BINARY_DIR}/hdf5_config_*.log.")
      endif()
    endif()

    message(STATUS "Building and installing hdf5 (${NUM_BUILD_THREADS} threads)...")
    execute_process(COMMAND make install -j${NUM_BUILD_THREADS}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hdf5
                    OUTPUT_VARIABLE hdf5_build_log ERROR_VARIABLE hdf5_build_err
                    RESULT_VARIABLE stat)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_build.log ${hdf5_build_log})
    if (NOT stat EQUAL 0)
      if (TRAVIS_CI)
        message(WARNING "Build of hdf5 library failed:")
        message(FATAL_ERROR ${hdf5_build_err})
      else()
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_build_errors.log ${hdf5_build_err})
        message(FATAL_ERROR "Build of hdf5 library failed. See ${CMAKE_CURRENT_BINARY_DIR}/hdf5_build_*.log.")
      endif()
    endif()

    # Installation considerations.
    if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
      # On Macs we need to change the rpath of shared libraries so we can install them elsewhere.
      if (APPLE AND BUILD_SHARED_LIBS)
        execute_process(COMMAND otool -DX ${PROJECT_BINARY_DIR}/lib/lib${HDF5_LIB_NAME}${LIB_SUFFIX} 
                        OUTPUT_VARIABLE libhdf5_id)
        string(STRIP ${libhdf5_id} libhdf5_id)
        execute_process(COMMAND install_name_tool -id @rpath/lib${HDF5_LIB_NAME}${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/${HDF5_LIB_NAME}${LIB_SUFFIX})
        execute_process(COMMAND install_name_tool -id @rpath/lib${HDF5_HL_LIB_NAME}${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/lib${HDF5_HL_LIB_NAME}${LIB_SUFFIX})
        execute_process(COMMAND install_name_tool -change ${libhdf5_id} @rpath/libhdf5${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/${HDF5_HL_LIB_NAME}${LIB_SUFFIX})
      endif()
    endif()

  endif()
  if (HDF5_LIBRARY STREQUAL "${PROJECT_BINARY_DIR}/lib/lib${HDF5_LIB_NAME}${LIB_SUFFIX}")
    set(POLYMEC_HDF5_LIBRARY "${CMAKE_INSTALL_PREFIX}/lib/lib${HDF5_LIB_NAME}${LIB_SUFFIX}" PARENT_SCOPE)
    set(POLYMEC_HDF5_HL_LIBRARY "${CMAKE_INSTALL_PREFIX}/lib/lib${HDF5_HL_LIB_NAME}${LIB_SUFFIX}" PARENT_SCOPE)
    set(POLYMEC_HDF5_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include" PARENT_SCOPE)
  else()
    set(POLYMEC_HDF5_LIBRARY ${HDF5_LIBRARY} PARENT_SCOPE)
    set(POLYMEC_HDF5_HL_LIBRARY ${HDF5_HL_LIBRARY} PARENT_SCOPE)
    set(POLYMEC_HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIR} PARENT_SCOPE)
  endif()

  # If our HDF5 library is named something weird, copy it to a file with the 
  # regular boring name. This prevents libraries like Silo and Exodus (which 
  # take a while to catch up to recent developments) from freaking out.
  if (NOT HDF5_LIB_NAME STREQUAL "hdf5" AND HDF5_LIBRARY MATCHES "${PROJECT_BINARY_DIR}/lib" AND NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libhdf5${LIB_SUFFIX})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${HDF5_LIBRARY} ${PROJECT_BINARY_DIR}/lib/libhdf5${LIB_SUFFIX})
  endif()

  set(POLYMEC_TP_LIBS ${HDF5_LIBRARIES};${POLYMEC_TP_LIBS})
endif()

# Install silo with HDF5 support.
if (NOT EXISTS ${SILO_LIBRARY})

  message(STATUS "Preparing silo...")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/silo ${CMAKE_CURRENT_BINARY_DIR}/silo
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()
  execute_process(COMMAND patch -p1 -i ${CMAKE_CURRENT_SOURCE_DIR}/silo.patch 
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/silo
                  OUTPUT_VARIABLE patchy
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Patching library source failed.")
  endif()

  message(STATUS "Configuring silo...")
  set(SILO_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR})
  set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --with-hdf5=${HDF5_INCLUDE_DIR},${HDF5_LIBRARY_DIR})
  set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --with-zlib=${Z_INCLUDE_DIR},${Z_LIBRARY_DIR})
  set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --disable-fortran --without-readline --disable-silex)
  set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --disable-dependency-tracking)
  set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --disable-hzip --disable-fpzip)
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --enable-optimization)
    set(SILO_C_FLAGS -O2)
  else()
    set(SILO_C_FLAGS -g)
  endif()
  if (BUILD_SHARED_LIBS AND LINUX)
    set(SILO_C_FLAGS "${SILO_C_FLAGS} -fPIC")
  endif()
  set(SILO_ENV CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CFLAGS=${SILO_C_FLAGS})
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
    if (TRAVIS_CI)
      message(WARNING "Configuration of silo library failed:")
      message(FATAL_ERROR ${silo_config_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/silo_config_errors.log ${silo_config_err})
      message(FATAL_ERROR "Configuration of silo library failed. See ${CMAKE_CURRENT_BINARY_DIR}/silo_config_*.log.")
    endif()
  endif()

  message(STATUS "Building and installing silo (${NUM_BUILD_THREADS} threads)...")
  execute_process(COMMAND make install -j${NUM_BUILD_THREADS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/silo
                  OUTPUT_VARIABLE silo_build_log ERROR_VARIABLE silo_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/silo_build.log ${silo_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Build of silo library failed:")
      message(FATAL_ERROR ${silo_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/silo_build_errors.log ${silo_build_err})
      message(FATAL_ERROR "Build of silo library failed. See ${CMAKE_CURRENT_BINARY_DIR}/silo_build_*.log.")
    endif()
  endif()

endif()
set(POLYMEC_TP_LIBS ${SILO_LIBRARIES};${POLYMEC_TP_LIBS})

# Build LIS, a linear iterative solver library that allows us to bundle a 
# built-in Krylov solver without a YUGE library like PETSc or HYPRE.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/liblis.a)
  if (BUILD_SHARED_LIBS)
    set(SHARED_OPTS --enable-shared --disable-static)
  else()
    set(SHARED_OPTS --enable-static --disable-shared)
  endif()
  set(LIS_CONFIG_OPTS --enable-longlong ${SHARED_OPTS} --prefix=${PROJECT_BINARY_DIR})
  if (HAVE_OPENMP)
    set(LIS_CONFIG_OPTS ${LIS_CONFIG_OPTS} --enable-omp)
  endif()
  if (HAVE_MPI)
    set(LIS_CONFIG_OPTS ${LIS_CONFIG_OPTS} --enable-mpi)
  endif()
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(LIS_CONFIG_OPTS ${LIS_CONFIG_OPTS} --enable-debug)
  endif()
  
  message(STATUS "Preparing lis...")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/lis ${CMAKE_CURRENT_BINARY_DIR}/lis
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Adjusting lis configure script for Debug build...")
    file(READ ${CMAKE_CURRENT_BINARY_DIR}/lis/configure lis_configure)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lis/configure.old "${lis_configure}")
    string(REPLACE "-D_DEBUG" "" lis_configure "${lis_configure}")
    string(REPLACE "-O3" "-g" lis_configure "${lis_configure}")
    string(REPLACE "-fomit-frame-pointer" "-fno-omit-frame-pointer" lis_configure "${lis_configure}")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lis/configure "${lis_configure}")
  endif()

  message(STATUS "Configuring lis...")
  if (HAVE_MPI)
    set(LIS_ENV MPICC=${CMAKE_C_COMPILER})
  else()
    set(LIS_ENV CC=${CMAKE_C_COMPILER})
  endif()
  set(LIS_CFLAGS ${CMAKE_C_FLAGS})
  string(REPLACE "-Werror-implicit-function-declaration" "" LIS_CFLAGS ${LIS_CFLAGS})
  string(REPLACE "-pedantic-errors" "" LIS_CFLAGS ${LIS_CFLAGS})
  set(LIS_ENV ${LIS_ENV} CFLAGS=${LIS_CFLAGS})
  execute_process(COMMAND env ${LIS_ENV} ./configure ${LIS_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lis
                  OUTPUT_VARIABLE lis_config_log ERROR_VARIABLE lis_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lis_config.log ${lis_config_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Configuration of lis library failed:")
      message(FATAL_ERROR ${lis_config_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lis_config_errors.log ${lis_config_err})
      message(FATAL_ERROR "Configuration of lis library failed. See ${CMAKE_CURRENT_BINARY_DIR}/lis_config_*.log.")
    endif()
  endif()

  # NOTE: LIS's 'make install' cannot be run on more than 1 thread reliably,
  # NOTE: so we have to separate the build/install processes.
  message(STATUS "Building and installing lis (${NUM_BUILD_THREADS} threads)...")
  execute_process(COMMAND make -j${NUM_BUILD_THREADS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lis
                  OUTPUT_VARIABLE lis_build_log ERROR_VARIABLE lis_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lis_build.log ${lis_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Build of lis library failed:")
      message(FATAL_ERROR ${lis_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lis_build_errors.log ${lis_build_err})
      message(FATAL_ERROR "Build of lis library failed. See ${CMAKE_CURRENT_BINARY_DIR}/lis_build_*.log.")
    endif()
  endif()

  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lis
                  OUTPUT_VARIABLE lis_install_log ERROR_VARIABLE lis_install_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lis_install.log ${lis_install_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Install of lis library failed:")
      message(FATAL_ERROR ${lis_install_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lis_install_errors.log ${lis_install_err})
      message(FATAL_ERROR "Install of lis library failed. See ${CMAKE_CURRENT_BINARY_DIR}/lis_install_*.log.")
    endif()
  endif()

endif()
set(POLYMEC_TP_LIBS lis;${POLYMEC_TP_LIBS})

# Build SuperLU, a (serial) direct sparse linear solver library that 
# we will use for preconditioning.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libsuperlu.a)
  message(STATUS "Preparing superlu...")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/superlu ${CMAKE_CURRENT_BINARY_DIR}/superlu
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message(STATUS "Adjusting superlu make.inc...")
  file(READ ${CMAKE_CURRENT_BINARY_DIR}/superlu/make.inc superlu_makeinc)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/superlu/make.inc.old ${superlu_makeinc})
  string(REPLACE "gcc" ${CMAKE_C_COMPILER} superlu_makeinc ${superlu_makeinc})
  string(REPLACE "$(HOME)/Dropbox/Codes/SuperLU/SuperLU" ${PROJECT_BINARY_DIR} superlu_makeinc ${superlu_makeinc})
  string(REPLACE "libsuperlu_5.0.a" "libsuperlu.a" superlu_makeinc ${superlu_makeinc})

  # SuperLU has some non-standard C, so we have to be careful about adjusting its 
  # C compiler flags. :-/
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
      set(superlu_debug_flags "-O0 -g3")
    else()
      set(superlu_debug_flags "-O0 -g")
    endif()
    string(REPLACE "-O3" "${superlu_debug_flags}" superlu_makeinc ${superlu_makeinc})
  endif()
  if (BUILD_SHARED_LIBS AND LINUX)
    string(REPLACE "-std=c99" "-std=c99 -fPIC" superlu_makeinc ${superlu_makeinc})
    string(REPLACE "NOOPTS       =" "NOOPTS       =-fPIC" superlu_makeinc ${superlu_makeinc})
  endif()
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/superlu/make.inc ${superlu_makeinc})

  message(STATUS "Building and installing superlu (1 thread)...")
  # Looks like SuperLU wants to be built serially, so no -jN here.
  execute_process(COMMAND make lib
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/superlu
                  OUTPUT_VARIABLE superlu_build_log ERROR_VARIABLE superlu_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/superlu_build.log ${superlu_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Build of superlu library failed:")
      message(FATAL_ERROR ${superlu_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/superlu_build_errors.log ${superlu_build_err})
      message(FATAL_ERROR "Build of superlu library failed. See ${CMAKE_CURRENT_BINARY_DIR}/superlu_build_*.log.")
    endif()
  endif()

#  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/SuperLU_4.3/libsuperlu.a DESTINATION ${PROJECT_BINARY_DIR}/lib)
  file(GLOB includes "${CMAKE_CURRENT_BINARY_DIR}/superlu/SRC/*.h")
  foreach(inc ${includes})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${inc} ${PROJECT_BINARY_DIR}/include)
  endforeach()
endif()
set(POLYMEC_TP_LIBS superlu;${POLYMEC_TP_LIBS})

# Build Sundials, a C library for integrating stiff, nonlinear ODEs and 
# differential-algebraic equations.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libsundials_cvode.a)
  message(STATUS "Configuring sundials...")
  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/sundials
                  RESULT_VARIABLE stat)
  set(SUNDIALS_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS} -DEXAMPLES_ENABLE=OFF -DFCMIX_ENABLE=OFF -DSUNDIALS_PRECISION=${POLYMEC_PRECISION})
  set(SUNDIALS_CMAKE_OPTS ${SUNDIALS_CMAKE_OPTS} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER})
  if (HAVE_MPI)
    set(SUNDIALS_CMAKE_OPTS ${SUNDIALS_CMAKE_OPTS} -DMPI_ENABLE=ON -DMPI_MPICC=${CMAKE_C_COMPILER})
  else()
    set(SUNDIALS_CMAKE_OPTS ${SUNDIALS_CMAKE_OPTS} -DMPI_ENABLE=OFF)
  endif()
  execute_process(COMMAND cmake ${SUNDIALS_CMAKE_OPTS} ${CMAKE_CURRENT_SOURCE_DIR}/sundials
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/sundials
                  OUTPUT_VARIABLE sundials_config_log ERROR_VARIABLE sundials_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/sundials_config.log ${sundials_config_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Configuration of sundials library failed:")
      message(FATAL_ERROR ${sundials_config_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/sundials_config_errors.log ${sundials_config_err})
      message(FATAL_ERROR "Configuration of sundials library failed. See ${CMAKE_CURRENT_BINARY_DIR}/sundials_config_*.log.")
    endif()
  endif()

  message(STATUS "Building and installing sundials (${NUM_BUILD_THREADS} threads)...")
  execute_process(COMMAND make install -j${NUM_BUILD_THREADS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/sundials
                  OUTPUT_VARIABLE sundials_build_log ERROR_VARIABLE sundials_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/sundials_build.log ${sundials_build_log})
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Build of sundials library failed.")
    if (TRAVIS_CI)
      message(WARNING "Build of sundials library failed:")
      message(FATAL_ERROR ${sundials_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/sundials_build_errors.log ${sundials_build_err})
      message(FATAL_ERROR "Build of sundials library failed. See ${CMAKE_CURRENT_BINARY_DIR}/sundials_build_*.log.")
    endif()
  endif()
endif()
if (NOT HAVE_MPI)
  set (NVEC_LIBS nvecserial)
else()
  set (NVEC_LIBS nvecparallel;nvecserial)
endif()
foreach(component arkode;cvode;cvodes;ida;idas;kinsol;${NVEC_LIBS})
  list(APPEND SUNDIALS_LIBS sundials_${component})
endforeach()
set(POLYMEC_TP_LIBS ${SUNDIALS_LIBS};${POLYMEC_TP_LIBS})

# Build lua, a simple interpreted language library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/liblua.a)
  if (${APPLE})
    set(LUA_ARCH macosx)
  else() 
    set(LUA_ARCH generic)
  endif()
  message(STATUS "Preparing lua...")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/lua ${CMAKE_CURRENT_BINARY_DIR}/lua
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message(STATUS "Adjusting lua Makefile...")
  file(READ ${CMAKE_CURRENT_BINARY_DIR}/lua/src/Makefile lua_makefile)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lua/src/Makefile.old ${lua_makefile})
  string(REPLACE "gcc" ${CMAKE_C_COMPILER} lua_makefile ${lua_makefile})
  string(REPLACE "CC=cc" "CC=${CMAKE_C_COMPILER}" lua_makefile ${lua_makefile})
  string(REPLACE "MYCFLAGS=" "MYCFLAGS=-fPIC ${CMAKE_C_FLAGS}" lua_makefile ${lua_makefile})
  string(REPLACE "MYLDFLAGS=" "MYLDFLAGS=${CMAKE_EXE_LINKER_FLAGS}" lua_makefile ${lua_makefile})
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    string(REPLACE "-O2 -Wall" "${CMAKE_C_FLAGS} -g3" lua_makefile ${lua_makefile})
  endif()
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lua/src/Makefile ${lua_makefile})

  message(STATUS "Building lua (${NUM_BUILD_THREADS} threads)...")
  execute_process(COMMAND make ${LUA_ARCH} CC=${CMAKE_C_COMPILER} CFLAGS=${CMAKE_C_FLAGS} LDFLAGS=${CMAKE_EXE_LINKER_FLAGS} -j${NUM_BUILD_THREADS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lua
                  OUTPUT_VARIABLE lua_build_log ERROR_VARIABLE lua_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lua_build.log ${lua_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Build of lua library failed:")
      message(FATAL_ERROR ${lua_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lua_build_errors.log ${lua_build_err})
      message(FATAL_ERROR "Build of lua library failed. See ${CMAKE_CURRENT_BINARY_DIR}/lua_build_*.log.")
    endif()
  endif()
endif()
file(GLOB lua_includes "${CMAKE_CURRENT_BINARY_DIR}/lua/src/*.h")
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/liblua.a)
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/lua/src/liblua.a DESTINATION ${PROJECT_BINARY_DIR}/lib)
  foreach(inc ${lua_includes})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${inc} ${PROJECT_BINARY_DIR}/include)
  endforeach()
endif()
set(POLYMEC_TP_LIBS lua;${POLYMEC_TP_LIBS})

# Build jrs_predicates, Jonathan Richard Shewchuk's (public domain) library 
# for arbitrary-precision geometric predicates. 
# Since random() is not part of the C standard, we can't lean on it in general.
add_library(jrs_predicates jrs_predicates/predicates.c)
if (LINUX EQUAL 1)
  # On Linux, rand() and random() are the same.
  set_target_properties(jrs_predicates PROPERTIES COMPILE_FLAGS "-Drandom=rand")
endif()
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS jrs_predicates DESTINATION lib)
endif()
set(POLYMEC_TP_LIBS jrs_predicates;${POLYMEC_TP_LIBS})

# If we're building a parallel code, build INRIA's Scotch library.
if (HAVE_MPI)
  if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libscotch.a)
    find_program(FLEX flex)
    find_program(BISON bison)

    # Scotch is a little funky -- it has empty structs when built without 
    # threading support (and uses non-standard threading mechanisms like 
    # pthread barriers when built with threading support!). So we have 
    # to be a little more accommodating, use the gnu99 standard instead 
    # of c11, and turn off some error checking. I'm not wild about this, 
    # but if Scotch does its job, it doesn't much matter.
    set(SCOTCH_C_FLAGS ${CMAKE_C_FLAGS})
    if (NOT CMAKE_C_COMPILER_ID STREQUAL "Intel")
      string(REPLACE "c11" "gnu99" SCOTCH_C_FLAGS ${SCOTCH_C_FLAGS})
    endif()
    string(REPLACE "-pedantic-errors" "" SCOTCH_C_FLAGS ${SCOTCH_C_FLAGS})
    string(REPLACE "-Werror-implicit-function-declaration" "" SCOTCH_C_FLAGS ${SCOTCH_C_FLAGS})

    # For now, we use the "old" timing mechanism to simplify builds using 
    # older GNU compilers.
    set(SCOTCH_C_FLAGS "${SCOTCH_C_FLAGS} -DCOMMON_TIMING_OLD")

    # Respect debug flags.
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set(SCOTCH_C_FLAGS "${SCOTCH_C_FLAGS} -g")
    endif()

    message(STATUS "Preparing scotch...")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/scotch ${CMAKE_CURRENT_BINARY_DIR}/scotch
                    RESULT_VARIABLE stat)
    if (NOT stat EQUAL 0)
      message(FATAL_ERROR "Copying library source failed.")
    endif()
    configure_file(
      "scotch.Makefile.inc.in"
      "${CMAKE_CURRENT_BINARY_DIR}/scotch/src/Makefile.inc"
    )

    message(STATUS "Building scotch (${NUM_BUILD_THREADS} threads)...")
    execute_process(COMMAND make ptscotch -j${NUM_BUILD_THREADS} 
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/scotch/src
                    OUTPUT_VARIABLE scotch_build_log ERROR_VARIABLE scotch_build_err
                    RESULT_VARIABLE stat)
    if (NOT stat EQUAL 0)
      if (TRAVIS_CI)
        message(WARNING "Build of scotch library failed:")
        message(FATAL_ERROR ${scotch_build_err})
      else()
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/scotch_build_errors.log ${scotch_build_err})
        message(FATAL_ERROR "Build of scotch library failed. See ${CMAKE_CURRENT_BINARY_DIR}/scotch_build_*.log.")
      endif()
    endif()
  endif()
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/scotch/include/scotch.h DESTINATION ${PROJECT_BINARY_DIR}/include)
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/scotch/include/ptscotch.h DESTINATION ${PROJECT_BINARY_DIR}/include)
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/scotch/lib/libscotch.a DESTINATION ${PROJECT_BINARY_DIR}/lib)
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/scotch/lib/libptscotch.a DESTINATION ${PROJECT_BINARY_DIR}/lib)
  set(POLYMEC_TP_LIBS ptscotch;scotch;${POLYMEC_TP_LIBS})
endif()

# Add all the libraries to the build system at large.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${POLYMEC_TP_C_FLAGS}" PARENT_SCOPE)
set(POLYMEC_LIBRARIES ${POLYMEC_TP_LIBS};${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_INCDIRS ${POLYMEC_INCDIRS};${POLYMEC_TP_INCDIRS} PARENT_SCOPE)

