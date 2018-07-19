include(ExternalProject)
include_directories(${PROJECT_BINARY_DIR}/include)

# We use the rule of thumb that the optimal number of build threads is 
# one more than the number of logical cores.
math(EXPR NUM_BUILD_THREADS "${NUMBER_OF_CORES} + 1")

# We don't hold third-party library code to the same standard of excellence  
# to which we aspire. :-)
string(REPLACE "-Werror " " " CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
string(REPLACE "-Weverything " " " CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
string(REPLACE "-Wmissing-prototypes " " " CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
string(REPLACE "-Wfloat-equal" " " CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
string(REPLACE "-fprofile-arcs" " " CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
string(REPLACE "-ftest-coverage" " " CMAKE_C_FLAGS ${CMAKE_C_FLAGS})

# Create canonical top-level directories in the build directory.
foreach(dir include;lib;bin;share)
  set(top_level_dir ${PROJECT_BINARY_DIR}/${dir})
  if (NOT EXISTS ${top_level_dir})
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${top_level_dir})
  endif()
endforeach()

#----------------------------------------------------------------------------
# Build linenoise -- A readline alternative that provides terminal history
# in about 1,000 lines of code.
#----------------------------------------------------------------------------
add_library(linenoise ${CMAKE_CURRENT_SOURCE_DIR}/linenoise/linenoise.c)
set_target_properties(linenoise PROPERTIES FOLDER Libraries)
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/linenoise/linenoise.h ${PROJECT_BINARY_DIR}/include)
install(TARGETS linenoise DESTINATION ${PROJECT_BINARY_DIR}/lib)
set(POLYMEC_TP_LIBS linenoise;${POLYMEC_TP_LIBS})

#----------------------------------------------------------------------------
# Build cmocka, a C unit testing library.
#----------------------------------------------------------------------------
set(CMOCKA_C_FLAGS ${CMAKE_C_FLAGS})
string(REPLACE "-pedantic-errors" "-Wno-format-pedantic" CMOCKA_C_FLAGS ${CMOCKA_C_FLAGS})
set(CMOCKA_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_C_FLAGS=${CMOCKA_C_FLAGS} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DWITH_STATIC_LIB=ON -DWITH_EXAMPLES=OFF)
if (CMAKE_C_COMPILER_ID STREQUAL "Intel" OR CMAKE_C_COMPILER_ID MATCHES "AppleClang")
  set(CMOCKA_CMAKE_OPTS ${CMOCKA_CMAKE_OPTS} -DHAVE_WORDS_BIGENDIAN=TRUE -DWORDS_BIGENDIAN=FALSE)
endif()
ExternalProject_Add(cmocka_proj
                    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/cmocka
                    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmocka
                    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/cmocka
                    INSTALL_DIR ${PROJECT_BINARY_DIR}
                    CMAKE_ARGS ${CMOCKA_CMAKE_OPTS}
                    LOG_CONFIGURE TRUE
                    LOG_BUILD TRUE
                    LOG_INSTALL TRUE)
add_library(cmocka STATIC IMPORTED GLOBAL)
add_dependencies(cmocka cmocka_proj)
set_target_properties(cmocka PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libcmocka.a)
set(POLYMEC_TP_LIBS cmocka;${POLYMEC_TP_LIBS})

#----------------------------------------------------------------------------
# Build libarena -- A fast C arena/memory pool implementation.
#----------------------------------------------------------------------------
file(GLOB arena_sources "${CMAKE_CURRENT_SOURCE_DIR}/libarena/src/*.c")
add_library(arena ${arena_sources})
if (CMAKE_C_COMPILER_ID MATCHES "Clang")
  set_target_properties(arena PROPERTIES COMPILE_OPTIONS "-Wno-everything")
endif()
set_target_properties(arena PROPERTIES FOLDER Libraries)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/include/arena)
file(GLOB arena_includes "${CMAKE_CURRENT_SOURCE_DIR}/libarena/src/*.h")
foreach(inc ${arena_includes})
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${inc} ${PROJECT_BINARY_DIR}/include/arena)
endforeach()
set(POLYMEC_TP_LIBS arena;${POLYMEC_TP_LIBS})
install(TARGETS arena DESTINATION lib)

#----------------------------------------------------------------------------
# Build zlib: simple compression.
#----------------------------------------------------------------------------
set(ZLIB_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --static)
ExternalProject_Add(zlib_proj
                    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/zlib
                    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/zlib
                    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/zlib
                    INSTALL_DIR ${PROJECT_BINARY_DIR}
                    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/zlib ${CMAKE_CURRENT_BINARY_DIR}/zlib 
                    COMMAND env CC=${CMAKE_C_COMPILER} CFLAGS=${CMAKE_C_FLAGS} ./configure ${ZLIB_CONFIG_OPTS}
                    LOG_CONFIGURE TRUE
                    BUILD_COMMAND make -j${NUM_BUILD_THREADS}
                    LOG_BUILD TRUE
                    INSTALL_COMMAND make install
                    LOG_INSTALL TRUE)
add_library(z STATIC IMPORTED GLOBAL)
add_dependencies(z zlib_proj)
set_target_properties(z PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libz.a)
set(POLYMEC_TP_LIBS z;${POLYMEC_TP_LIBS})

#----------------------------------------------------------------------------
# Build the HDF5 parallel I/O library. 
#----------------------------------------------------------------------------
if (NOT EXISTS ${HDF5_LIBRARY})
  set(HDF5_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER})
  set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}) 
  set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DZLIB_INCLUDE_DIRS=${PROJECT_BINARY_DIR}/include -DZLIB_LIBARIES=z -DZLIB_USE_EXTERNAL=OFF -DHDF5_ENABLE_PARALLEL=${HAVE_MPI})
  set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DHDF5_BUILD_CPP_LIB=OFF)
  set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DHDF5_BUILD_FORTRAN=OFF)
  set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DHDF5_BUILD_TOOLS=ON -DBUILD_TESTING=OFF -DHDF5_BUILD_EXAMPLES=OFF)
  ExternalProject_Add(hdf5_proj
                      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/hdf5
                      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/hdf5
                      BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/hdf5
                      INSTALL_DIR ${PROJECT_BINARY_DIR}
                      CMAKE_ARGS ${HDF5_CMAKE_OPTS}
                      DEPENDS zlib_proj
                      LOG_CONFIGURE TRUE
                      LOG_BUILD TRUE
                      LOG_INSTALL TRUE)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/share/cmake/hdf5 
          DESTINATION share/cmake/hdf5)
endif()
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(HDF5_LIBRARY_BASENAME "hdf5_debug")
  set(HDF5_HL_LIBRARY_BASENAME "hdf5_hl_debug")
else()
  set(HDF5_LIBRARY_BASENAME "hdf5")
  set(HDF5_HL_LIBRARY_BASENAME "hdf5_hl")
endif()
set(HDF5_LIBRARY "lib${HDF5_LIBRARY_BASENAME}${LIB_SUFFIX}")
set(HDF5_HL_LIBRARY "lib${HDF5_HL_LIBRARY_BASENAME}${LIB_SUFFIX}")
if (BUILD_SHARED_LIBS)
  add_library(${HDF5_LIBRARY_BASENAME} SHARED IMPORTED GLOBAL)
  add_library(${HDF5_HL_LIBRARY_BASENAME} SHARED IMPORTED GLOBAL)
else()
  add_library(${HDF5_LIBRARY_BASENAME} STATIC IMPORTED GLOBAL)
  add_library(${HDF5_HL_LIBRARY_BASENAME} STATIC IMPORTED GLOBAL)
endif()
add_dependencies(${HDF5_LIBRARY_BASENAME} hdf5_proj)
set_target_properties(${HDF5_LIBRARY_BASENAME} PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/${HDF5_LIBRARY})
add_dependencies(${HDF5_HL_LIBRARY_BASENAME} hdf5_proj)
set_target_properties(${HDF5_HL_LIBRARY_BASENAME} PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/${HDF5_HL_LIBRARY})

# On Macs we need to change the rpath of shared libraries so we can install them elsewhere.
# FIXME: Broken!
#if (APPLE AND BUILD_SHARED_LIBS)
#  execute_process(COMMAND otool -DX ${HDF5_LIBRARY_PATH} OUTPUT_VARIABLE libhdf5_id)
#  string(STRIP ${libhdf5_id} libhdf5_id)
#  execute_process(COMMAND install_name_tool -id @rpath/${HDF5_LIBRARY} ${PROJECT_BINARY_DIR}/lib/${HDF5_LIBRARY})
#  execute_process(COMMAND install_name_tool -id @rpath/${HDF5_HL_LIBRARY} ${PROJECT_BINARY_DIR}/lib/${HDF5_HL_LIBRARY})
#  execute_process(COMMAND install_name_tool -change ${libhdf5_id} @rpath/libhdf5${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/${HDF5_HL_LIBRARY})
#endif()

if (HDF5_LIBRARY STREQUAL "${PROJECT_BINARY_DIR}/lib/libhdf5${LIB_SUFFIX}")
  set(POLYMEC_HDF5_LIBRARY "${CMAKE_INSTALL_PREFIX}/lib/${HDF5_LIBRARY}" PARENT_SCOPE)
  set(POLYMEC_HDF5_HL_LIBRARY "${CMAKE_INSTALL_PREFIX}/lib/${HDF5_HL_LIBRARY}" PARENT_SCOPE)
  set(POLYMEC_HDF5_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include" PARENT_SCOPE)
else()
  set(POLYMEC_HDF5_LIBRARY ${HDF5_LIBRARY} PARENT_SCOPE)
  set(POLYMEC_HDF5_HL_LIBRARY ${HDF5_HL_LIBRARY} PARENT_SCOPE)
  set(POLYMEC_HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIR} PARENT_SCOPE)
endif()
set(POLYMEC_TP_LIBS ${HDF5_HL_LIBRARY_BASENAME};${HDF5_LIBRARY_BASENAME};${POLYMEC_TP_LIBS})

#----------------------------------------------------------------------------
# Install silo with HDF5 support.
#----------------------------------------------------------------------------
set(USILO_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER})
set(USILO_CMAKE_OPTS ${USILO_CMAKE_OPTS} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}) 
set(USILO_CMAKE_OPTS ${USILO_CMAKE_OPTS} -DZLIB_INCLUDE_DIR=${Z_INCLUDE_DIR} -DZLIB_LIBRARY_DIR=${Z_LIBRARY_DIR})
set(USILO_CMAKE_OPTS ${USILO_CMAKE_OPTS} -DHDF5_INCLUDE_DIR=${HDF5_INCLUDE_DIR} -DHDF5_LIBRARY_DIR=${HDF5_LIBRARY_DIR})
ExternalProject_Add(usilo_proj
                    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/usilo
                    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/usilo
                    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/usilo
                    INSTALL_DIR ${PROJECT_BINARY_DIR}
                    CMAKE_ARGS ${USILO_CMAKE_OPTS}
                    DEPENDS hdf5_proj zlib_proj
                    LOG_CONFIGURE TRUE
                    LOG_BUILD TRUE
                    LOG_INSTALL TRUE)
if (BUILD_SHARED_LIBS)
  add_library(usilo SHARED IMPORTED GLOBAL)
else()
  add_library(usilo STATIC IMPORTED GLOBAL)
endif()
add_dependencies(usilo usilo_proj)
set_target_properties(usilo PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libusilo${LIB_SUFFIX})
set(POLYMEC_TP_LIBS usilo;${POLYMEC_TP_LIBS})

#----------------------------------------------------------------------------
# Build Sundials, a C library for integrating stiff, nonlinear ODEs and 
# differential-algebraic equations.
#----------------------------------------------------------------------------
set(SUNDIALS_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS} -DEXAMPLES_ENABLE=OFF -DFCMIX_ENABLE=OFF -DSUNDIALS_PRECISION=${POLYMEC_PRECISION})
set(SUNDIALS_CMAKE_OPTS ${SUNDIALS_CMAKE_OPTS} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER})
if (HAVE_MPI)
  set(SUNDIALS_CMAKE_OPTS ${SUNDIALS_CMAKE_OPTS} -DMPI_ENABLE=ON -DMPI_MPICC=${CMAKE_C_COMPILER})
else()
  set(SUNDIALS_CMAKE_OPTS ${SUNDIALS_CMAKE_OPTS} -DMPI_ENABLE=OFF)
endif()
ExternalProject_Add(sundials_proj
                    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/sundials
                    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/sundials
                    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/sundials
                    INSTALL_DIR ${PROJECT_BINARY_DIR}
                    CMAKE_ARGS ${SUNDIALS_CMAKE_OPTS}
                    LOG_CONFIGURE TRUE
                    LOG_BUILD TRUE
                    LOG_INSTALL TRUE)
if (NOT HAVE_MPI)
  set (NVEC_LIBS nvecserial)
else()
  set (NVEC_LIBS nvecparallel;nvecserial)
endif()
foreach(component arkode;cvode;cvodes;ida;idas;kinsol;${NVEC_LIBS})
  set(lib sundials_${component})
  add_library(${lib} STATIC IMPORTED GLOBAL)
  add_dependencies(${lib} sundials_proj)
  set_target_properties(${lib} PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/lib${lib}.a)
  list(APPEND SUNDIALS_LIBS ${lib})
endforeach()
set(POLYMEC_TP_LIBS ${SUNDIALS_LIBS};${POLYMEC_TP_LIBS})

#----------------------------------------------------------------------------
# Build lua, a simple interpreted language library.
#----------------------------------------------------------------------------
set(LUA_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -DLUA_USE_POSIX")
ExternalProject_Add(lua_proj
                    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/lua
                    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lua
                    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/lua
                    INSTALL_DIR ${PROJECT_BINARY_DIR}
                    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/lua ${CMAKE_CURRENT_BINARY_DIR}/lua 
                    LOG_CONFIGURE TRUE
                    BUILD_COMMAND make generic CC=${CMAKE_C_COMPILER} CFLAGS=${LUA_C_FLAGS} LDFLAGS=${CMAKE_EXE_LINKER_FLAGS} -j${NUM_BUILD_THREADS}
                    DEPENDS linenoise
                    LOG_BUILD TRUE
                    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/lua/src/liblua.a ${PROJECT_BINARY_DIR}/lib
                    LOG_INSTALL TRUE)
file(GLOB lua_includes "${CMAKE_CURRENT_SOURCE_DIR}/lua/src/*.h")
foreach(inc ${lua_includes})
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${inc} ${PROJECT_BINARY_DIR}/include)
endforeach()
add_library(lua STATIC IMPORTED GLOBAL)
set_target_properties(lua PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/liblua.a)
add_dependencies(lua lua_proj)
set(POLYMEC_TP_LIBS lua;${POLYMEC_TP_LIBS})

#----------------------------------------------------------------------------
# Build jrs_predicates, Jonathan Richard Shewchuk's (public domain) library 
# for arbitrary-precision geometric predicates. 
#----------------------------------------------------------------------------
# Since random() is not part of the C standard, we can't lean on it in general.
add_library(jrs_predicates jrs_predicates/predicates.c)
if (LINUX)
  # On Linux, rand() and random() are the same.
  set_target_properties(jrs_predicates PROPERTIES COMPILE_FLAGS 
                        "-Drandom=rand ${CMAKE_C_FLAGS}")
else()
  set_target_properties(jrs_predicates PROPERTIES COMPILE_FLAGS 
                        "${CMAKE_C_FLAGS}")
endif()
set_target_properties(jrs_predicates PROPERTIES FOLDER Libraries)
set(POLYMEC_TP_LIBS jrs_predicates;${POLYMEC_TP_LIBS})
install(TARGETS jrs_predicates DESTINATION lib)

#----------------------------------------------------------------------------
# Build specfun, a Fortran library implementing special functions.
#----------------------------------------------------------------------------
add_library(specfun specfun/special_functions.f90)
set(POLYMEC_TP_LIBS specfun;${POLYMEC_TP_LIBS})
install(TARGETS specfun DESTINATION lib)

#----------------------------------------------------------------------------
# If we're using MPI, build INRIA's Scotch library.
#----------------------------------------------------------------------------
if (HAVE_MPI)
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

  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/scotch)
  configure_file(
    "scotch.Makefile.inc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/scotch/src/Makefile.inc"
  )
  ExternalProject_Add(scotch_proj
                      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/scotch
                      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/scotch
                      BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/scotch/src
                      INSTALL_DIR ${PROJECT_BINARY_DIR}
                      CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/scotch ${CMAKE_CURRENT_BINARY_DIR}/scotch 
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND make ptscotch -j${NUM_BUILD_THREADS}
                      LOG_BUILD TRUE
                      INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/scotch/include/scotch.h ${PROJECT_BINARY_DIR}/include
                      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/scotch/include/ptscotch.h ${PROJECT_BINARY_DIR}/include
                      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/scotch/lib/libscotch.a ${PROJECT_BINARY_DIR}/lib
                      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/scotch/lib/libptscotch.a ${PROJECT_BINARY_DIR}/lib
                      LOG_INSTALL TRUE)
  if (BUILD_SHARED_LIBS)
    add_library(scotch SHARED IMPORTED GLOBAL)
    add_library(ptscotch SHARED IMPORTED GLOBAL)
  else()
    add_library(scotch STATIC IMPORTED GLOBAL)
    add_library(ptscotch STATIC IMPORTED GLOBAL)
  endif()
  add_dependencies(scotch scotch_proj)
  set_target_properties(scotch PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libscotch.a)
  add_dependencies(ptscotch scotch_proj)
  set_target_properties(ptscotch PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libptscotch.a)
  set(POLYMEC_TP_LIBS ptscotch;scotch;${POLYMEC_TP_LIBS})
endif()

# This target expresses a single dependency on all third-party libs.
add_custom_target(all_3rdparty_libs ALL 
                  DEPENDS ${POLYMEC_TP_LIBS})

# Add all the libraries to the build system at large.
set(POLYMEC_LIBRARIES ${POLYMEC_TP_LIBS};${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_INCDIRS ${POLYMEC_INCDIRS};${POLYMEC_TP_INCDIRS} PARENT_SCOPE)

