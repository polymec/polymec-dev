# Build a configuration header file from our options.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/polymec.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/polymec.h"
)

# Generate the thirdparty Lua module from the metadata in 3rdparty/thirdparty.lua.
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lua_thirdparty.c
                   COMMAND ${PROJECT_BINARY_DIR}/bin/lua ${PROJECT_SOURCE_DIR}/tools/generate_thirdparty_module.lua ${PROJECT_SOURCE_DIR}/3rdparty/thirdparty.lua ${CMAKE_CURRENT_BINARY_DIR}/lua_thirdparty.c
                   DEPENDS ${PROJECT_SOURCE_DIR}/3rdparty/thirdparty.lua;${PROJECT_SOURCE_DIR}/tools/generate_thirdparty_module.lua;lua_proj)

include(add_polymec_library)
add_polymec_library(polymec_core polymec.c arch.c allocators.c logging.c timer.c
                    memory_info.c point.c point2.c tensor2.c hilbert.c kd_tree.c
                    octree.c rng.c parallel_sort.c permutations.c adj_graph.c
                    sp_func.c st_func.c polynomial.c enumerable.c
                    least_squares.c options.c partitioning.c text_buffer.c
                    serializer.c string_utils.c file_utils.c array_utils.c
                    lookup1.c linear_algebra.c exchanger.c blob_exchanger.c
                    lua_types.c lua_core.c lua_replace_tostring.c lua_array.c
                    lua_ndarray.c lua_driver.c
                    ${CMAKE_CURRENT_BINARY_DIR}/lua_thirdparty.c)
set_target_properties(polymec_core PROPERTIES LINKER_LANGUAGE C)

# polymec.c depends on the creation of polymec_version.h, the creation
# of the high-level headers, and the building of 3rd-party libraries.
add_dependencies(polymec_core update_version_h generate_headers all_3rdparty_libs)

set(POLYMEC_LIBRARIES polymec_core;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_LIBRARIES polymec_core;${POLYMEC_LIBRARIES})
add_subdirectory(tests)

# Install library.
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS polymec_core DESTINATION lib)

  # Install headers, including generated files.
  file(GLOB HEADERS *.h)
  LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/polymec.h)
  LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/polymec_version.h)
  install(FILES ${HEADERS} DESTINATION include/polymec/core)
endif()
