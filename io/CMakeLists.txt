include(add_polymec_library)
add_polymec_library(polymec_io silo_file.c silo_file_unimesh.c
                    silo_file_colmesh.c lua_io.c)
add_dependencies(polymec_io all_3rdparty_libs)

set(POLYMEC_LIBRARIES polymec_io;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_LIBRARIES polymec_io;${POLYMEC_LIBRARIES})

add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS polymec_io DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/polymec/io)
endif()
