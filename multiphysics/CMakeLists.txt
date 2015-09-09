add_library(polymec_multiphysics 
            physics_kernel.c physics_state.c physics_controller.c)

set(POLYMEC_LIBRARIES polymec_multiphysics;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_LIBRARIES polymec_multiphysics;${POLYMEC_LIBRARIES})
add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  # Install library.
  install(TARGETS polymec_multiphysics DESTINATION lib)

  # Install headers, including generated files.
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/polymec/multiphysics)
endif()
