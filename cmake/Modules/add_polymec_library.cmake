function(add_polymec_library lib)
  add_library(${lib} ${ARGN})
  if (BUILD_SHARED_LIBS)
    target_link_libraries(${lib} ${POLYMEC_LIBRARIES})
  endif()
endfunction(add_polymec_library)

