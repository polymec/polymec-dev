function(add_polymec_library lib)
  add_library(${lib} ${ARGN})
  add_dependencies(${lib} all_3rdparty_libs)
  set_target_properties(${lib} PROPERTIES FOLDER Libraries)
  if (BUILD_SHARED_LIBS)
    target_link_libraries(${lib} ${POLYMEC_LIBRARIES})
  endif()
endfunction(add_polymec_library)

