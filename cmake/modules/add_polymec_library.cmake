function(add_polymec_library lib)
  add_library(${lib} ${ARGN})
  set_target_properties(${lib} PROPERTIES FOLDER Libraries)
  target_link_libraries(${lib} ${POLYMEC_LIBRARIES})
  add_dependencies(${lib} all_3rdparty_libs)
endfunction(add_polymec_library)

