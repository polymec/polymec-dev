function(add_polymec_executable_with_libs exe libs)
  add_executable(${exe} ${ARGN})
  set_target_properties(${exe} PROPERTIES FOLDER Executables)
  target_link_libraries(${exe} ${libs})
endfunction(add_polymec_executable_with_libs)

function(add_polymec_executable exe)
  add_polymec_executable_with_libs(${exe} "${POLYMEC_LIBRARIES}" ${ARGN})
endfunction(add_polymec_executable)

