function(add_polymec_executable exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} ${POLYMEC_LIBRARIES})
endfunction(add_polymec_executable)

