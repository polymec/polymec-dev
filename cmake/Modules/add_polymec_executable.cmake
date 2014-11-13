# This function adds an executable program that links against polymec.
function(add_polymec_executable exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} ${POLYMEC_LIBRARIES})
endfunction(add_polymec_executable)

