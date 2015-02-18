# This macro makes adjustments to compiler flags for certain machines.
macro(adjust_compilers)

  include(CMakeForceCompiler)
  if (CMAKE_ADJUSTED_C_COMPILER)
    CMAKE_FORCE_C_COMPILER(${CMAKE_ADJUSTED_C_COMPILER} ${CMAKE_ADJUSTED_C_COMPILER_ID})
  endif()

  if (CMAKE_ADJUSTED_CXX_COMPILER)
    CMAKE_FORCE_CXX_COMPILER(${CMAKE_ADJUSTED_CXX_COMPILER} ${CMAKE_ADJUSTED_CXX_COMPILER_ID})
  endif()

endmacro()
