# Executable.
add_executable(polymesher polymesher.c)
set_target_properties(polymesher PROPERTIES OUTPUT_NAME polymesher)
target_link_libraries(polymesher ${POLYMEC_LIBS})

