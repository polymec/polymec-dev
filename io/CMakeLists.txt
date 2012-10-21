set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DARBI_C_COMPILER=${CMAKE_C_COMPILER} -DARBI_INCLUDE_DIR=${PROJECT_SOURCE_DIR}")
add_library(arbi_io silo_io.c vtk_plot_io.c 
            generate_face_node_conn.c generate_cell_node_conn.c)
set(ARBI_LIBS arbi_io;${ARBI_LIBS} PARENT_SCOPE)
set(ARBI_LIBS arbi_io;${ARBI_LIBS})
add_subdirectory(tests)

