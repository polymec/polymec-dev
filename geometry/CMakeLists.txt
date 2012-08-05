add_library(arbi_geometry planar_voronoi.c voronoi.c
            create_box_mesh.c create_cyl_mesh.c create_sphere_mesh.c)
set(ARBI_LIBS arbi_geometry;${ARBI_LIBS} PARENT_SCOPE)
set(ARBI_LIBS arbi_geometry;${ARBI_LIBS})

