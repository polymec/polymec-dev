add_library(arbi_geometry voronoi.c planar_voronoi.c create_cvt.c
            plane.c cylinder.c sphere.c intersection.c cubic_lattice.c
            create_cubic_lattice_mesh.c)
            #create_box_mesh.c create_cyl_mesh.c create_sphere_mesh.c)
set(ARBI_LIBS arbi_geometry;${ARBI_LIBS} PARENT_SCOPE)
set(ARBI_LIBS arbi_geometry;${ARBI_LIBS})

add_subdirectory(tests)
