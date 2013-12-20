if (${HAVE_TETGEN})
  set(MESH_GENERATION_SOURCES create_voronoi_mesh.c #create_cvt.c create_cvt_with_lloyd_iteration.c 
      create_boundary_generators.c generate_random_points.c point_factory.c 
      mesh_factory.c cvt_gen_dist.c prob_cvt_gen_dist.c #force_cvt_gen_dist.c 
      witkin_heckbert_sampling.c mgw_sampling.c)
endif()

add_library(polymec_geometry mapping.c plane.c cylinder.c sphere.c intersection.c 
            union.c difference.c scaled.c rect_prism.c 
            giftwrap_hull.c polygon2.c polygon.c
            ${MESH_GENERATION_SOURCES}
            interpreter_register_geometry_functions.c 
            interpreter_register_spfuncs.c)

set(POLYMEC_LIBS polymec_geometry;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_geometry;${POLYMEC_LIBS})

add_subdirectory(tests)

install(TARGETS polymec_geometry DESTINATION lib)
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polymec/geometry)
