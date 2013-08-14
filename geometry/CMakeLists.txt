# To QHull or to Tetgen... that is the question.
# NOTE: Currently Tetgen produces an unbounded tessellation, whereas QHull 
# NOTE: produces only the bounded cells.
set(TETGEN 0)
if (TETGEN EQUAL 1)
  set(VORONOI_SOURCES voronoi_tessellator_tetgen.cpp;delaunay_triangulator_tetgen.cpp)
else()
  set(VORONOI_SOURCES voronoi_tessellator_qhull.c)
endif()

add_library(polymec_geometry mapping.c plane.c cylinder.c sphere.c intersection.c 
            union.c difference.c scaled.c rect_prism.c 
            giftwrap_hull.c polygon2.c polygon.c
            cubic_lattice.c create_cubic_lattice_mesh.c cvt_gen_dist.c
            prob_cvt_gen_dist.c #force_cvt_gen_dist.c 
            witkin_heckbert_sampling.c mgw_sampling.c
            create_voronoi_mesh.c ${VORONOI_SOURCES}
            create_cvt.c create_cvt_with_lloyd_iteration.c
            generate_random_points.c point_factory.c mesh_factory.c
            interpreter_register_geometry_functions.c 
            interpreter_register_spfuncs.c)

set(POLYMEC_LIBS polymec_geometry;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_geometry;${POLYMEC_LIBS})

add_subdirectory(tests)
