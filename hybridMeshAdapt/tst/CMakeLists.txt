add_executable(GMDS_SIMPLEX_INDEXING
        main_indexing_SURFACE_CURVE_CORNER.cpp)

add_executable(GMDS_SIMPLEX_TEST
        SimplexMeshTestAdd_NTt.h
        SimplexMeshTestSuite.h
        SimplexMeshTestReorient.h
        ModeleCAD1.h
        ModeleCAD0.h
        ModeleCAD2.h
        ModeleCAD3.h
        ModeleCAD4.h
        ModeleCAD5.h
        main_test.cpp)



target_link_libraries(GMDS_SIMPLEX_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IO}
        gtest)

target_link_libraries(GMDS_SIMPLEX_INDEXING PUBLIC
      ${GMDS_LIB}
      ${LIB_GMDS_IO}
      ${LIB_GMDS_FRAME_3D}
      gtest)

target_include_directories(GMDS_SIMPLEX_TEST PRIVATE .)

set_target_properties(GMDS_SIMPLEX_INDEXING PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)

set_target_properties(GMDS_SIMPLEX_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_SIMPLEX_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
