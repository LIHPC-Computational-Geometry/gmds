add_executable(GMDS_SIMPLEX_INSERTION_TEST
        main_INSERTION_TEST.cpp)

add_executable(GMDS_SIMPLEX_HEX_GENERATION
        main_HEX_GENERATION.cpp)

add_executable(GMDS_SIMPLEX_POINT_GENERATION
        main_POINT_GENERATION.cpp)

add_executable(GMDS_SIMPLEX_INDEXING
        main_INDEXING_SURFACE_CURVE_CORNER.cpp)

add_executable(GMDS_SIMPLEX_INDEXING_CORRECTION
        main_INDEXING_CORRECTION.cpp)

add_executable(GMDS_SIMPLEX_METRIC_BASED_ADAPTATION
        main_METRIC_BASED_ADAPTATION.cpp)

add_executable(GMDS_SIMPLEX_METRIC_FF_POINTGEN
        main_METRIC_FF_PointGeneration.cpp)

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
        GTest::gtest)

target_link_libraries(GMDS_SIMPLEX_METRIC_FF_POINTGEN PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IO}
        GTest::gtest)

target_link_libraries(GMDS_SIMPLEX_METRIC_BASED_ADAPTATION PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IO}
        GTest::gtest)

target_link_libraries(GMDS_SIMPLEX_INDEXING PUBLIC
      ${GMDS_LIB}
      ${LIB_GMDS_IO}
      ${LIB_GMDS_FRAME_3D}
      )

target_link_libraries(GMDS_SIMPLEX_INDEXING_CORRECTION PUBLIC
${GMDS_LIB}
${LIB_GMDS_IO}
${LIB_GMDS_FRAME_3D}
      )

target_link_libraries(GMDS_SIMPLEX_POINT_GENERATION PUBLIC
      ${GMDS_LIB}
      ${LIB_GMDS_IO}
      ${LIB_GMDS_FRAME_3D}
      )

target_link_libraries(GMDS_SIMPLEX_HEX_GENERATION PUBLIC
      ${GMDS_LIB}
      ${LIB_GMDS_IO}
      GTest::gtest)

target_link_libraries(GMDS_SIMPLEX_INSERTION_TEST PUBLIC
      ${GMDS_LIB}
      ${LIB_GMDS_IO}
      GTest::gtest)


set_target_properties(GMDS_SIMPLEX_INDEXING PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)

#==============================================================================
gtest_discover_tests(GMDS_SIMPLEX_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
