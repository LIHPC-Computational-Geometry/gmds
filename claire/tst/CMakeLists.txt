add_executable(GMDS_CLAIRE_TEST
        ClaireTestSuite.h
        LevelSetTestSuite.h
        DistanceMapTestSuite.h
        GradientComputation2DTestSuite.h
        GradientComputation3DTestSuite.h
        PointFollowingVectorField2DTestSuite.h
        main_test.cpp ../inc/gmds/claire/LevelSetFromIntToOut.h)

target_link_libraries(GMDS_CLAIRE_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IG_ALGO}
        gtest)

target_include_directories(GMDS_CLAIRE_TEST PRIVATE .)

set_target_properties(GMDS_CLAIRE_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_CLAIRE_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
