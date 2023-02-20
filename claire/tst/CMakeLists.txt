add_executable(GMDS_CLAIRE_TEST
        ClaireTestSuite.h
        LevelSetTestSuite.h
        DistanceMapTestSuite.h
        GradientComputationTestSuite.h
        PointFollowingVectorFieldTestSuite.h
        AeroPipelineTestSuite.h
        AeroMeshQualityTestSuite.h
        FastLocalizeTestSuite.h
        main_test.cpp ../inc/gmds/claire/LevelSetCombined.h
        UtilsTestSuite.h
        AeroBoundariesTestSuite.h
        SU2WriterTestSuite.h
        RefinementBetaTestSuite.h
        SmoothLineSweepingTestSuite.h
        DiffusionEquation2DTestSuite.h)

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
