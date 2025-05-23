add_executable(GMDS_CLAIRE_TEST
        AeroTestSuite.h
        LevelSetTestSuite.h
        DistanceMapTestSuite.h
        GradientComputationTestSuite.h
        PointFollowingVectorFieldTestSuite.h
        AeroPipelineTestSuite.h
        AeroMeshQualityTestSuite.h
        FastLocalizeTestSuite.h
        main_test.cpp ../inc/gmds/aero/LevelSetCombined.h
        UtilsTestSuite.h
        AeroBoundariesTestSuite.h
        SU2WriterTestSuite.h
        RefinementBetaTestSuite.h
        SmoothLineSweepingTestSuite.h
        DiffusionEquation2DTestSuite.h
        NodeNeighbourhoodOnFront_3DTestSuite.h
        MeshAlignmentTestSuite.h
        Front_3DTestSuite.h
        FrontTestSuite.h
        MFEMMeshWriterTestSuite.h
        LayerStructureManager_3DTestSuite.h
        Blocking3DTestSuite.h
        TransfiniteInterpolation_3DTestSuite.h
        ComputeBezierDegree_3DTestSuite.h
        ComputeBezierCurveCtrlPtstoInterpolateCurveTestSuite.h
        RefinementBetaBlock3DTestSuite.h
        ComputeBezierCtrlPtstoInterpolateSurfaceTestSuite.h)
#==============================================================================
target_link_libraries(GMDS_CLAIRE_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IG_ALGO}
        ${LIB_GMDS_BLOCKING}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_CLAIRE_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
