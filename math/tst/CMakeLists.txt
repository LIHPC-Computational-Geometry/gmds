add_executable(GMDS_MATH_TEST
        ChartTestSuite.h
        CrossTestSuite.h
        Cross2DTestSuite.h
        DiscretizationScheme1DTestSuite.h
        MathTestSuite.h
        PointTestSuite.h
        QuaternionTestSuite.h
        OrientationTestSuite.h
        TransfiniteInterpolationTestSuite.h
        main_test.cpp
        BezierSurfaceTestSuite.h
        BezierCurveTestSuite.h)
#==============================================================================
target_link_libraries(GMDS_MATH_TEST PUBLIC
        ${GMDS_LIB}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_MATH_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
