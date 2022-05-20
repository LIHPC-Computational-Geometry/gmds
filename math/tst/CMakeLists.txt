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
        main_test.cpp LineTestSuite.h)

target_link_libraries(GMDS_MATH_TEST PUBLIC
        ${GMDS_LIB}
        gtest)

target_include_directories(GMDS_MATH_TEST PRIVATE .)

set_target_properties(GMDS_MATH_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_MATH_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
