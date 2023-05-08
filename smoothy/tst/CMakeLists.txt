add_executable(GMDS_SMOOTHY_TEST
        LaplacianSmootherTestSuite.h
        AngleBasedQuadSmootherTestSuite.h
        ElliptingSmoothingTestSuite.h
        main_test.cpp)
#==============================================================================
target_link_libraries(GMDS_SMOOTHY_TEST PUBLIC
        ${GMDS_LIB} ${LIB_GMDS_QUALITY}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_SMOOTHY_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
