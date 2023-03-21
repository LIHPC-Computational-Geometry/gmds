add_executable(GMDS_SMOOTHY_TEST
        LaplacianSmootherTestSuite.h
        AngleBasedQuadSmootherTestSuite.h
        ElliptingSmoothingTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_SMOOTHY_TEST PUBLIC
        ${GMDS_LIB} ${LIB_GMDS_QUALITY}
        gtest)

target_include_directories(GMDS_SMOOTHY_TEST PRIVATE .)

set_target_properties(GMDS_SMOOTHY_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_SMOOTHY_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
