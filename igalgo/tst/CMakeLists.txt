add_executable(GMDS_IGALGO_TEST
        BoundaryOpTestSuite.h
        BoundaryExtractorTestSuite.h
        GridBuilderTestSuite.h
        main_test.cpp)
#==============================================================================
target_link_libraries(GMDS_IGALGO_TEST PUBLIC
        ${GMDS_LIB} ${LIB_GMDS_IO}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_IGALGO_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
