add_executable(GMDS_IGALGO_TEST
        BoundaryOpTestSuite.h
        BoundaryExtractorTestSuite.h
        GridBuilderTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_IGALGO_TEST PUBLIC
        ${GMDS_LIB} ${LIB_GMDS_IO}
        gtest)

target_include_directories(GMDS_IGALGO_TEST PRIVATE .)

set_target_properties(GMDS_IGALGO_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_IGALGO_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
