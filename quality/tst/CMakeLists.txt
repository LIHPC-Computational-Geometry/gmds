add_executable(GMDS_QUALITY_TEST
        QuadEvalTestSuite.h
        HexEvalTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_QUALITY_TEST PUBLIC
        ${GMDS_LIB}
        gtest)

target_include_directories(GMDS_QUALITY_TEST PRIVATE .)

set_target_properties(GMDS_QUALITY_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_QUALITY_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
