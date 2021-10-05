add_executable(GMDS_PADDING_TEST
        PaddingTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_PADDING_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IG_ALGO}
        gtest)

target_include_directories(GMDS_PADDING_TEST PRIVATE .)

set_target_properties(GMDS_PADDING_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_PADDING_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
