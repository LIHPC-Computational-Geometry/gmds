add_executable(GMDS_UTILS_TEST
        UtilsTestSuite.h
        ParamTestSuite.h
        ArrayTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_UTILS_TEST PUBLIC
        ${GMDS_LIB}
        gtest)

target_include_directories(GMDS_UTILS_TEST PRIVATE .)

set_target_properties(GMDS_UTILS_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_UTILS_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
