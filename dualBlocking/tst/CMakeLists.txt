add_executable(GMDS_DUBLO_TEST
        DuBloTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_DUBLO_TEST PUBLIC
        ${GMDS_LIB}
        gtest)

target_include_directories(GMDS_DUBLO_TEST PRIVATE .)

set_target_properties(GMDS_DUBLO_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_DUBLO_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
