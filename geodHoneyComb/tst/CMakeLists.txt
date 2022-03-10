add_executable(GMDS_GEOD_HONEY_COMB_TEST
        GeodHoneyCombTestSuite.h
        main_test.cpp )

target_link_libraries(GMDS_GEOD_HONEY_COMB_TEST PUBLIC
        ${GMDS_LIB}
        gtest)

target_include_directories(GMDS_GEOD_HONEY_COMB_TEST PRIVATE .)

set_target_properties(GMDS_GEOD_HONEY_COMB_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_GEOD_HONEY_COMB_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
