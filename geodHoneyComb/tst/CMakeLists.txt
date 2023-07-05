add_executable(GMDS_GEOD_HONEY_COMB_TEST
        GeodHoneyCombTestSuite.h
        RegularIsocahedronTestSuite.h
        main_test.cpp )
#==============================================================================
target_link_libraries(GMDS_GEOD_HONEY_COMB_TEST PUBLIC
        ${GMDS_LIB}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_GEOD_HONEY_COMB_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
