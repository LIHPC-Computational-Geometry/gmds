add_executable(GMDS_POLYBLOCK_TEST
        PolyBlockTestSuite.h
        PolyBlock2DTestSuite.h
        main_test.cpp)
#==============================================================================
target_link_libraries(GMDS_POLYBLOCK_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IG_ALGO}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_POLYBLOCK_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
