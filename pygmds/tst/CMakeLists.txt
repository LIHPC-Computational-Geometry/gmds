add_executable(GMDS_POLYBLOCK_TEST
        PolyBlockTestSuite.h
        PolyBlock2DTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_POLYBLOCK_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IG_ALGO}
        gtest)

target_include_directories(GMDS_POLYBLOCK_TEST PRIVATE .)

set_target_properties(GMDS_POLYBLOCK_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_POLYBLOCK_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
