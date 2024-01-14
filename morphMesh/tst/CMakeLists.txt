add_executable(GMDS_MORPHMESH_TEST
        FastLocalizeTestSuite.h
        main_test.cpp
        )
#==============================================================================
target_link_libraries(GMDS_MORPHMESH_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IO}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_MORPHMESH_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
