#==============================================================================
add_executable(GMDS_MCTSBLOCK_TEST
        BlockingTestSuite.h
        BlockingClassifierTestSuite.h
        MCTSTestSuite.h
        GraphTestSuite.h
        DataMCTSTestSuite.h
        main_test.cpp
        )
#==============================================================================
target_link_libraries(GMDS_MCTSBLOCK_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IG}
        CGAL::CGAL_Core
        GTest::gtest)

#==============================================================================
gtest_discover_tests(GMDS_MCTSBLOCK_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
