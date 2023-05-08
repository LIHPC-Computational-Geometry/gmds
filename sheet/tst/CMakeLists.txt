add_executable(GMDS_SHEET_TEST
        SheetPillowTestSuite.h
        SheetCollapseTestSuite.h
        SheetSelectorTestSuite.h
        main_test.cpp)
#==============================================================================
target_link_libraries(GMDS_SHEET_TEST PUBLIC
        ${GMDS_LIB} ${LIB_GMDS_IO} ${LIB_GMDS_SMOOTHY}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_SHEET_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
