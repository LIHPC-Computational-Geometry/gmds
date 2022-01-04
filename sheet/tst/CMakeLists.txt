add_executable(GMDS_SHEET_TEST
        SheetPillowTestSuite.h
        SheetCollapseTestSuite.h
        SheetSelectorTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_SHEET_TEST PUBLIC
        ${GMDS_LIB} ${LIB_GMDS_IO} ${LIB_GMDS_SMOOTHY}
        gtest)

target_include_directories(GMDS_SHEET_TEST PRIVATE .)

set_target_properties(GMDS_SHEET_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_SHEET_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
