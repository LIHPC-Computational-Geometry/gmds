add_executable(GMDS_BLOCKING_TEST
        BlockingTestSuite.h
        CGALTestSuite.h
        InputMarkedDartsTestSuite.h
        SheetInsertTestSuite.h
        WriterDartsVTKTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_BLOCKING_TEST PUBLIC
        ${GMDS_LIB}
        CGAL::CGAL_Core
        gtest)

if (WITH_CGAL_QT5)
    target_link_libraries(${GMDS_LIB} PUBLIC
            CGAL::CGAL_Basic_viewer
            )
endif()

target_include_directories(GMDS_BLOCKING_TEST PRIVATE .)

set_target_properties(GMDS_BLOCKING_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_BLOCKING_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
