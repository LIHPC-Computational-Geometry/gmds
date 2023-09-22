add_executable(GMDS_BLOCKING_TEST
        BlockingTestSuite.h
        CurvedBlockingTestSuite.h
        CurvedBlockingClassificationTestSuite.h
        CGALTestSuite.h
        InputMarkedDartsTestSuite.h
        SheetInsertTestSuite.h
        WriterDartsVTKTestSuite.h
        main_test.cpp
        ExecutionActionsTestSuite.h
        ClassificationTestSuite.h
)
#==============================================================================

if(WITH_CGNS)
    set(GMDS_BLOCKING_TEST_FILES ${GMDS_BLOCKING_TEST_FILES} CGNSWriterTestSuite.h)
    target_compile_definitions(GMDS_BLOCKING_TEST PRIVATE USE_CGNS)
endif()
#==============================================================================
target_link_libraries(GMDS_BLOCKING_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IG}
        CGAL::CGAL_Core
        GTest::gtest)

if(WITH_CGNS)
target_link_libraries(GMDS_BLOCKING_TEST PUBLIC cgns::cgns)
endif()

if (WITH_CGAL_QT5)
    target_link_libraries(GMDS_BLOCKING_TEST PUBLIC
            CGAL::CGAL_Basic_viewer
            )
endif()
#==============================================================================
gtest_discover_tests(GMDS_BLOCKING_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
