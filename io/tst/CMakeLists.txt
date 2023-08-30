set(TEST_SUITES
        ReaderTestSuite.h
        WriterTestSuite.h
        MdlReaderTestSuite.h
        )

if(WITH_MESHB)
    list(APPEND TEST_SUITES MeshBWriterTestSuite.h)
endif()

if(WITH_LIMA)
    list(APPEND TEST_SUITES LimaWriteAndReadTestSuite.h)
endif()
add_executable(GMDS_IO_TEST ${TEST_SUITES}
        main_test.cpp)
#==============================================================================
target_link_libraries(GMDS_IO_TEST PUBLIC
        ${GMDS_LIB}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_IO_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
