set(TEST_SUITES
        ReaderTestSuite.h
        WriterTestSuite.h
        )

if(WITH_MESHB)
    list(APPEND TEST_SUITES MeshBWriterTestSuite.h)
endif()

if(WITH_LIMA)
    list(APPEND TEST_SUITES LimaWriteAndReadTestSuite.h)
endif()
add_executable(GMDS_IO_TEST ${TEST_SUITES}
        main_test.cpp)

target_link_libraries(GMDS_IO_TEST PUBLIC
        ${GMDS_LIB}
        gtest)

target_include_directories(GMDS_IO_TEST PRIVATE .)

set_target_properties(GMDS_IO_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_IO_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
