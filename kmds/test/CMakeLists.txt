set(UT_KMDS_SRC_FILES
		ConnectivityTest.h
		ContainersTest.h
		EContainerTest.h
		FaceTest.h
		FakeTypesTest.h
		FContainersTest.h
		GrowingViewTest.h
		InitToolsTest.h
		KokkosUnorderedMapTest.h
		MeshTest.h
		NContainerTest.h
		RContainersTest.h
        RegionTest.h
		VariableTest.h
		WriterTest.h
        GraphTest.h)

add_executable(UT_kmds main.cpp ${UT_KMDS_SRC_FILES})

target_link_libraries(UT_kmds GTest::gtest kmds ${LIB_GMDS_IG})

#add_custom_command(
#   TARGET UT_kmds
#   POST_BUILD
#   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:UT_kmds> ${CMAKE_BINARY_DIR}/Testing/UT_kmds
#   COMMENT "Copy kmds unit tests"
#   VERBATIM
#   )
#
#
#add_test (NAME    UT_kmds
#	  COMMAND  $<TARGET_FILE:UT_kmds>
#	  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})


#add_subdirectory(Samples)
gtest_discover_tests(UT_kmds
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
