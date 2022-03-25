set(UT_ELG3D_SRC_FILES

		AssignCellsTest.h
		BadPillowingTest.h
		DefeaturingTest.h
		DummyTest.h
		EigenTest.h
		ExtractGeomTest.h
		FracPresTest.h
		FracPresEnforcementTest.h
		GCOTest.h
		InihTest.h
		InitDataTest.h
		InterfaceNodesPosSmoothVFTest.h
		MaterialAssignmentTest.h
		OptimizationSmoothTest.h
		Refinement.h
		SubsetProblemTest.h
        MaterialInterfacesTest.h ManifoldDetectionTest.h InterfaceNodesPosTest.h BoundingBoxGeomAssociationTest.h MaterialGradientComputationTest.h SmartLaplacianTest.h MoveToNewPosTest.h PillowingTest.h)

add_executable(UT_elg3d main.cpp ${UT_ELG3D_SRC_FILES})

target_link_libraries(UT_elg3d
		PUBLIC gmock elg3d kmds Kokkos::kokkoscore GCO INIH
		INTERFACE glpk)



#add_custom_command(
#   TARGET UT_elg3d
#   POST_BUILD
#   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:UT_elg3d> ${CMAKE_BINARY_DIR}/Testing/UT_elg3d
#   COMMENT "Copy elg3d unit tests"
#   VERBATIM
#   )


add_test (NAME    UT_elg3d
	  COMMAND  $<TARGET_FILE:UT_elg3d>
	  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})


#add_subdirectory(Samples)
