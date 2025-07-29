#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_AERO})
set(GMDS_LIB_PREFIX gmds/aero)

set(GMDS_INC
		${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/aero/AbstractAeroBoundaries.h
        inc/gmds/aero/AbstractAeroPipeline.h
        inc/gmds/aero/AbstractLevelSet.h
        inc/gmds/aero/AbstractPatternEdge.h
        inc/gmds/aero/AbstractPatternNode.h
        inc/gmds/aero/AbstractSmoothLineSweeping_2D.h
        inc/gmds/aero/AdvectedPointRK4_2D.h
        inc/gmds/aero/AdvectedPointRK4_3D.h
        inc/gmds/aero/AeroBoundaries_2D.h
        inc/gmds/aero/AeroBoundaries_3D.h
        inc/gmds/aero/AeroEllipticSmoothing_2D.h
        inc/gmds/aero/AeroException.h
        inc/gmds/aero/AeroExtrusion_2D.h
        inc/gmds/aero/AeroExtrusion_3D.h
        inc/gmds/aero/AeroMeshQuality.h
        inc/gmds/aero/AeroPipeline_2D.h
        inc/gmds/aero/AeroPipeline_3D.h
        inc/gmds/aero/Blocking3D.h
        inc/gmds/aero/ComputeBezierCtrlPtstoInterpolateSurface.h
        inc/gmds/aero/ComputeBezierCurveCtrlPtstoInterpolateCurve.h
        inc/gmds/aero/ComputeBezierDegree_3D.h
        inc/gmds/aero/ControlPointsSmoothing_2D.h
        inc/gmds/aero/DiffusionEquation2D.h
        inc/gmds/aero/DistanceMap.h
        inc/gmds/aero/FastLocalize.h
        inc/gmds/aero/Front.h
        inc/gmds/aero/Front_3D.h
        inc/gmds/aero/FrontEdgesNodesClassification_3D.h
        inc/gmds/aero/GradientComputation_2D.h
        inc/gmds/aero/GradientComputation_3D.h
        inc/gmds/aero/Grid_Smooth2D.h
        inc/gmds/aero/IntervalAssignment_2D.h
        inc/gmds/aero/IntervalAssignment_3D.h
        inc/gmds/aero/LayerStructureManager_3D.h
        inc/gmds/aero/LeastSquaresGradientComputation.h
        inc/gmds/aero/LevelSetCombined.h
        inc/gmds/aero/LevelSetEloi.h
        inc/gmds/aero/LevelSetExtended.h
        inc/gmds/aero/MeshAlignment_2D.h
        inc/gmds/aero/MFEMMeshWriter.h
        inc/gmds/aero/NodeNeighbourhoodOnFront_3D.h
        inc/gmds/aero/Params.h
        inc/gmds/aero/PatternEdgeCorner.h
        inc/gmds/aero/PatternEdgeEnd.h
        inc/gmds/aero/PatternEdgeReversal.h
        inc/gmds/aero/PatternFace.h
        inc/gmds/aero/PatternNode1Corner2End.h
        inc/gmds/aero/PatternNode2Corner1End.h
        inc/gmds/aero/PatternNode2Corner1Reversal.h
        inc/gmds/aero/PatternNode2Corner2End.h
        inc/gmds/aero/PatternNode2End1Reversal.h
        inc/gmds/aero/PatternNode3Corner.h
        inc/gmds/aero/PatternNode3Corner3End.h
        inc/gmds/aero/PatternNode3End.h
        inc/gmds/aero/RefinementBeta.h
        inc/gmds/aero/RefinementBetaBlock3D.h
        inc/gmds/aero/RefinementBetaBlocking.h
        inc/gmds/aero/RefinementBetaBlocking3D.h
        inc/gmds/aero/Smooth2D.h
        inc/gmds/aero/SmoothingPaving_2D.h
        inc/gmds/aero/SmoothLineSweepingOrtho.h
        inc/gmds/aero/SmoothLineSweepingYao.h
        inc/gmds/aero/SU2Writer.h
        inc/gmds/aero/SU2Writer_3D.h
        inc/gmds/aero/TransfiniteInterpolation_3D.h
        inc/gmds/aero/Utils.h
		inc/gmds/aero/BLKWriter.h
		)
if(WITH_CGNS)
	set(GMDS_SRC ${GMDS_SRC} inc/gmds/aero/CGNSWriter3D.h
	)
endif()

set(GMDS_SRC
        src/AbstractAeroBoundaries.cpp
        src/AbstractAeroPipeline.cpp
        src/AbstractLevelSet.cpp
        src/AbstractPatternEdge.cpp
        src/AbstractPatternNode.cpp
        src/AbstractSmoothLineSweeping_2D.cpp
        src/AdvectedPointRK4_2D.cpp
        src/AdvectedPointRK4_3D.cpp
        src/AeroBoundaries_2D.cpp
        src/AeroBoundaries_3D.cpp
        src/AeroEllipticSmoothing_2D.cpp
        src/AeroException.cpp
        src/AeroExtrusion_2D.cpp
        src/AeroExtrusion_3D.cpp
        src/AeroMeshQuality.cpp
        src/AeroPipeline_2D.cpp
        src/AeroPipeline_3D.cpp
        src/Blocking3D.cpp
        src/ComputeBezierCtrlPtstoInterpolateSurface.cpp
        src/ComputeBezierCurveCtrlPtstoInterpolateCurve.cpp
        src/ComputeBezierDegree_3D.cpp
        src/ControlPointsSmoothing_2D.cpp
        src/DiffusionEquation2D.cpp
        src/DistanceMap.cpp
        src/FastLocalize.cpp
        src/Front.cpp
        src/Front_3D.cpp
        src/FrontEdgesNodesClassification_3D.cpp
        src/GradientComputation_2D.cpp
        src/GradientComputation_3D.cpp
        src/Grid_Smooth2D.cpp
        src/IntervalAssignment_2D.cpp
        src/IntervalAssignment_3D.cpp
        src/LeastSquaresGradientComputation.cpp
        src/LayerStructureManager_3D.cpp
        src/LevelSetCombined.cpp
        src/LevelSetEloi.cpp
        src/LevelSetExtended.cpp
        src/MeshAlignment_2D.cpp
        src/MFEMMeshWriter.cpp
        src/NodeNeighbourhoodOnFront_3D.cpp
        src/PatternEdgeCorner.cpp
        src/PatternEdgeEnd.cpp
        src/PatternEdgeReversal.cpp
        src/PatternFace.cpp
        src/PatternNode1Corner2End.cpp
        src/PatternNode2Corner1End.cpp
        src/PatternNode2Corner1Reversal.cpp
        src/PatternNode2Corner2End.cpp
        src/PatternNode2End1Reversal.cpp
        src/PatternNode3Corner.cpp
        src/PatternNode3Corner3End.cpp
        src/PatternNode3End.cpp
        src/RefinementBeta.cpp
        src/RefinementBetaBlock3D.cpp
        src/RefinementBetaBlocking.cpp
        src/RefinementBetaBlocking3D.cpp
        src/Smooth2D.cpp
        src/SmoothingPaving_2D.cpp
        src/SmoothLineSweepingOrtho.cpp
        src/SmoothLineSweepingYao.cpp
        src/SU2Writer.cpp
        src/SU2Writer_3D.cpp
        src/TransfiniteInterpolation_3D.cpp
        src/Utils.cpp
        src/main.cpp
		src/BLKWriter.cpp
)
if(WITH_CGNS)
	set(GMDS_SRC ${GMDS_SRC} src/CGNSWriter3D.cpp
	)
endif()
#==============================================================================
add_library(${GMDS_LIB} ${GMDS_INC} ${GMDS_SRC})
#==============================================================================
include(GenerateExportHeader)
generate_export_header(${GMDS_LIB}
		EXPORT_FILE_NAME ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
		EXPORT_MACRO_NAME ${GMDS_LIB}_API)
#==============================================================================
# TARGET DEFINITION
#==============================================================================
include(GNUInstallDirs)

#LIBRARY TO INSTALL
target_link_libraries(${GMDS_LIB} PUBLIC
		${LIB_GMDS_IG}
		${LIB_GMDS_IO}
		${LIB_GMDS_IG_ALGO}
		${LIB_GMDS_CAD}
		${LIB_GMDS_QUALITY}
		${LIB_GMDS_CADFAC}
		${LIB_GMDS_SMOOTHY}
		${LIB_GMDS_BLOCKING}
		GMDSSmoothy
		ANN
		)

target_compile_features(${GMDS_LIB} PUBLIC cxx_std_14)

# INCLUDE TO INSTALL
target_include_directories(${GMDS_LIB} PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
		$<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
		)
set_target_properties(${GMDS_LIB} PROPERTIES PUBLIC_HEADER "${GMDS_INC}")

install(TARGETS ${GMDS_LIB}
		EXPORT GMDS_SUITE
		PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${GMDS_LIB_PREFIX}
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
		ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
		RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

#==============================================================================
if(WITH_TEST)
	add_subdirectory(tst)
endif(WITH_TEST)

#==============================================================================
# EXECUTABLE
#==============================================================================
add_executable(Aero src/main.cpp)
target_link_libraries(Aero PRIVATE ${GMDS_LIB})
target_compile_features(Aero PUBLIC cxx_std_14)
install(TARGETS Aero)

add_executable(aero_dev src/main_dev.cpp)
target_link_libraries(aero_dev PRIVATE ${GMDS_LIB})
target_compile_features(aero_dev PUBLIC cxx_std_14)
