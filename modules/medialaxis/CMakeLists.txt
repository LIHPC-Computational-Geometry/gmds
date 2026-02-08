#==============================================================================
# LIBRARY DEFINITION (SOURCE FILES)
#==============================================================================
# Explicitly used the name given in this preamble
set(GMDS_LIB ${LIB_GMDS_MEDIAL_AXIS})
set(GMDS_LIB_PREFIX gmds/medialaxis)

set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/medialaxis/MedialAxis2D.h
        inc/gmds/medialaxis/MedialAxisMath.h
        inc/gmds/medialaxis/MedialAxis2DBuilder.h
        inc/gmds/medialaxis/MedialAxis3D.h
        inc/gmds/medialaxis/MedialAxis3DBuilder.h
        inc/gmds/medialaxis/CrossField.h
        inc/gmds/medialaxis/NonConformalHalfEdge.h
        inc/gmds/medialaxis/QuantizationSolver.h
        inc/gmds/medialaxis/QuantizationGraph.h
        inc/gmds/medialaxis/ConformalMeshBuilder.h
        inc/gmds/medialaxis/BlockStructureSimplifier.h
        inc/gmds/medialaxis/MinDelaunayCleaner.h
        inc/gmds/medialaxis/Conformalizer.h
        inc/gmds/medialaxis/MedaxBasedTMeshBuilder.h
        inc/gmds/medialaxis/TrianglesRemover.h
)
set(GMDS_SRC
        src/MedialAxis2D.cpp
        src/MedialAxisMath.cpp
        src/MedialAxis2DBuilder.cpp
        src/MedialAxis3D.cpp
        src/MedialAxis3DBuilder.cpp
        src/CrossField.cpp
        src/NonConformalHalfEdge.cpp
        src/QuantizationSolver.cpp 
        src/QuantizationGraph.cpp
        src/ConformalMeshBuilder.cpp
        src/BlockStructureSimplifier.cpp
        src/MinDelaunayCleaner.cpp
        src/Conformalizer.cpp
        src/MedaxBasedTMeshBuilder.cpp
        src/TrianglesRemover.cpp
)
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
        ${LIB_GMDS_CAD}
        ${LIB_GMDS_IG})

#==============================================================================
# NOTHING TO UPDATE BELOW
#==============================================================================

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
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gmds
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gmds)

#==============================================================================
if(WITH_TEST)
        add_subdirectory(tst)
endif(WITH_TEST)
#==============================================================================
# EXECUTABLE
#==============================================================================
add_executable(medialaxis_exec src/main_medax_build.cpp)
target_link_libraries(medialaxis_exec PRIVATE ${GMDS_LIB})
target_compile_features(medialaxis_exec PUBLIC cxx_std_14)
install(TARGETS medialaxis_exec)

add_executable(medialaxis_crossfield_exec src/main_crossfield.cpp)
target_link_libraries(medialaxis_crossfield_exec PRIVATE ${GMDS_LIB})
target_compile_features(medialaxis_crossfield_exec PUBLIC cxx_std_14)
install(TARGETS medialaxis_crossfield_exec)

add_executable(medialaxis_block_decomp_exec src/main_block_decomp.cpp)
target_link_libraries(medialaxis_block_decomp_exec PRIVATE ${GMDS_LIB})
target_compile_features(medialaxis_block_decomp_exec PUBLIC cxx_std_14)
install(TARGETS medialaxis_block_decomp_exec)

add_executable(quantization_exec src/main_quantization.cpp)
target_link_libraries(quantization_exec PRIVATE ${GMDS_LIB})
target_compile_features(quantization_exec PUBLIC cxx_std_14)
install(TARGETS quantization_exec)