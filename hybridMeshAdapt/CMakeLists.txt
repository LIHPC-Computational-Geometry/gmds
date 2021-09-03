#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_HYBRIDMESHADAPT})
set(GMDS_LIB_PREFIX gmds/hybridMeshAdapt/)

set(GMDS_INC
        inc/gmds/hybridMeshAdapt/SimplicesNode.h
        inc/gmds/hybridMeshAdapt/SimplexMesh.h
        inc/gmds/hybridMeshAdapt/SimplicesCell.h
        inc/gmds/hybridMeshAdapt/AccessBehavior.h
        inc/gmds/hybridMeshAdapt/CommonInfo.h
        inc/gmds/hybridMeshAdapt/SimplicesTriangle.h
        inc/gmds/hybridMeshAdapt/CavityOperator.h
        inc/gmds/hybridMeshAdapt/ISimplexMeshIOService.h
        inc/gmds/hybridMeshAdapt/PointInsertion.h
        inc/gmds/hybridMeshAdapt/PointSmoothing.h
        inc/gmds/hybridMeshAdapt/EdgeCollapse.h
        inc/gmds/hybridMeshAdapt/EdgeInsertion.h
        inc/gmds/hybridMeshAdapt/StructuredGrid.h
        inc/gmds/hybridMeshAdapt/DelaunayPointInsertion.h
        inc/gmds/hybridMeshAdapt/Metric.h
        inc/gmds/hybridMeshAdapt/Octree.h
        )

set(GMDS_SRC
        src/SimplicesNode.cpp
        src/SimplexMesh.cpp
        src/SimplicesNode.cpp
        src/SimplicesCell.cpp
        src/AccessBehavior.cpp
        #src/CommonInfo.cpp
        src/SimplicesTriangle.cpp
        src/ISimplexMeshIOService.cpp
        src/CavityOperator.cpp
        src/ICriterion.cpp
        src/PointInsertion.cpp
        src/PointSmoothing.cpp
        src/EdgeCollapse.cpp
        src/EdgeInsertion.cpp
        src/StructuredGrid.cpp
        src/DelaunayPointInsertion.cpp
        src/Metric.cpp
        src/Octree.cpp
        )
#==============================================================================
add_library(${GMDS_LIB} ${GMDS_INC} ${GMDS_SRC})
#==============================================================================
# TARGET DEFINITION
#==============================================================================
include(GNUInstallDirs)

#LIBRARY TO INSTALL
target_link_libraries(${GMDS_LIB} PUBLIC
        ${LIB_GMDS_IG}
        ${LIB_GMDS_IO})

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
add_executable(HybridMeshAdapter src/main.cpp)
target_link_libraries(HybridMeshAdapter PRIVATE ${GMDS_LIB})
target_compile_features(HybridMeshAdapter PUBLIC cxx_std_14)
