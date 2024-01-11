#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_RLBLOCKING})
set(GMDS_LIB_PREFIX gmds/rlBlocking)
set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/rlBlocking/BlockingQuality.h
        inc/gmds/rlBlocking/LinkerBlockingGeom.h
        inc/gmds/rlBlocking/ValidBlocking.h
        inc/gmds/rlBlocking/MCTSAlgorithm.h
        inc/gmds/rlBlocking/MCTSState.h
        inc/gmds/rlBlocking/MCTSTree.h
        inc/gmds/rlBlocking/MCTSMove.h
        inc/gmds/rlBlocking/MCTSMovePolycube.h
        inc/gmds/rlBlocking/MCTSStatePolycube.h
        inc/gmds/rlBlocking/MCTSAgent.h
        )
set(GMDS_SRC
        src/BlockingQuality.cpp
        src/LinkerBlockingGeom.cpp
        src/ValidBlocking.cpp
        src/MCTSTree.cpp
        src/MCTSAlgorithm.cpp
        src/MCTSMovePolycube.cpp
        src/MCTSStatePolycube.cpp
        src/MCTSAgent.cpp
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
        ${LIB_GMDS_IG}
        ${LIB_GMDS_IO}
        ${LIB_GMDS_IG_ALGO}
        ${LIB_GMDS_QUALITY}
        ${LIB_GMDS_BLOCK_MESHER}
        ${LIB_GMDS_BLOCKING}
        )

target_compile_features(${GMDS_LIB} PUBLIC cxx_std_14)

# INCLUDE TO INSTALL
target_include_directories(${GMDS_LIB} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
        )
set_target_properties(${GMDS_LIB} PROPERTIES PUBLIC_HEADER "${GMDS_INC}")


#==============================================================================
if(WITH_TEST)
    add_subdirectory(tst)
endif(WITH_TEST)

#==============================================================================
# EXECUTABLE
#==============================================================================
add_executable(RlBlocking src/main_rlBlocking.cpp)
target_link_libraries(RlBlocking PRIVATE ${GMDS_LIB} )
target_compile_features(RlBlocking PUBLIC cxx_std_14)
