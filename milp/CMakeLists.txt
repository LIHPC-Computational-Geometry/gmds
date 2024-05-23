
#==============================================================================
# LIBRARY DEFINITION (SOURCE FILES)
#==============================================================================
# Explicitly used the name given in this preamble
set(GMDS_LIB ${LIB_GMDS_MILP})
set(GMDS_LIB_PREFIX gmds/milp)

set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/milp/milp.h
)
set(GMDS_SRC
        src/milp.cpp
        src/main_exemple_ortools.cpp
        src/InitializeLog.h
        src/InitializeLog.cpp
        #src/convertion.cpp
        #src/convertion.h
)
#==============================================================================
add_library(${GMDS_LIB} ${GMDS_INC} ${GMDS_SRC})
#==============================================================================
include(GenerateExportHeader)
generate_export_header(${GMDS_LIB}
        EXPORT_FILE_NAME ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        EXPORT_MACRO_NAME ${GMDS_LIB}_API)

find_package(ortools)
#==============================================================================
# TARGET DEFINITION
#==============================================================================
include(GNUInstallDirs)
#LIBRARY TO INSTALL
target_link_libraries(${GMDS_LIB} PUBLIC
        ${LIB_GMDS_CAD}
        ${LIB_GMDS_IG}
        ortools::ortools
)


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

add_executable(milp_exec src/main.cpp
        #src/convertion.cpp
        #src/convertion.h
)
target_link_libraries(milp_exec PRIVATE ${GMDS_LIB})
target_compile_features(milp_exec PUBLIC cxx_std_14)
install(TARGETS milp_exec)

add_executable(ortools_exec src/main_ortools.cpp
#        src/convertion.cpp
#        src/convertion.h
)
target_link_libraries(ortools_exec PRIVATE ${GMDS_LIB} ortools::ortools)
target_compile_features(ortools_exec PUBLIC cxx_std_14)
install(TARGETS ortools_exec)


add_executable(pixelvf2d_exec src/main_pixelVF2d.cpp inc/gmds/milp/milp_InterfaceNodesPosSmoothVF.h src/milp_InterfaceNodesPosSmoothVF.cpp)
target_link_libraries(pixelvf2d_exec PRIVATE ${GMDS_LIB} elg3d  )
target_include_directories(pixelvf2d_exec PRIVATE ../Elg3D/inc)
target_compile_features(pixelvf2d_exec PUBLIC cxx_std_14)
install(TARGETS pixelvf2d_exec)

add_executable(pixelvf2d_exec_ortools src/main_pixelVF2d_ortools.cpp inc/gmds/milp/milp_InterfaceNodesPosSmoothVF.h src/milp_InterfaceNodesPosSmoothVF.cpp inc/gmds/milp/ortools_InterfaceNodesPosSmoothVF.h src/ortools_InterfaceNodesPosSmoothVF.cpp)
target_link_libraries(pixelvf2d_exec_ortools PRIVATE ${GMDS_LIB} elg3d  )
target_include_directories(pixelvf2d_exec_ortools PRIVATE ../Elg3D/inc)
target_compile_features(pixelvf2d_exec_ortools PUBLIC cxx_std_14)
install(TARGETS pixelvf2d_exec_ortools)