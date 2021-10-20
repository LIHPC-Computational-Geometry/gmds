#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_DUAL_BLOCKING})
set(GMDS_LIB_PREFIX gmds/dualBlocking)

set(GMDS_INC
        inc/gmds/dualBlocking/DualBlockingSession.h
        inc/gmds/dualBlocking/DualSurfaceCreator.h
        inc/gmds/dualBlocking/BoundarySurfaceCreator.h
        inc/gmds/dualBlocking/DualSheetCreator.h
        inc/gmds/dualBlocking/DualSheet.h
        inc/gmds/dualBlocking/EdgeDiscrAlgo.h
        )
set(GMDS_SRC
        src/DualSheetCreator.cpp
        src/DualSheet.cpp
        src/DualBlockingSession.cpp
        src/DualSurfaceCreator.cpp
        src/BoundarySurfaceCreator.cpp
        src/EdgeDiscrAlgo.cpp )
#==============================================================================
add_library(${GMDS_LIB} ${GMDS_INC} ${GMDS_SRC})
#==============================================================================
# TARGET DEFINITION
#==============================================================================
include(GNUInstallDirs)


if(NOT WIN32)
    if((NOT DEFINED GLPK_LIB) OR (NOT DEFINED GLPK_INC))
        message(FATAL_ERROR "SinggraphBuild module requires GLPK")
    endif()
endif()

#LIBRARY TO INSTALL
target_link_libraries(${GMDS_LIB} PUBLIC
        ${LIB_GMDS_IG}
        ${LIB_GMDS_IG_ALGO}
        ${LIB_GMDS_IO}
        ${LIB_GMDS_CAD}
        ${LIB_GMDS_SMOOTHY}
        Predicates_psm
        INTERFACE glpk)

target_compile_features(${GMDS_LIB} PUBLIC cxx_std_14)

# INCLUDE TO INSTALL
target_include_directories(${GMDS_LIB} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
        )
set_target_properties(${GMDS_LIB} PROPERTIES PUBLIC_HEADER "${GMDS_INC}")

#install(TARGETS ${GMDS_LIB}
    #    EXPORT GMDS_SUITE
    #    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${GMDS_LIB_PREFIX}
    #    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gmds
    #    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gmds)




#==============================================================================
if(WITH_TEST)
    add_subdirectory(tst)
endif(WITH_TEST)

#==============================================================================
add_executable(dualblock
        src/main.cpp)
target_link_libraries(dualblock PUBLIC ${LIB_GMDS_DUAL_BLOCKING})
#==============================================================================
