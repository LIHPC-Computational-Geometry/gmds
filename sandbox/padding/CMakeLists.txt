#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_PADDING})
set(GMDS_LIB_PREFIX gmds/padding)

set(GMDS_INC
        inc/gmds/padding/SelectivePadding.h
        )
set(GMDS_SRC
        src/SelectivePadding.cpp
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
        ${LIB_GMDS_IO}
        glpk)

target_compile_features(${GMDS_LIB} PUBLIC cxx_std_14)

# INCLUDE TO INSTALL
target_include_directories(${GMDS_LIB} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
        )
set_target_properties(${GMDS_LIB} PROPERTIES PUBLIC_HEADER "${GMDS_INC}")

#install(TARGETS ${GMDS_LIB}
 #       EXPORT GMDS_SUITE
  #      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${GMDS_LIB_PREFIX}
   #     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gmds
    #    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gmds)


#==============================================================================
if(WITH_TEST)
    add_subdirectory(tst)
endif(WITH_TEST)

#==============================================================================
# EXECUTABLE
#==============================================================================
add_executable(PaddingExe src/main.cpp)
target_link_libraries(PaddingExe PRIVATE ${GMDS_LIB})
target_compile_features(PaddingExe PUBLIC cxx_std_14)
