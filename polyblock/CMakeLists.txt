#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_POLYBLOCK})
set(GMDS_LIB_PREFIX gmds/polyblock)

set(GMDS_INC
		inc/gmds/polyblock/Gregson2011.h
		inc/gmds/polyblock/PolycubeToolbox.h
		inc/gmds/polyblock/PolyBlock2D.h
		)
set(GMDS_SRC
		src/Gregson2011.cpp
		src/PolycubeToolbox.cpp
		src/PolyBlock2D.cpp)
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
		${LIB_GMDS_IG_ALGO}
		OpenNL GCO)

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
add_executable(PolyBlocker src/main.cpp)
target_link_libraries(PolyBlocker PRIVATE ${GMDS_LIB})
target_compile_features(PolyBlocker PUBLIC cxx_std_14)
