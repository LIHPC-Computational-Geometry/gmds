#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_CLAIRE})
set(GMDS_LIB_PREFIX gmds/claire)

set(GMDS_INC
		${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
		inc/gmds/claire/Smooth2D.h
		)
set(GMDS_SRC
		src/main.cpp
		src/Smooth2D.cpp
        inc/gmds/claire/Grid_Smooth2D.h src/Grid_Smooth2D.cpp inc/gmds/claire/LevelSet2D.h src/LevelSet2D.cpp)
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
add_executable(Claire src/main.cpp)
target_link_libraries(Claire PRIVATE ${GMDS_LIB})
target_compile_features(Claire PUBLIC cxx_std_14)
