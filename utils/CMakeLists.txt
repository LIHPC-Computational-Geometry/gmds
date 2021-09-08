#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_UTILS})
set(GMDS_LIB_PREFIX gmds/utils)

set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/utils/BitVector.h
        inc/gmds/utils/CommonFlags.h
        inc/gmds/utils/CommonTypes.h
        inc/gmds/utils/Exception.h
        inc/gmds/utils/IndexedVector.h
        inc/gmds/utils/LocalCellTopology.h
        inc/gmds/utils/Log.h
        inc/gmds/utils/LogStream.h
        inc/gmds/utils/Marks32.h
        inc/gmds/utils/Parameters.h
        inc/gmds/utils/SmartVector.h
        inc/gmds/utils/Variable.h
        inc/gmds/utils/VariableManager.h
        inc/gmds/utils/FileSystem.h
        inc/gmds/utils/RandomGenerator.h
        )
set(GMDS_SRC
        src/CommonTypes.cpp
        src/Log.cpp
        src/BitVector.cpp
        src/VariableManager.cpp
        src/FileSystem.cpp
        src/LocalCellTopology.cpp
        src/Parameters.cpp
        src/RandomGenerator.cpp
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
target_link_libraries(${GMDS_LIB} PRIVATE INIH)

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
