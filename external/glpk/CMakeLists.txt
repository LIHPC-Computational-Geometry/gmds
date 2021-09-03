add_library(glpk SHARED IMPORTED GLOBAL)

set(GLPK_INC ${CMAKE_SOURCE_DIR}/external/glpk/src/glpk/src CACHE PATH "Glpk headers location")

if(WIN32)
    set(GLPK_LIB ${CMAKE_SOURCE_DIR}/external/glpk/src/glpk/w64)
    if(NOT EXISTS "${GLPK_LIB}/glpk_4_65.dll")
        ADD_CUSTOM_TARGET(
            buildGlpkMSVC
            COMMAND cd ${GLPK_LIB}
            COMMAND cmd /c ${GLPK_LIB}/Build_GLPK_with_VC14_DLL.bat
        )
    endif()

    file(COPY "${GLPK_LIB}/glpk_4_65.dll" DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug)
    file(COPY "${GLPK_LIB}/glpk_4_65.dll" DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release)

    add_dependencies(glpk buildGlpkMSVC)
    set_target_properties(glpk PROPERTIES
                 IMPORTED_LOCATION "${GLPK_LIB}/glpk_4_65.dll" 
                 IMPORTED_IMPLIB "${GLPK_LIB}/glpk_4_65.lib" 
                 INTERFACE_INCLUDE_DIRECTORIES ${GLPK_INC})
    install(FILES "${GLPK_LIB}/glpk_4_65.dll" DESTINATION  ${CMAKE_INSTALL_BINDIR})
    install(FILES "${GLPK_LIB}/glpk_4_65.lib" DESTINATION  ${CMAKE_INSTALL_LIBDIR})

else() # unix

    # Issue to compile on MacOS so give the path to lib and inc
    set(GLPK_LIB ${CMAKE_BINARY_DIR}/glpk/lib/libglpk.so CACHE PATH "Glpk lib location")  
    if(NOT EXISTS "${GLPK_LIB}")
        MESSAGE(STATUS "GLPK NOT DEFINED > require to build it as an external project")
        ExternalProject_Add(buildGlpk
                PREFIX ${CMAKE_BINARY_DIR}/external/glpk
                URL ${CMAKE_SOURCE_DIR}/external/glpk-4.65.tar.gz
                PATCH_COMMAND autoreconf -if
                CONFIGURE_COMMAND /${CMAKE_BINARY_DIR}/external/glpk/src/buildGlpk/configure --enable-silent-rules
                CC=${CMAKE_C_COMPILER}
                CXX=${CMAKE_CXX_COMPILER}
                "CXXFLAGS=${EXTERNAL_CXX_FLAGS}"
                "LDFLAGS=${EXTERNAL_LD_FLAGS}"
                --prefix=${CMAKE_BINARY_DIR}/glpk
                ${SHARED_LIB_CONFIG_COMMAND}
                BUILD_COMMAND make V=0
                INSTALL_COMMAND make install)
    endif()
    add_dependencies(glpk buildGlpk)
    set_target_properties(glpk PROPERTIES
             IMPORTED_LOCATION ${GLPK_LIB}
             INTERFACE_INCLUDE_DIRECTORIES ${GLPK_INC})
    install(FILES ${GLPK_LIB} DESTINATION  ${CMAKE_INSTALL_BINDIR})
endif()

install(FILES "${GLPK_INC}/glpk.h" DESTINATION  ${CMAKE_INSTALL_INCLUDEDIR})