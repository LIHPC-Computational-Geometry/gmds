add_library(glpk SHARED IMPORTED GLOBAL)

# CMake doesn't allow to install IMPORTED libraries as TARGETS
# -> https://stackoverflow.com/questions/41175354/can-i-install-shared-imported-library

if(WIN32) # always build shared on windows for now
    set(GLPK_INC ${CMAKE_SOURCE_DIR}/external/glpk/src/glpk/src)
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
    file(COPY "${GLPK_LIB}/glpk_4_65.dll" DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/RelWithDebInfo)

    add_dependencies(glpk buildGlpkMSVC)
    set_target_properties(glpk PROPERTIES
                 IMPORTED_LOCATION "${GLPK_LIB}/glpk_4_65.dll" 
                 IMPORTED_IMPLIB "${GLPK_LIB}/glpk_4_65.lib" 
                 INTERFACE_INCLUDE_DIRECTORIES ${GLPK_INC})
    install(FILES "${GLPK_LIB}/glpk_4_65.dll" DESTINATION  ${CMAKE_INSTALL_BINDIR})
    install(FILES "${GLPK_LIB}/glpk_4_65.lib" DESTINATION  ${CMAKE_INSTALL_LIBDIR})
    install(FILES "${GLPK_INC}/glpk.h" DESTINATION  ${CMAKE_INSTALL_INCLUDEDIR})

else() # unix and macos

    if(BUILD_SHARED_LIBS)
    	if(APPLE)
            set(GLPK_LIBNAME "libglpk.dylib")
	    else()
	        set(GLPK_LIBNAME "libglpk.so")
        endif()
    else()
        set(GLPK_LIBNAME "libglpk.a")
    endif()
    set(GLPK_INC ${CMAKE_SOURCE_DIR}/external/glpk/src/glpk/src)
    message(${GLPK_INC})

    set_target_properties(glpk PROPERTIES
             IMPORTED_LOCATION ${GLPK_LIB}/${GLPK_LIBNAME}
             INTERFACE_INCLUDE_DIRECTORIES ${GLPK_INC})

    install(FILES ${GLPK_LIB}/${GLPK_LIBNAME} DESTINATION  ${CMAKE_INSTALL_LIBDIR})
    install(FILES "${GLPK_INC}/glpk.h" DESTINATION  ${CMAKE_INSTALL_INCLUDEDIR})

    message(STATUS "============= GLPK  ===============")
    message(STATUS "lib path = ${GLPK_LIB}")
    message(STATUS "inc path = ${GLPK_INC}")
    message(STATUS "===================================")

endif()
