#==============================================================================
add_library(INIH STATIC
        inc/dictionary.h
        inc/iniparser.h
        src/dictionary.cpp
        src/iniparser.cpp
        )

#==============================================================================
# TARGET DEFINITION
#==============================================================================
# CPP INTERFACE TO EXPORT
target_include_directories(INIH  PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
        )


install(TARGETS INIH
        EXPORT GMDS_SUITE
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION include
        )
