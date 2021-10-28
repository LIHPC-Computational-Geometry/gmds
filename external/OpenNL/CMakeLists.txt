#==============================================================================
# OpenNL LIBRARY
#==============================================================================
add_library(OpenNL STATIC
        OpenNL_psm.h
        OpenNL_psm.cpp)

#==============================================================================
# TARGET DEFINITION
#==============================================================================
# CPP INTERFACE TO EXPORT
target_include_directories(OpenNL  PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
        )


install(TARGETS OpenNL
        EXPORT GMDS_SUITE
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION include
        )
