#==============================================================================
add_library(Predicates_psm STATIC
        Predicates_psm.h
        Predicates_psm.cpp)
#==============================================================================
# TARGET DEFINITION
#==============================================================================
# CPP INTERFACE TO EXPORT
target_include_directories(Predicates_psm  PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
        )



install(TARGETS Predicates_psm
        EXPORT GMDS_SUITE
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION include
        )