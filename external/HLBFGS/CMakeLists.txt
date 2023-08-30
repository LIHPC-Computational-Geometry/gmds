
#==============================================================================
add_library(hlbfgs STATIC
        HLBFGS.cpp
        HLBFGS.h
        HLBFGS_BLAS.cpp
        HLBFGS_BLAS.h
        ICFS.cpp
        ICFS.h
        LineSearch.cpp
        LineSearch.h
        Lite_Sparse_Matrix.cpp
        Lite_Sparse_Matrix.h
        Sparse_Entry.h
        )

#==============================================================================
# TARGET DEFINITION
set_target_properties(hlbfgs PROPERTIES SOVERSION 1.2)
#==============================================================================
# CPP INTERFACE TO EXPORT
target_include_directories(hlbfgs  PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
        )


install(TARGETS hlbfgs
        EXPORT GMDS_SUITE
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION include
        )
