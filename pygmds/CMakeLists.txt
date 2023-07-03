cmake_minimum_required(VERSION 3.10)
#==============================================================================

find_package(pybind11 REQUIRED)

pybind11_add_module(gmds
        pygmds/binding_math.cpp
        pygmds/binding_mesh.cpp
        pygmds/binding_geometry.cpp
        pygmds/binding_blocking.cpp
        pygmds/gmds_facade.cpp
        )

target_link_libraries(gmds PRIVATE
        ${LIB_GMDS_IG}
        ${LIB_GMDS_CADFAC}
        LIB_GMDS_BLOCKING
        ${LIB_GMDS_IO}
        )
target_compile_definitions(gmds PUBLIC cxx_std_14)

install(TARGETS gmds
        COMPONENT python
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
#==============================================================================
if(WITH_TEST)
    add_subdirectory(tst)
endif(WITH_TEST)