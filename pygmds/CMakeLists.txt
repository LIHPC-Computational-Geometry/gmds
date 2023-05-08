cmake_minimum_required(VERSION 3.10)
#==============================================================================
find_package(pybind11 REQUIRED)

pybind11_add_module(pygmds src/gmds_wrapper.cpp)

target_link_libraries(pygmds PRIVATE GMDSIg GMDSIo)
target_compile_definitions(pygmds PUBLIC cxx_std_14)

install(TARGETS pygmds
        COMPONENT python
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )


#==============================================================================
if(WITH_TEST)
    add_subdirectory(tst)
endif(WITH_TEST)