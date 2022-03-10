cmake_minimum_required(VERSION 3.12)

include(CMakePrintHelpers)

project(examplegmds
  LANGUAGES CXX
  VERSION "1.0.0"
  DESCRIPTION "Example of how to link with gmds"
)

cmake_print_variables(examplegmds_VERSION)

find_package(gmds)

add_executable(examplegmds examplegmds.cpp)

# other components can be added when necessary, such as GMDSIo or GMDSCad
target_link_libraries(examplegmds
	PRIVATE GMDSIg
)
