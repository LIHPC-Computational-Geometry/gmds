file(GLOB SOURCES *.cpp)
add_library(hlbfgs STATIC ${SOURCES})
set_target_properties(hlbfgs PROPERTIES SOVERSION 1.2)
