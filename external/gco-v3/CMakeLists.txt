#==============================================================================
add_library(GCO
		inc/block.h
		inc/energy.h
		inc/GCoptimization.h
		inc/graph.h
		inc/LinkedBlockList.h
		inc/GCoptimization.cpp
		inc/graph.cpp
		inc/LinkedBlockList.cpp
		inc/maxflow.cpp
		)

#==============================================================================
# TARGET DEFINITION
#==============================================================================
# CPP INTERFACE TO EXPORT
target_include_directories(GCO  PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
		$<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
		)


install(TARGETS GCO
		EXPORT GMDS_SUITE
		DESTINATION ${CMAKE_INSTALL_LIBDIR}
		INCLUDES DESTINATION include
		)