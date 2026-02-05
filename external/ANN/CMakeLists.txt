#==============================================================================
add_library(ANN STATIC
	include/ANN/ANN.h
        include/ANN/ANNperf.h
        include/ANN/ANNx.h
        src/ANN.cpp
        src/bd_fix_rad_search.cpp
        src/bd_pr_search.cpp
        src/bd_search.cpp
        src/bd_tree.cpp
        src/bd_tree.h
        src/brute.cpp
        src/kd_dump.cpp
        src/kd_fix_rad_search.cpp
        src/kd_fix_rad_search.h
        src/kd_pr_search.cpp
        src/kd_pr_search.h
        src/kd_search.cpp
        src/kd_search.h
        src/kd_split.cpp
        src/kd_split.h
        src/kd_tree.cpp
        src/kd_tree.h
        src/kd_util.cpp
        src/kd_util.h
        src/perf.cpp
        src/pr_queue.h
        src/pr_queue_k.h
        )

set_target_properties(ANN PROPERTIES PUBLIC_HEADER "include/ANN/ANN.h;include/ANN/ANNperf.h;include/ANN/ANNx.h")

#==============================================================================
# TARGET DEFINITION
#==============================================================================
# CPP INTERFACE TO EXPORT
target_include_directories(ANN  PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
        )


install(TARGETS ANN
        EXPORT GMDS_SUITE
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION include
	PUBLIC_HEADER DESTINATION include/ANN
        )
