cmake_minimum_required(VERSION 2.6)

add_executable(ann_test ann_test.cpp rand.cpp)
target_link_libraries(ann_test ann)
install(TARGETS ann_test DESTINATION bin)
