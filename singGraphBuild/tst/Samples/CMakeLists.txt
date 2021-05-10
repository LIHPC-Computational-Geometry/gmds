cmake_minimum_required(VERSION 2.8.10)

file(GLOB_RECURSE TEST_FILES *.* .)

file(COPY ${TEST_FILES} DESTINATION ${CMAKE_BINARY_DIR}/Samples/)

