cmake_minimum_required (VERSION 3.8)

project ("CMakeProject")

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CONFIGURATION_TYPES Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_compile_options(-mavx)

add_executable (CMakeProject "test.cpp")

target_include_directories(CMakeProject PRIVATE "../package/boost_1_82_0" ${CMAKE_CURRENT_SOURCE_DIR})
