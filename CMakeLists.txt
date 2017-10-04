cmake_minimum_required(VERSION 2.8)
project(decomp)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated -fPIC")

add_subdirectory(decomp decomp_build)
add_subdirectory(merger merger_build)
