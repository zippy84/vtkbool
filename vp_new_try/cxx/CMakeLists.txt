cmake_minimum_required(VERSION 2.8)
project(test)

cmake_policy(SET CMP0037 OLD)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated -fPIC")

# pkg-config --cflags --libs jsoncpp

# https://cmake.org/cmake/help/v3.0/module/FindPkgConfig.html
find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)

include_directories(${JSONCPP_INCLUDE_DIRS} .)
add_executable(test test2.cxx Tools.cxx RmTrivials.cxx VisPoly.cxx)
target_link_libraries(test ${JSONCPP_LIBRARIES})
