cmake_minimum_required(VERSION 2.8)
project(test1)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated -fPIC")

set(OMIT_TEST False)
add_subdirectory(../../vp_new_try/cxx vp_build)
#set_target_properties(test PROPERTIES EXCLUDE_FROM_ALL TRUE)

include_directories(../../vp_new_try/cxx)

find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)

include_directories(${JSONCPP_INCLUDE_DIRS})

add_executable(test1 Decomposer.cxx test.cxx)
target_link_libraries(test1 test_ ${JSONCPP_LIBRARIES})
