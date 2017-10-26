cmake_minimum_required(VERSION 2.8)
project(decomp)

if (MSVC)
    set(CMAKE_CXX_FLAGS "-std=c++11 /EHsc")
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
else ()
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated -fPIC -Wunused")
endif ()

add_subdirectory(../vp_new_try vp_build)

include_directories(../vp_new_try)
add_library(decomp STATIC Decomposer.cxx)
target_link_libraries(decomp vp)

if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    cmake_policy(SET CMP0037 OLD)

    # pkg-config --cflags --libs jsoncpp
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(JSONCPP jsoncpp)

    include_directories(${JSONCPP_INCLUDE_DIRS})

    add_executable(test test.cxx)
    target_link_libraries(test decomp ${JSONCPP_LIBRARIES})

endif ()
