cmake_minimum_required(VERSION 2.8)
project(vp)

if (MSVC)
    set(CMAKE_CXX_FLAGS "-std=c++11 /EHsc")
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
else ()
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated -fPIC -Wunused")
endif ()

add_compile_options()

add_library(vp STATIC Tools.cxx RmTrivials.cxx VisPoly.cxx)

if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})

    cmake_policy(SET CMP0037 OLD)

    # pkg-config --cflags --libs jsoncpp

    # https://cmake.org/cmake/help/v3.0/module/FindPkgConfig.html
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(JSONCPP jsoncpp)

    include_directories(${JSONCPP_INCLUDE_DIRS})

    #add_executable(test test.cxx)
    add_executable(test test2.cxx)
    target_link_libraries(test ${JSONCPP_LIBRARIES} vp)

endif ()
