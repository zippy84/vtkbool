cmake_minimum_required(VERSION 2.8)
project(merger)

if (MSVC)
    set(CMAKE_CXX_FLAGS "-std=c++11 /EHsc")
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
else ()
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated -fPIC -Wunused")
endif ()

include_directories(../vp_new_try)

if (NOT DEFINED VTK_DIR)
    find_package(VTK REQUIRED COMPONENTS vtkFiltersSources vtkIOLegacy vtkFiltersExtraction vtkFiltersGeometry vtkFiltersModeling NO_MODULE)
    include(${VTK_USE_FILE})
endif ()

add_library(merger STATIC Merger.cxx ../vp_new_try/Tools.cxx)
target_link_libraries(merger ${VTK_LIBRARIES})

if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    cmake_policy(SET CMP0037 OLD)

    # pkg-config --cflags --libs jsoncpp
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(JSONCPP jsoncpp)

    include_directories(${JSONCPP_INCLUDE_DIRS})

    add_executable(test test.cxx)
    target_link_libraries(test merger ${JSONCPP_LIBRARIES})

endif ()
