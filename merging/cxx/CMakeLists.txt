cmake_minimum_required(VERSION 2.8)
project(test)

cmake_policy(SET CMP0037 OLD)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated -fPIC")

# pkg-config --cflags --libs jsoncpp

# https://cmake.org/cmake/help/v3.0/module/FindPkgConfig.html
find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)

find_package(VTK REQUIRED COMPONENTS vtkFiltersSources vtkIOLegacy vtkFiltersExtraction vtkFiltersGeometry vtkFiltersModeling vtkRenderingFreeType NO_MODULE)
include(${VTK_USE_FILE})

include_directories(${JSONCPP_INCLUDE_DIRS} . ../../vp_new_try/cxx)
add_executable(test test.cxx Merger.cxx ../../vp_new_try/cxx/Tools.cxx)
target_link_libraries(test ${JSONCPP_LIBRARIES} ${VTK_LIBRARIES})
