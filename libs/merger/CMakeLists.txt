# Copyright 2012-2019 Ronald Römer
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.1)
project(merger)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)

if(MSVC)
    add_compile_options(/EHsc)
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
else()
    add_compile_options(-Wno-deprecated -fPIC -Wunused)
endif()

if(NOT VTK_FOUND AND NOT ParaView_FOUND)
    find_package(VTK REQUIRED COMPONENTS vtkFiltersSources vtkIOLegacy vtkFiltersExtraction vtkFiltersGeometry vtkFiltersModeling vtkRenderingFreeType vtkWrappingPythonCore NO_MODULE)
    include(${VTK_USE_FILE})
endif()

include_directories(../vp ../aabb)

add_library(merger STATIC Merger.cxx ../vp/Tools.cxx)
target_link_libraries(merger ${VTK_LIBRARIES})

if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    cmake_policy(SET CMP0037 OLD)

    # pkg-config --cflags --libs jsoncpp
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(JSONCPP jsoncpp)

    include_directories(${JSONCPP_INCLUDE_DIRS} ../third/inja)

    add_executable(test test.cxx)
    target_link_libraries(test merger ${JSONCPP_LIBRARIES})

endif()
