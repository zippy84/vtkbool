# Copyright 2012-2018 Ronald Römer
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

cmake_minimum_required(VERSION 2.8)
project(merger)

if(MSVC)
    set(CMAKE_CXX_FLAGS "-std=c++11 /EHsc")
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
else()
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated -fPIC -Wunused")
endif()

find_package(VTK REQUIRED COMPONENTS vtkFiltersSources vtkIOLegacy vtkFiltersExtraction vtkFiltersGeometry vtkFiltersModeling NO_MODULE)
include(${VTK_USE_FILE})

if(VTK_FOUND)
    include_directories(../vp)

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

endif()
