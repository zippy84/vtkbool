#   Copyright 2012-2018 Ronald Römer
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

cmake_minimum_required(VERSION 2.8)
project(vtkbool)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated -fPIC -Wunused")

option(BUILD_PARAVIEW_PLUGIN "" Off)
option(BUILD_AS_DEBUG "" Off)
option(BUILD_TESTING "" Off)

set(srcs
    vtkPolyDataBooleanFilter.cxx
    vtkPolyDataContactFilter.cxx
    Utilities.cxx)

include_directories(".")

if(BUILD_AS_DEBUG)
    add_definitions(-DDEBUG)
endif()

if(BUILD_PARAVIEW_PLUGIN)
    set(BUILD_TESTING Off CACHE BOOL "" FORCE)

    add_subdirectory(paraview_plugin ${CMAKE_BINARY_DIR}/paraview_plugin)

else()

    #unset(VTK_DIR CACHE)

    find_package(VTK REQUIRED COMPONENTS vtkFiltersSources vtkIOLegacy vtkFiltersExtraction vtkFiltersGeometry vtkFiltersModeling vtkRenderingFreeType vtkWrappingPythonCore NO_MODULE)

    if(VTK_FOUND)
        include(${VTK_USE_FILE})

        add_subdirectory(dev)

        include_directories(dev/merger dev/vp_new_try)

        add_library(${PROJECT_NAME} SHARED ${srcs})
        target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} decomp merger)

        if(BUILD_TESTING)
            enable_testing()

            add_executable(testing testing/testing.cxx)
            target_link_libraries(testing ${PROJECT_NAME} ${VTK_LIBRARIES})

            foreach(i RANGE 16)
                add_test(Test_${i} testing ${i})
            endforeach()

        endif()

        if(VTK_WRAP_PYTHON)
            find_package(PythonLibs 2.7 REQUIRED)

            if(PYTHONLIBS_FOUND)

                include_directories(${PYTHON_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})

                include(vtkWrapPython)

                vtk_wrap_python3(${PROJECT_NAME}Python ${PROJECT_NAME}Python_srcs vtkPolyDataBooleanFilter.cxx)

                add_library(${PROJECT_NAME}PythonD ${${PROJECT_NAME}Python_srcs})
                target_link_libraries(${PROJECT_NAME}PythonD ${PROJECT_NAME} vtkWrappingPythonCore vtkCommonExecutionModelPythonD)

                python_add_module(${PROJECT_NAME}Python ${PROJECT_NAME}PythonInit.cxx)

            endif()
        endif()

    endif(VTK_FOUND)

endif()
