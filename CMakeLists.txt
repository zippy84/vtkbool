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

cmake_minimum_required(VERSION 3.1)
project(vtkbool VERSION 2.5)

if(MSVC)
    set(CMAKE_CXX_FLAGS "-std=c++11 /EHsc")
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
else()
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated -fPIC -Wunused")
endif()

option(BUILD_PARAVIEW_PLUGIN "" OFF)
option(BUILD_AS_DEBUG "" OFF)
option(BUILD_TESTING "" OFF)

mark_as_advanced(BUILD_AS_DEBUG)

set(srcs
    vtkPolyDataBooleanFilter.cxx
    vtkPolyDataContactFilter.cxx
    Utilities.cxx)

include_directories(".")

if(BUILD_AS_DEBUG)
    add_definitions(-DDEBUG)
endif()

if(BUILD_PARAVIEW_PLUGIN)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

    add_subdirectory(paraview_plugin ${CMAKE_BINARY_DIR}/paraview_plugin)

else()

    #unset(VTK_DIR CACHE)

    find_package(VTK REQUIRED COMPONENTS vtkFiltersSources vtkIOLegacy vtkFiltersExtraction vtkFiltersGeometry vtkFiltersModeling vtkRenderingFreeType vtkWrappingPythonCore NO_MODULE)

    if(VTK_FOUND)
        include(${VTK_USE_FILE})

        add_subdirectory(libs)

        include_directories(libs/merger libs/decomp libs/vp)

        add_library(${PROJECT_NAME} SHARED ${srcs})
        target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} decomp merger)

        string(REPLACE "/" "\\" vtk_lib ${VTK_INSTALL_PREFIX})

        if(BUILD_TESTING)
            enable_testing()

            add_executable(testing testing/testing.cxx)
            target_link_libraries(testing ${PROJECT_NAME} ${VTK_LIBRARIES})

            foreach(i RANGE 16)
                add_test(NAME Test_${i} COMMAND testing ${i})

                if(WIN32)
                    set_property(TEST Test_${i} PROPERTY ENVIRONMENT "PATH=${vtk_lib}\\bin")
                endif()

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

                python_add_module(${PROJECT_NAME}Python PyInit.cxx ${PROJECT_NAME}PythonInit.cxx)
                target_link_libraries(${PROJECT_NAME}Python ${PROJECT_NAME}PythonD)

                find_package(PythonInterp 2.7 REQUIRED)

                if(PYTHONINTERP_FOUND)
                    add_test(NAME Test_Py COMMAND ${PYTHON_EXECUTABLE}
                        ${CMAKE_SOURCE_DIR}/testing/test_py_module.py
                        ${VTK_INSTALL_PREFIX}/lib/python2.7/site-packages
                        ${CMAKE_BINARY_DIR}
                        $<CONFIG>)

                    if(WIN32)
                        set_property(TEST Test_Py PROPERTY ENVIRONMENT "PATH=${vtk_lib}\\bin")
                    else()
                        set_property(TEST Test_Py PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${VTK_INSTALL_PREFIX}/lib")
                    endif()
                endif()

            endif()
        endif()

    endif(VTK_FOUND)

endif()
