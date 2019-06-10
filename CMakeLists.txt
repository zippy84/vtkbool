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
project(vtkbool VERSION 2.6)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    set(CMAKE_CXX_FLAGS "/EHsc")
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
else()
    set(CMAKE_CXX_FLAGS "-Wno-deprecated -fPIC -Wunused")
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

    find_package(VTK REQUIRED COMPONENTS vtkFiltersSources vtkIOLegacy vtkFiltersExtraction vtkFiltersGeometry vtkFiltersModeling vtkRenderingFreeType OPTIONAL_COMPONENTS vtkWrappingPythonCore NO_MODULE)

    if(VTK_FOUND)
        include(${VTK_USE_FILE})

        add_subdirectory(libs)

        include_directories(libs/merger libs/decomp libs/vp libs/aabb)

        add_library(${PROJECT_NAME} SHARED ${srcs})
        target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} decomp merger)

        if (VTK_INSTALL_PREFIX)
            set(_vtk ${VTK_INSTALL_PREFIX})
        else()
            set(_vtk ${VTK_DIR})
        endif()

        if(WIN32)
            string(REPLACE "/" "\\" _vtk ${_vtk_win})
        endif()

        if(BUILD_TESTING)
            enable_testing()

            add_executable(testing testing/testing.cxx)
            target_link_libraries(testing ${PROJECT_NAME} ${VTK_LIBRARIES})

            foreach(i RANGE 19)
                add_test(NAME Test_${i} COMMAND testing ${i})

                if(WIN32)
                    set_property(TEST Test_${i} PROPERTY ENVIRONMENT "PATH=${_vtk_win}\\bin")
                endif()

            endforeach()

        endif()

        if(VTK_WRAP_PYTHON)

            set(ENV{LD_LIBRARY_PATH} "${_vtk}/lib")

            execute_process(COMMAND "${_vtk}/bin/vtkpython" -c "import sys; sys.stdout.write('.'.join(map(str, sys.version_info[:2])))" OUTPUT_VARIABLE py_version)

            if(CMAKE_VERSION VERSION_LESS "3.12")
                find_package(PythonInterp "${py_version}")
                find_package(PythonLibs "${py_version}")

                if(PYTHONLIBS_FOUND)

                    include_directories(${PYTHON_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})

                    include(vtkWrapPython)

                    vtk_wrap_python3(${PROJECT_NAME}Python ${PROJECT_NAME}Python_srcs vtkPolyDataBooleanFilter.cxx)

                    add_library(${PROJECT_NAME}PythonD ${${PROJECT_NAME}Python_srcs})
                    target_link_libraries(${PROJECT_NAME}PythonD ${PROJECT_NAME} vtkWrappingPythonCore vtkCommonExecutionModelPythonD)

                    python_add_module(${PROJECT_NAME}Python PyInit.cxx ${PROJECT_NAME}PythonInit.cxx)
                    target_link_libraries(${PROJECT_NAME}Python PRIVATE ${PROJECT_NAME}PythonD)

                    if(BUILD_TESTING)
                        add_test(NAME Test_Py COMMAND ${PYTHON_EXECUTABLE}
                            ${CMAKE_SOURCE_DIR}/testing/test_py_module.py
                            ${_vtk}/lib/python${py_version}/site-packages
                            ${CMAKE_BINARY_DIR}
                            $<CONFIG>)

                        if(WIN32)
                            set_property(TEST Test_Py PROPERTY ENVIRONMENT "PATH=${_vtk_win}\\bin")
                        else()
                            set_property(TEST Test_Py PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${_vtk}/lib")
                        endif()
                    endif()

                endif()

            else()
                if (py_version MATCHES "^3")
                    find_package(Python3 COMPONENTS Interpreter Development)

                    if(Python3_Development_FOUND)

                        include_directories(${PYTHON_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})

                        include(vtkWrapPython)

                        vtk_wrap_python3(${PROJECT_NAME}Python ${PROJECT_NAME}Python_srcs vtkPolyDataBooleanFilter.cxx)

                        add_library(${PROJECT_NAME}PythonD ${${PROJECT_NAME}Python_srcs})
                        target_link_libraries(${PROJECT_NAME}PythonD ${PROJECT_NAME} vtkWrappingPythonCore vtkCommonExecutionModelPythonD)

                        Python3_add_library(${PROJECT_NAME}Python MODULE PyInit.cxx ${PROJECT_NAME}PythonInit.cxx)
                        target_link_libraries(${PROJECT_NAME}Python PRIVATE ${PROJECT_NAME}PythonD)

                        if(BUILD_TESTING)
                            add_test(NAME Test_Py COMMAND ${Python3_EXECUTABLE}
                                ${CMAKE_SOURCE_DIR}/testing/test_py_module.py
                                ${VTK_PYTHONPATH}
                                ${CMAKE_BINARY_DIR}
                                $<CONFIG>)

                            if(WIN32)
                                set_property(TEST Test_Py PROPERTY ENVIRONMENT "PATH=${_vtk_win}\\bin")
                            else()
                                set_property(TEST Test_Py PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${_vtk}/lib")
                            endif()
                        endif()
                    endif()

                elseif(py_version EQUAL "2.7")
                    find_package(Python2 COMPONENTS Interpreter Development)

                    if(Python2_Development_FOUND)

                        include_directories(${PYTHON_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})

                        include(vtkWrapPython)

                        vtk_wrap_python3(${PROJECT_NAME}Python ${PROJECT_NAME}Python_srcs vtkPolyDataBooleanFilter.cxx)

                        add_library(${PROJECT_NAME}PythonD ${${PROJECT_NAME}Python_srcs})
                        target_link_libraries(${PROJECT_NAME}PythonD ${PROJECT_NAME} vtkWrappingPythonCore vtkCommonExecutionModelPythonD)

                        Python2_add_library(${PROJECT_NAME}Python MODULE PyInit.cxx ${PROJECT_NAME}PythonInit.cxx)
                        target_link_libraries(${PROJECT_NAME}Python PRIVATE ${PROJECT_NAME}PythonD)

                        if(BUILD_TESTING)
                            add_test(NAME Test_Py COMMAND ${Python2_EXECUTABLE}
                                ${CMAKE_SOURCE_DIR}/testing/test_py_module.py
                                ${VTK_PYTHONPATH}
                                ${CMAKE_BINARY_DIR}
                                $<CONFIG>)

                            if(WIN32)
                                set_property(TEST Test_Py PROPERTY ENVIRONMENT "PATH=${_vtk}\\bin")
                            else()
                                set_property(TEST Test_Py PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${_vtk}/lib")
                            endif()
                        endif()
                    endif()

                endif()

            endif()

        endif()

    endif(VTK_FOUND)

endif()
