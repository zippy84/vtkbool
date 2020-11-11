# Copyright 2012-2020 Ronald Römer
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
project(vtkbool VERSION 2.9)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)

if(MSVC)
    add_compile_options(/EHsc)
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
else()
    add_compile_options(-Wall -Wextra -Wno-sign-compare -fPIC)
endif()

option(VTKBOOL_PARAVIEW "" OFF)
option(VTKBOOL_DEBUG "" OFF)
option(VTKBOOL_TESTING "" OFF)
option(VTKBOOL_COVERAGE "" OFF)

mark_as_advanced(VTKBOOL_DEBUG)
mark_as_advanced(VTKBOOL_TESTING)
mark_as_advanced(VTKBOOL_COVERAGE)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND VTKBOOL_COVERAGE)
    add_compile_options(-fprofile-arcs -ftest-coverage)
    set(CMAKE_EXE_LINKER_FLAGS "-fprofile-arcs")
    set(CMAKE_SHARED_LINKER_FLAGS "-fprofile-arcs")
endif()

set(srcs
    vtkPolyDataBooleanFilter.cxx
    vtkPolyDataContactFilter.cxx
    Utilities.cxx)

include_directories(".")

if(VTKBOOL_DEBUG)
    add_definitions(-DDEBUG)
endif()

if(VTKBOOL_PARAVIEW)
    set(VTKBOOL_TESTING OFF CACHE BOOL "" FORCE)

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

        if(VTK_INSTALL_PREFIX)
            set(_vtk ${VTK_INSTALL_PREFIX})
        else()
            set(_vtk ${VTK_DIR})
        endif()

        if(WIN32)
            string(REPLACE "/" "\\" _vtk_win ${_vtk})
        endif()

        if(VTKBOOL_TESTING)
            enable_testing()

            add_executable(testing testing/testing.cxx)
            target_link_libraries(testing ${PROJECT_NAME} ${VTK_LIBRARIES})

            foreach(i RANGE 20)
                add_test(NAME Test_${i} COMMAND testing ${i})

                if(WIN32)
                    set_property(TEST Test_${i} PROPERTY ENVIRONMENT "PATH=${_vtk_win}\\bin")
                endif()

            endforeach()

        endif()

        if(VTK_WRAP_PYTHON)

            if(${VTK_MAJOR_VERSION} GREATER 7)
                set(MODULE_HIERARCHY_NAME ${PROJECT_NAME}Hierarchy)
                set(${PROJECT_NAME}_WRAP_DEPENDS vtkCommonExecutionModel)

                set(${PROJECT_NAME}_LINK_DEPENDS vtkCommonExecutionModel)

                include(vtkWrapHierarchy)
                vtk_wrap_hierarchy(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR} vtkPolyDataBooleanFilter.h)
                set(KIT_HIERARCHY_FILE ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_HIERARCHY_NAME}.txt)
                set(LIB_HIERARCHY_STAMP ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_HIERARCHY_NAME}.stamp.txt)

                if(${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION} EQUAL "8.1")
                    configure_file(
                        ${VTK_CMAKE_DIR}/vtkWrapperInit.data.in
                        ${MODULE_HIERARCHY_NAME}.data
                        @ONLY
                    )

                    set(OTHER_HIERARCHY_FILES)

                    foreach(dep ${${PROJECT_NAME}_LINK_DEPENDS})
                        list(APPEND OTHER_HIERARCHY_FILES "${${dep}_WRAP_HIERARCHY_FILE}")
                    endforeach()

                    add_custom_command(OUTPUT ${KIT_HIERARCHY_FILE}
                            ${LIB_HIERARCHY_STAMP}
                        COMMAND ${VTK_WRAP_HIERARCHY_EXE} -o ${KIT_HIERARCHY_FILE}
                                ${MODULE_HIERARCHY_NAME}.data ${OTHER_HIERARCHY_FILES}
                        COMMAND ${CMAKE_COMMAND} -E touch ${LIB_HIERARCHY_STAMP}
                        DEPENDS ${VTK_WRAP_HIERARCHY_EXE}
                            ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_HIERARCHY_NAME}.data
                            ${OTHER_HIERARCHY_FILES}
                        VERBATIM)
                endif()

            endif()

            set(ENV{LD_LIBRARY_PATH} "${_vtk}/lib")

            execute_process(COMMAND "${_vtk}/bin/vtkpython" -c "import sys; sys.stdout.write('.'.join(map(str, sys.version_info[:2])))" OUTPUT_VARIABLE _py)

            if(_py)
                if(CMAKE_VERSION VERSION_LESS "3.12")
                    find_package(PythonInterp ${_py} EXACT)
                    find_package(PythonLibs ${_py} EXACT)

                    if(PYTHONLIBS_FOUND)

                        include_directories(${PYTHON_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})

                        include(vtkWrapPython)

                        vtk_wrap_python3(${PROJECT_NAME}Python ${PROJECT_NAME}Python_srcs vtkPolyDataBooleanFilter.cxx)

                        add_library(${PROJECT_NAME}PythonD ${${PROJECT_NAME}Python_srcs})
                        target_link_libraries(${PROJECT_NAME}PythonD ${PROJECT_NAME} vtkWrappingPythonCore vtkCommonExecutionModelPythonD)

                        python_add_module(${PROJECT_NAME}Python PyInit.cxx ${PROJECT_NAME}PythonInit.cxx)
                        target_link_libraries(${PROJECT_NAME}Python PRIVATE ${PROJECT_NAME}PythonD)

                        if(VTKBOOL_TESTING)
                            add_test(NAME Test_Py COMMAND ${PYTHON_EXECUTABLE}
                                ${CMAKE_SOURCE_DIR}/testing/test_py_module.py
                                ${_vtk}/lib/python${_py}/site-packages
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
                    find_package(Python ${_py} EXACT COMPONENTS Interpreter Development)

                    if(Python_Development_FOUND)

                        include_directories(${PYTHON_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})

                        include(vtkWrapPython)

                        vtk_wrap_python3(${PROJECT_NAME}Python ${PROJECT_NAME}Python_srcs vtkPolyDataBooleanFilter.cxx)

                        add_library(${PROJECT_NAME}PythonD ${${PROJECT_NAME}Python_srcs})
                        target_link_libraries(${PROJECT_NAME}PythonD ${PROJECT_NAME} vtkWrappingPythonCore vtkCommonExecutionModelPythonD)

                        Python_add_library(${PROJECT_NAME}Python MODULE PyInit.cxx ${PROJECT_NAME}PythonInit.cxx)
                        target_link_libraries(${PROJECT_NAME}Python PRIVATE ${PROJECT_NAME}PythonD)

                        if(VTKBOOL_TESTING)
                            add_test(NAME Test_Py COMMAND ${Python_EXECUTABLE}
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

                endif()

            endif()

        endif()

    endif(VTK_FOUND)

endif()
