
# vtkbool

![](https://raw.github.com/zippy84/vtkbool/master/cover.png)

## About

This is an extension of the graphics library VTK. The goal of the extension is to equip the library with boolean operations on polygonal meshes. I started the project at the end of my studies in mechanical engineering at the University of Applied Sciences ([HTWK](http://htwk-leipzig.de/)) in Leipzig. I used VTK to develop a program, which I had to create for a paper. At this time I would have wished, that this feature already exists. There was several implementations from third parties, but after some tests, I came to the conclusion, that none of them worked correct. I decided to start with my own implementation. This library is the result of my efforts.

## Features

- no extra libraries required
- triangulation is not needed
- non-convex polygons are allowed
- only the involved polygons are modified
- meshes can be stacked (coplanar polygons are right handled)
- 4 operation types (union, intersection, difference and difference2 - difference with interchanged operands)
- all types of polygonal cells are supported (triangles, quads, polygons, triangle-strips)
- the meshes don’t need to be closed
- CellData is passed (attached by the rules of vtkAppendPolyData)
- original cell ids are added to CellData (*OrigCellIdsA* and *OrigCellIdsB* as vtkIntArray)
- contact lines are available
- non-convex polygons are decomposed
- the filter is able to embed holes
- it is also a plugin for ParaView
- Python wrapped

## Limitations

- correctness depends on the polygon-defining normals (use vtkPolyDataNormals if you have problems with incorrect orientations)
- PointData is not preserved - you have to do your own mapping

## Todo

- the filter needs a mesh optimizer for sharp-angled triangles

## Requirements

- CMake >= 3.1
- VTK >= 6.1
- GCC >= 4.7 or another compiler with C++11-support (MSVC for example)

### Optional

- ParaView >= 5.0
- Python 2.7 or 3.x

## Library

To include vtkbool into your program, you have to compile it as a library. All you need is an installation of VTK with header files. If you have installed VTK over your package manager, CMake is able to find the required files. Otherwise you have to set **VTK\_DIR** manually. It must be a path like */home/zippy/VTK6/lib/cmake/vtk-6.3* or *C:/Users/zippy/VTK6/lib/cmake/vtk-6.3*.

The usage of the library is very simple. Look at *testing.cxx* and you can see how. Upon creating the instance of the boolean-filter and connecting the two inputs with the pipeline, you can choose between different operation types, in different manners. You can set the operation mode by calling one of the named methods:

- `SetOperModeToUnion`
- `SetOperModeToIntersection`
- `SetOperModeToDifference`
- `SetOperModeToDifference2`

The alternative is the more generic `SetOperMode`. The method must be called with the number of the desired operation, an integer between 0 and 3, with the same meaning as mentioned before. After updating the pipeline, the result is stored in the first output, typically accessable with `GetOutputPort()`. The second output, `GetOutputPort(1)`, contains the lines of contact between the inputs. The inputs must be outputs of filters or sources returning vtkPolyData. The outputs from this filter are of the same type.

Other options are `MergeRegs` and `DecPolys`.

The first option is used in testing and is deactivated by default. It ignors the `OperMode` and mergs all divided regions into the output. If you don't want to combine the regions yourself, don't use it.

The second option controls whether non-convex polygons will be decomposed into convex polygons. Only the created polygons will be decomposed. The option is activated by default and it has to stay activated, if you want to triangulate the mesh with `vtkTriangleFilter`.

### Example

Create a directory somewhere in your file system, download vtkbool and unpack it into that. Then create the following two files:

**exp.cxx**

```C++
#include <vtkCubeSource.h>
#include <vtkCylinderSource.h>

#include "vtkPolyDataBooleanFilter.h"

int Point::_tag = 0; // important

int main (int argc, char *argv[]) {
    vtkCubeSource *cu = vtkCubeSource::New();
    cu->SetYLength(.5);

    vtkCylinderSource *cyl = vtkCylinderSource::New();
    cyl->SetResolution(32);
    cyl->SetHeight(.5);
    cyl->SetCenter(0, .5, 0);

    vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
    bf->SetInputConnection(0, cu->GetOutputPort());
    bf->SetInputConnection(1, cyl->GetOutputPort());
    bf->SetOperModeToDifference();
    bf->Update();

    bf->Delete();
    cyl->Delete();
    cu->Delete();

    return 0;
}
```
**CMakeLists.txt**

```CMake
cmake_minimum_required(VERSION 3.1)
project(exp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(vtkbool-master)

add_subdirectory(vtkbool-master)

find_package(VTK REQUIRED COMPONENTS vtkFiltersSources NO_MODULE)
include(${VTK_USE_FILE})

include_directories(vtkbool-master/libs/vp)

add_executable(exp exp.cxx)
target_link_libraries(exp ${VTK_LIBRARIES} vtkbool)
```
Also create a directory named build. If you are on Linux, open a terminal and enter this directory. Run `ccmake ..`, follow the instructions, and finally type `make`.

## ParaView Plugin

To build the plugin you have to compile ParaView from source. Download the current version from <http://www.paraview.org> and follow the compilation instructions. As soon as ParaView is compiled, it may take a while, you can build the plugin by activating the **BUILD_PARAVIEW** option within CMake. In CMake you also have to point to **ParaView_DIR** if CMake can't found it and it is not installed in a common location like */usr/lib* or */usr/local/lib*. Make sure **PARAVIEW_INSTALL_DEVELOPMENT_FILES** is set.

When everything has been compiled successfully, you can install the plugin. For that purpose I made screenshots from the necessary steps. You can find them under *paraview\_plugin/install/*. There is also a screenshot of how to use it.

## Python

The Python module will be generated automatically, if three conditions are met:

- vtkbool is configured as a library
- Python 2 or 3 is installed with header files
- VTK itself is wrapped to Python

After a successful compilation, the module can be used as follows:

```python
import sys
sys.path.append('/path/to/your/build/directory')

import vtk
import vtkboolPython

cube = vtk.vtkCubeSource()

sphere = vtk.vtkSphereSource()
sphere.SetCenter(.5, .5, .5)
sphere.SetThetaResolution(20)
sphere.SetPhiResolution(20)

boolean = vtkboolPython.vtkPolyDataBooleanFilter()
boolean.SetInputConnection(0, cube.GetOutputPort())
boolean.SetInputConnection(1, sphere.GetOutputPort())
boolean.SetOperModeToDifference()

# write the result, if you want ...

writer = vtk.vtkPolyDataWriter()
writer.SetInputConnection(boolean.GetOutputPort())
writer.SetFileName('result.vtk')

writer.Update()
```

There are two complex examples in *python/examples*.

## Donating

If you like this project and you want to support it, you can donate with [PayPal](https://paypal.me/zippy84).

## Copyright

2012-2019 Ronald Römer

## License

Apache License, Version 2.0


********************************************************

Use CMake to build and package vtkbool.

Get vtkbool codes from repository.
from vtbool folder right click to open vs2019

select vtkbool\CMakeLists.txt and right click "Generate Cache for vtkbool".
    1> [CMake] -- Build files have been written to: C:/Workarea/vtkbool/out/build/x64-Debug
	
select vtkbool\CMakeLists.txt (two times) and right click "Build" to make sure the original vtkbool codes could be built successfully! 

  [13/13] cmd.exe /C "cmd.exe /C ""C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" -E __create_def C:\Workarea\vtkbool\out\build\x64-Debug\CMakeFiles\vtkbool.dir\exports.def C:\Workarea\vtkbool\out\build\x64-Debug\CMakeFiles\vtkbool.dir\exports.def.objs && cd C:\Workarea\vtkbool\out\build\x64-Debug" && "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" -E vs_link_dll --intdir=CMakeFiles\vtkbool.dir --rc=C:\PROGRA~2\WI3CF2~1\10\bin\100183~1.0\x64\rc.exe --mt=C:\PROGRA~2\WI3CF2~1\10\bin\100183~1.0\x64\mt.exe --manifests  -- C:\PROGRA~2\MICROS~2\2019\PROFES~1\VC\Tools\MSVC\1422~1.279\bin\Hostx64\x64\link.exe /nologo @CMakeFiles\vtkbool.rsp  /out:vtkbool.dll /implib:vtkbool.lib /pdb:vtkbool.pdb /dll /version:0.0 /machine:x64  /debug /INCREMENTAL  /DEF:CMakeFiles\vtkbool.dir\exports.def   && cmd.exe /C "cd /D C:\Workarea\vtkbool\out\build\x64-Debug && powershell -noprofile -executionpolicy Bypass -file C:/Workarea/vcpkg/scripts/buildsystems/msbuild/applocal.ps1 -targetBinary C:/Workarea/vtkbool/out/build/x64-Debug/vtkbool.dll -installedDir C:/Workarea/vcpkg/installed/x64-windows/debug/bin -OutVariable out""
     Creating library vtkbool.lib and object vtkbool.exp
     Creating library vtkbool.lib and object vtkbool.exp
Build succeeded.	 
	 
The vtkbool has self codes and three dependencies codes.

Add install target command at vtkbool and dependencies CMakeLists.txt

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# for dependency install target at CMakeLists.txt

set(INSTALL_INCLUDE_DIR include CACHE PATH
	"Installation directory for header files")
	
set_target_properties(dependency PROPERTIES
  PUBLIC_HEADER "dependency.h")

install(TARGETS dependency
  # IMPORTANT: Add the foo library to the "export-set"
  EXPORT vtkboolTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/vtkbool/libs/dependency"
    COMPONENT dev)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# for vtkbool install target  at CMakeLists.txt

set_target_properties(vtkbool PROPERTIES
  PUBLIC_HEADER vtkPolyDataBooleanFilter.h vtkPolyDataContactFilter.h Utilities.h)

install(TARGETS vtkbool
  # IMPORTANT: Add the vtkbool library to the "export-set"
  EXPORT vtkboolTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}" 
    COMPONENT dev)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# for vtkbool pakcage  at CMakeLists.txt

	set(VTKBOOL_MAJOR_VERSION 0)
	set(VTKBOOL_MINOR_VERSION 1)
	set(VTKBOOL_PATCH_VERSION 0)
	set(VTKBOOL_VERSION
	  ${VTKBOOL_MAJOR_VERSION}.${VTKBOOL_MINOR_VERSION}.${VTKBOOL_PATCH_VERSION})

	# Offer the user the choice of overriding the installation directories
	set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
	set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
	set(INSTALL_INCLUDE_DIR include CACHE PATH
	  "Installation directory for header files")
	if(WIN32 AND NOT CYGWIN)
	  set(DEF_INSTALL_CMAKE_DIR CMake)
	else()
	  set(DEF_INSTALL_CMAKE_DIR lib/CMake/vtkbool)
	endif()
	set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
	  "Installation directory for CMake files")

	# Make relative paths absolute (needed later on)
	foreach(p LIB BIN INCLUDE CMAKE)
	  set(var INSTALL_${p}_DIR)
	  if(NOT IS_ABSOLUTE "${${var}}")
		set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
	  endif()
	endforeach()


	# set up include-directories
	include_directories(
	  "${PROJECT_SOURCE_DIR}"   # to find vtkbool .h files
	  "${PROJECT_BINARY_DIR}")  # 


    #add_subdirectory()

		# The interesting stuff goes here
	# ===============================

	# Add all targets to the build-tree export set
	
	export(TARGETS vtkbool dependency
	  FILE "${PROJECT_BINARY_DIR}/vtkboolTargets.cmake")

	# Export the package for use from the build-tree
	# (this registers the build-tree with a global CMake-registry)
	export(PACKAGE vtkbool)

		# Create the vtkboolConfig.cmake and vtkboolVersionConfig files
	file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
	   "${INSTALL_INCLUDE_DIR}")
	# ... for the build tree
	set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/libs/vp" "${PROJECT_BINARY_DIR}/libs/vp")
	configure_file(vtkbool.cmake.in
	  "${PROJECT_BINARY_DIR}/vtkboolConfig.cmake" @ONLY)
	# ... for the install tree
	set(CONF_INCLUDE_DIRS "\${VTKBOOL_CMAKE_DIR}/${REL_INCLUDE_DIR}")
	configure_file(vtkbool.cmake.in
	  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/vtkboolConfig.cmake" @ONLY)
	# ... for both
	configure_file(vtkboolVersion.cmake.in
	  "${PROJECT_BINARY_DIR}/vtkboolVersionConfig.cmake" @ONLY)


	# Install the vtkbool.cmake and vtkboolVersion.cmake
	install(FILES
	  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/vtkboolConfig.cmake"
	  "${PROJECT_BINARY_DIR}/vtkboolVersionConfig.cmake"
	  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

	#Install the export set for use with the install-tree
	install(EXPORT vtkboolTargets DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	export(TARGETS vtkbool dependency FILE "${PROJECT_BINARY_DIR}/vtkboolTargets.cmake")
	    will export vtkbool and dependency to create vtkboolTargets.cmake.
			
    export(PACKAGE vtkbool) will registers the build-tree with a global CMake-registry at
	local computer Registry Editor 
	  Computer\HKEY_CURRENT_USER\Software\Kitware\CMake\Packages\VTKBool
	      Data: C:/Workarea/vtkbool/out/build/x64-Debug
	
    From this location the user will using find_package command to get vtkboolConfig.cmake file to get vtkbool package information.
	
	configure_file(vtkbool.cmake.in "${PROJECT_BINARY_DIR}/vtkboolConfig.cmake" @ONLY)
	    will export vtkboolConfig.make based on vtkbool.cmake.in file information.
		
		# - Config file for the VTKBOOL package
		# It defines the following variables
		#  VTKBOOL_INCLUDE_DIRS - include directories for MeshFix
		#  VTKBOOL_LIBRARIES    - libraries to link against

		# Compute paths
		get_filename_component(VTKBOOL_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
		set(VTKBOOL_INCLUDE_DIRS "@CONF_INCLUDE_DIRS@")
		include("${VTKBOOL_CMAKE_DIR}/VTKBoolTargets.cmake")
		set(VTKBOOL_LIBRARIES VTKBOOL)
		
    		
		install(FILES ...) and install(EXPORT ...) will create install lib and header files to the target machine.
		
    After all done and saved,
	
	select vtkbool\CMakeLists.txt and right click "Generate Cache for vtkbool".
		1> [CMake] -- Build files have been written to: C:/Workarea/vtkbool/out/build/x64-Debug
		
	select vtkbool\CMakeLists.txt (two times) and right click "Build" 

	  [13/13] cmd.exe /C "cmd.exe /C ""C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" -E __create_def C:\Workarea\vtkbool\out\build\x64-Debug\CMakeFiles\vtkbool.dir\exports.def C:\Workarea\vtkbool\out\build\x64-Debug\CMakeFiles\vtkbool.dir\exports.def.objs && cd C:\Workarea\vtkbool\out\build\x64-Debug" && "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" -E vs_link_dll --intdir=CMakeFiles\vtkbool.dir --rc=C:\PROGRA~2\WI3CF2~1\10\bin\100183~1.0\x64\rc.exe --mt=C:\PROGRA~2\WI3CF2~1\10\bin\100183~1.0\x64\mt.exe --manifests  -- C:\PROGRA~2\MICROS~2\2019\PROFES~1\VC\Tools\MSVC\1422~1.279\bin\Hostx64\x64\link.exe /nologo @CMakeFiles\vtkbool.rsp  /out:vtkbool.dll /implib:vtkbool.lib /pdb:vtkbool.pdb /dll /version:0.0 /machine:x64  /debug /INCREMENTAL  /DEF:CMakeFiles\vtkbool.dir\exports.def   && cmd.exe /C "cd /D C:\Workarea\vtkbool\out\build\x64-Debug && powershell -noprofile -executionpolicy Bypass -file C:/Workarea/vcpkg/scripts/buildsystems/msbuild/applocal.ps1 -targetBinary C:/Workarea/vtkbool/out/build/x64-Debug/vtkbool.dll -installedDir C:/Workarea/vcpkg/installed/x64-windows/debug/bin -OutVariable out""
		 Creating library vtkbool.lib and object vtkbool.exp
		 Creating library vtkbool.lib and object vtkbool.exp
	Build succeeded.	 
	
    select vtkbool\CMakeLists.txt (two times) and right click "Install"
	
>------ Build started: Project: CMakeLists, Configuration:  ------
  [0/1] cmd.exe /C "cd /D C:\Workarea\vtkbool\out\build\x64-Debug && "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" -P cmake_install.cmake"
  -- Install configuration: "Debug"
  -- Up-to-date: C:/Workarea/vtkbool/out/install/x64-Debug/lib/vtkbool.lib
  -- Up-to-date: C:/Workarea/vtkbool/out/install/x64-Debug/bin/vtkbool.dll
  -- Up-to-date: C:/Workarea/vtkbool/out/install/x64-Debug/include/vtkPolyDataBooleanFilter.h
  -- Up-to-date: C:/Workarea/vtkbool/out/install/x64-Debug/CMake/vtkboolConfig.cmake
  -- Up-to-date: C:/Workarea/vtkbool/out/install/x64-Debug/CMake/vtkboolVersionConfig.cmake
  -- Up-to-date: C:/Workarea/vtkbool/out/install/x64-Debug/CMake/vtkboolTargets.cmake
  -- Up-to-date: C:/Workarea/vtkbool/out/install/x64-Debug/CMake/vtkboolTargets-debug.cmake
  -- Up-to-date: C:/Workarea/vtkbool/out/install/x64-Debug/lib/decomp.lib
  -- Up-to-date: C:/Workarea/vtkbool/out/install/x64-Debug/include/libs/decomp/Decomposer.h
  -- Up-to-date: C:/Workarea/vtkbool/out/install/x64-Debug/lib/vp.lib
  -- Up-to-date: C:/Workarea/vtkbool/out/install/x64-Debug/include/libs/vp/Tools.h
  -- Up-to-date: C:/Workarea/vtkbool/out/install/x64-Debug/lib/merger.lib
  -- Up-to-date: C:/Workarea/vtkbool/out/install/x64-Debug/include/libs/merger/Merger.h
Install succeeded.


    Open another application and add find_package(vtkbool REQUIRED) command to look for vtkbool at register Computer\HKEY_CURRENT_USER\Software\Kitware\CMake\Packages\VTKBool location file vtkboolConfig.camke
	
	find_package(vtkbool REQUIRED)
	include_directories(${VTKBOOL_INCLUDE_DIRS})
	include_directories(${VTKBOOL_INCLUDE_DIRS}/libs/vp)

	add_executable(application MACOSX_BUNDLE application.cxx) 
    target_link_libraries(application vtkbool)

    select application\CMakeLists.txt and right click "Generate Cache for application".
	select application\CMakeLists.txt (two times) and right click "Build" 
	