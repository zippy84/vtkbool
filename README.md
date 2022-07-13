# vtkbool [![travis](https://app.travis-ci.com/zippy84/vtkbool.svg?branch=master)](https://app.travis-ci.com/zippy84/vtkbool) [![codecov](https://codecov.io/gh/zippy84/vtkbool/branch/master/graph/badge.svg?token=EUV9QKEW1M)](https://codecov.io/gh/zippy84/vtkbool)

![](/cover.png)

## About

This is an extension of the graphics library VTK. The goal of the extension is to equip the library with boolean operations on polygonal meshes. I started the project at the end of my studies in mechanical engineering at the University of Applied Sciences ([HTWK](http://htwk-leipzig.de/)) in Leipzig. I used VTK to develop a program, which I had to create for a paper. At this time I would have wished, that this feature already exists. There was several implementations from third parties, but after some tests, I came to the conclusion, that none of them worked correct. I decided to start with my own implementation. This library is the result of my efforts.

## Features

- no extra libraries required
- 4 operation types available (union, intersection, difference and difference2 - difference with interchanged operands)
- triangulation is not needed
- all types of polygonal cells are supported (triangles, quads, polygons, triangle-strips)
- triangle-strips and quads will be transformed into triangles (quads only if points are not on the same plane)
- non-convex polygons are allowed
- meshes can be stacked (coplanar polygons are right handled)
- the meshes don’t need to be closed
- CellData is passed (attached by the rules of vtkAppendPolyData)
- contact-lines are available in the 3th output
- the filter is able to embed holes
- compileable as ParaView plugin
- Python wrapped

## Limitations

- the filter assumes well defined triangles, quads and polygons
- PointData is not preserved - you have to do your own mapping (use *OrigCellIdsA* and *OrigCellIdsB*)

## Requirements

- CMake >= 3.12
- VTK >= 9.0
- C++11 compiler

### Optional

- ParaView >= 5.0
- Python 3.x

## Library

To include vtkbool into your program, you have to compile it as a library. All you need is an installation of VTK with header files. If you have installed VTK over your package manager, CMake is able to find the required files. Otherwise you have to set **VTK\_DIR** manually. It must be a path like */home/zippy/VTK9/lib/cmake/vtk-9.1* or *C:/Users/zippy/VTK9/lib/cmake/vtk-9.1*.

The usage of the library is very simple. Look at the example in the section below. You can set the operation mode by calling one of the named methods:

- `SetOperModeToNone`
- `SetOperModeToUnion`
- `SetOperModeToIntersection`
- `SetOperModeToDifference`
- `SetOperModeToDifference2`

The alternative is the more generic `SetOperMode`. The method must be called with the number of the desired operation, an integer between 0 and 4, with the same meaning as mentioned before. The default is Union.

### Example

Create a directory somewhere in your file system, download vtkbool and unpack it into that. Then create the following two files:

**test.cxx**

```C++
#include <vtkSmartPointer.h>
#include <vtkCubeSource.h>
#include <vtkCylinderSource.h>
#include <vtkPolyDataWriter.h>

#include "vtkPolyDataBooleanFilter.h"

int main (int argc, char *argv[]) {
    vtkSmartPointer<vtkCubeSource> cube = vtkSmartPointer<vtkCubeSource>::New();
    cube->SetYLength(.5);

    vtkSmartPointer<vtkCylinderSource> cyl = vtkSmartPointer<vtkCylinderSource>::New();
    cyl->SetResolution(32);
    cyl->SetHeight(.5);
    cyl->SetCenter(0, .5, 0);

    vtkSmartPointer<vtkPolyDataBooleanFilter> bf = vtkSmartPointer<vtkPolyDataBooleanFilter>::New();
    bf->SetInputConnection(0, cube->GetOutputPort());
    bf->SetInputConnection(1, cyl->GetOutputPort());
    bf->SetOperModeToDifference();

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputConnection(bf->GetOutputPort());
    writer->SetFileName("result.vtk");
    writer->Update();

    return 0;
}
```

**CMakeLists.txt**

```CMake
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(test)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# find_package(VTK REQUIRED COMPONENTS FiltersSources IOLegacy)

find_package(VTK REQUIRED COMPONENTS FiltersSources IOLegacy FiltersExtraction FiltersGeometry FiltersModeling FiltersFlowPaths WrappingPythonCore)

if(VTK_FOUND)
    include_directories(vtkbool-master)
    add_subdirectory(vtkbool-master)

    add_executable(test test.cxx)
    target_link_libraries(test PRIVATE vtkBool ${VTK_LIBRARIES})

    vtk_module_autoinit(
        TARGETS test
        MODULES ${VTK_LIBRARIES}
    )
endif(VTK_FOUND)

```

Also create a directory named build. If you are on Linux, open a terminal and enter this directory. Run `ccmake ..`, follow the instructions, and finally type `make`.

## ParaView Plugin

To build the plugin you have to compile ParaView from source. Download the current version from <http://www.paraview.org> and follow the compilation instructions. As soon as ParaView is compiled, it may take a while, you can build the plugin by activating the **VTKBOOL_PARAVIEW** option within CMake. In CMake you also have to point to **ParaView_DIR** if CMake can't found it and it is not installed in a common location like */usr/lib* or */usr/local/lib*. Make sure **PARAVIEW_INSTALL_DEVELOPMENT_FILES** is set.

When everything has been compiled successfully, you can install the plugin.

## Python

The Python module will be generated automatically, if three conditions are met:

- vtkbool is configured as a library
- Python 3 is installed with header files
- VTK itself is wrapped to Python

After a successful compilation, the module can be used as follows:

```python
import sys
sys.path.append('/path/to/your/build/directory') # also look at the python files of the testing directory

from vtkmodules.vtkFiltersSources import vtkCubeSource, vtkSphereSource
from vtkmodules.vtkIOLegacy import vtkPolyDataWriter
from vtkBool import vtkPolyDataBooleanFilter

cube = vtkCubeSource()

sphere = vtkSphereSource()
sphere.SetCenter(.5, .5, .5)
sphere.SetThetaResolution(20)
sphere.SetPhiResolution(20)

boolean = vtkPolyDataBooleanFilter()
boolean.SetInputConnection(0, cube.GetOutputPort())
boolean.SetInputConnection(1, sphere.GetOutputPort())
boolean.SetOperModeToDifference()

# write the result, if you want ...

writer = vtkPolyDataWriter()
writer.SetInputConnection(boolean.GetOutputPort())
writer.SetFileName('result.vtk')

writer.Update()
```

## Copyright

2012-2022 Ronald Römer

## License

Apache License, Version 2.0
