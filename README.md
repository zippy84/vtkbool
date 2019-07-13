# vtkbool [![Travis CI Build Status](https://travis-ci.com/zippy84/vtkbool.svg?branch=master)](https://travis-ci.com/zippy84/vtkbool) [![Appveyor Build Status](https://ci.appveyor.com/api/projects/status/github/zippy84/vtkbool?branch=master&svg=true)](https://ci.appveyor.com/project/zippy84/vtkbool) [![Coverage Status](https://codecov.io/gh/zippy84/vtkbool/branch/master/graph/badge.svg)](https://codecov.io/gh/zippy84/vtkbool)

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

To build the plugin you have to compile ParaView from source. Download the current version from <http://www.paraview.org> and follow the compilation instructions. As soon as ParaView is compiled, it may take a while, you can build the plugin by activating the **VTKBOOL_PARAVIEW** option within CMake. In CMake you also have to point to **ParaView_DIR** if CMake can't found it and it is not installed in a common location like */usr/lib* or */usr/local/lib*. Make sure **PARAVIEW_INSTALL_DEVELOPMENT_FILES** is set.

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
