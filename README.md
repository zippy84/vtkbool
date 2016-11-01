# vtkbool

![](https://raw.github.com/zippy84/vtkbool/master/cover.png)

## About

This is an extension of the graphics library VTK. The goal of the extension is to equip the library with boolean operations on polygonal meshes. I started the project at the end of my studies in mechanical engineering at the University of Applied Sciences ([HTWK](http://htwk-leipzig.de/)) in Leipzig. I used VTK to develop a program, which I had to create for a paper. At this time I would have wished, that this feature already exists. There was several implementations from third parties, but after some tests, I came to the conclusion, that none of them worked correct. The once who worked was bound to other big libraries like [CGAL](http://www.cgal.org/). I decided to start with my own implementation. This library is the result of my efforts.

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
- non-convex polygons can be decomposed
- it is also a plugin for ParaView
- Python wrapped

## Limitations

- the filter produces non-convex polygons (the visualization can be broken)
- PointData is not preserved - you have to do your own mapping
- correctness depends on the polygon-defining normals (use vtkPolyDataNormals if you have problems with incorrect orientations)
- when polygons are cut inside, holes are formed - the filter is not able to embed holes

## Todo

- integration of holes
- acceleration with OpenMP

## Build Requirements

- CMake >= 2.8
- VTK >= 6.1
- GCC >= 4.7 or another compiler with C++11-support

### Optional

- ParaView >= 5.0
- Python 2.x

## ParaView Plugin

To build the plugin you have to compile ParaView from source. Download the current version from <http://www.paraView.org> and follow the compilation instructions. As soon as ParaView is compiled, it may take a while, you can build the plugin with the following commands:

	git clone git://github.com/zippy84/vtkbool.git vtkbool
	cd vtkbool
	mkdir build
	cd build
	ccmake ..
	make

Once you have called CMake and typed c for the first time, you will get an error. The error says, that CMake could not found ParaView. You have to specify the directory in which ParaView has been built. Leave the output with e and navigate to the line with **ParaView\_DIR**, hit Return and point to that directory. A second Return will leave the input. Press c to configure the build process. When no further error occurs, generate the Makefile with g.

When everything has been compiled successfully, you can install the plugin. For that purpose I made screenshots from the necessary steps. You can find them under *paraView\_plugin/install/*. There is also a screenshot of how to use it.

## Library

To include vtkbool into your program, you have to compile it as a library. The approach is the same as in the chapter before, except that you don't need ParaView. All you need is an installation of VTK with header files. If you have installed VTK over your package manager, CMake is able to find the required files. Otherwise you have to set **VTK\_DIR** manually. Make sure to disable **BUILD\_PARAVIEW\_PLUGIN** in CMake.

The usage of the library is very simple. Look at *testing.cxx* and you can see how. Upon creating the instance of the boolean-filter and connecting the two inputs with the pipeline, you can choose between different operation types, in different manners. You can set the operation mode by calling one of the named methods:

- SetOperModeToUnion
- SetOperModeToIntersection
- SetOperModeToDifference
- SetOperModeToDifference2

The alternative is the more generic `SetOperMode`. The method must be called with the number of the desired operation, an integer between 0 and 3, with the same meaning as mentioned before. After updating the pipeline, the result is stored in the first output, typically accessable with `GetOutputPort()`. The second output, `GetOutputPort(1)`, contains the lines of contact between the inputs. The inputs must be outputs of filters or sources returning vtkPolyData. The outputs from this filter are of the same type.

## Python

The Python module will be generated automatically, if three conditions are met:

- vtkbool is built as a library
- Python 2.x is installed with header files
- VTK itself is wrapped to Python

After a successful compilation, the module can be used as follows:

```python
import sys
sys.path.append('/path/to/your/build/directory')

import vtk
import libvtkboolPython as vtkbool

cube = vtk.vtkCubeSource()

sphere = vtk.vtkSphereSource()
sphere.SetCenter(.5, .5, .5)
sphere.SetThetaResolution(20)
sphere.SetPhiResolution(20)

boolean = vtkbool.vtkPolyDataBooleanFilter()
boolean.SetInputConnection(0, cube.GetOutputPort())
boolean.SetInputConnection(1, sphere.GetOutputPort())
boolean.SetOperModeToDifference()

# write the result, if you want ...

writer = vtk.vtkPolyDataWriter()
writer.SetInputConnection(boolean.GetOutputPort())
writer.SetFileName('result.vtk')

writer.Update()
```

## Examples

For demonstration there are examples for the use in ParaView. They can be found in *paraView\_plugin/demos*. To open a pvsm-file, select **Load State...** in the file-menu of ParaView.

![](https://raw.github.com/zippy84/vtkbool/master/examples.png)

## Copyright

2012-2016 Ronald Römer

## License

Apache License, Version 2.0
