# vtkbool

![](https://raw.github.com/zippy84/vtkbool/master/cover.png)

## About

This is an extension of the graphics library VTK. The goal of the extension is to equip the library with boolean operations on polygonal meshes. I started the project at the end of my studies in mechanical engineering at the University of Applied Sciences ([HTWK](http://htwk-leipzig.de/)) in Leipzig. I used VTK to develop a program, which I had to create for a paper. At this time I would have wished, that this feature already exists. There was several implementations from third parties, but after some tests, I came to the conclusion, that none of them worked correct. The once who worked was bound to other big libraries like CGAL. I decided to start with my own implementation. That is the result of my efforts.

## Features

- no extra libraries required
- triangulation is not needed
- concave polygons are allowed
- only the involved polygons are modified
- meshes can be stacked (coplanar polygons are right handled)
- 4 operation types (union, intersection, difference and difference2 - difference with interchanged operands)
- all types of polygonal cells are supported (triangles, quads, polygons, triangle-strips)
- the meshes don’t need to be closed
- CellData is passed (attached by the rules of vtkAppendPolyData)
- original cell ids are added to CellData (*OrigCellIdsA* and *OrigCellIdsB* as vtkIntArray)
- contact lines are available
- it is also a plugin for Paraview
- Python wrapped
- ready for VTK 6

## Limitations

- the filter produces concave polygons (the visualization can be broken)
- PointData is not preserved - you have to do your own mapping
- correctness depends on the polygon-defining normals (use vtkPolyDataNormals if you have problems with incorrect orientations)
- when polygons are cut inside, holes are formed - the filter is not able to embed  holes

## Todo

- integration of holes
- acceleration with OpenMP

## Build Requirements

CMake >= 2.8, VTK >= 6.1, Paraview >= 5.0, gcc >= 4.7 (including libstdc++)

Other versions may also work. It is hard to find out, since when a feature is supported.

## Paraview plugin

To build the plugin you have to compile Paraview from source. Download the current version from <http://www.paraview.org> and follow the compilation instructions. As soon as Paraview is compiled, it may take a while, you can build the plugin with the following commands:

	git clone git://github.com/zippy84/vtkbool.git vtkbool
	cd vtkbool
	mkdir build
	cd build
	ccmake ..
	make

Once you have called CMake and typed c for the first time, you will get an error. The error says, that CMake could not found Paraview. You have to specify the directory in which Paraview has been built. Leave the output with e and navigate to the line with **ParaView\_DIR**, hit Return and point to that directory. A second Return will leave the input. Press c to configure the build process. When no further error occurs, generate the Makefile with g.

When everything has been compiled successfully, you can install the plugin. For that purpose I made screenshots from the necessary steps. You can find them under *paraview\_plugin/install/*. There is also a screenshot of how to use it.

## Library

To include vtkbool into your program, you have to compile it as a library. The approach is the same as in the chapter before, except that you don't need Paraview. All you need is an installation of VTK with header files. If you have installed VTK over your package manager, CMake is able to find the required files. Otherwise you have to set **VTK\_DIR** manually. Make sure to disable **BUILD\_PARAVIEW\_PLUGIN** in cmake.

The usage of the library is very simple. Look at *tests.cxx* and you can see how. Upon creating the instance of the boolean-filter and connecting the two inputs with the pipeline, you can choose between different operation types, in different manners. You can set the operation mode by calling one of the named methods:

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

For the purpose of demonstration and testing, I created examples in Paraview. The examples can be found in *paraview\_plugin/tests*. To open a contained pvsm-file, select Load State... in the file-menu of Paraview. So far, the examples contain all possible circumstances that can appear, when two meshes are combined. Tell me when this is not true!

![](https://raw.github.com/zippy84/vtkbool/master/examples.png)

## Copyright

2012-2016 Ronald Römer

## License

Apache License, Version 2.0
