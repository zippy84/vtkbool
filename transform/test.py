#!/usr/bin/env python
# *-* coding: UTF-8 *-*

import sys
sys.path.extend(['/home/zippy/vtkbool/build/lib/python3.13/site-packages/vtkbool'])

from vtkmodules.vtkIOLegacy import vtkPolyDataReader, vtkPolyDataWriter
from vtkmodules.vtkCommonMath import vtkMatrix4x4
from vtkmodules.vtkCommonTransforms import vtkTransform
from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter

from vtkbool import vtkPolyDataBooleanFilter

mandibleReader = vtkPolyDataReader()
mandibleReader.SetFileName('bone.vtk')
mandibleReader.Update()

boxReader = vtkPolyDataReader()
boxReader.SetFileName('box.vtk')
boxReader.Update()

moveBoxTransform = vtkTransform()
moveBoxTransform.Translate(20.076142063093318, -8.288617300972092, -53.29349952031101)
moveBoxMatrix = moveBoxTransform.GetMatrix()

bf = vtkPolyDataBooleanFilter()
bf.SetInputData(0, mandibleReader.GetOutput())
bf.SetInputData(1, boxReader.GetOutput())

# identityMatrix = vtkMatrix4x4()
# bf.SetMatrix(0, identityMatrix)

for i in range(10):
    moveBoxTransform.Translate(0, -10, 0)
    moveBoxMatrix = moveBoxTransform.GetMatrix()

    print(moveBoxMatrix.GetElement(1, 3))

    # transformFilter = vtkTransformPolyDataFilter()
    # transformFilter.SetInputData(boxReader.GetOutput())
    # transformFilter.SetTransform(moveBoxTransform)
    # transformFilter.Update()

    # _writer = vtkPolyDataWriter()
    # _writer.SetFileName(f'out/pdB_{i}.vtk')
    # _writer.SetInputConnection(transformFilter.GetOutputPort())
    # _writer.Update()

    bf.SetMatrix(1, moveBoxMatrix)

    writer = vtkPolyDataWriter()
    writer.SetFileName(f'out/transformed_{i}.vtk')
    writer.SetInputConnection(bf.GetOutputPort())

    writer.Update()
