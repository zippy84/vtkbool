#!/usr/bin/env python
# *-* coding: UTF-8 *-*

# Copyright 2012-2023 Ronald Römer
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

import sys
sys.path.extend(['/home/zippy/vtkbool/build/lib/python3.10/site-packages/vtkbool'])

from vtkBool import vtkPolyDataBooleanFilter

from vtkmodules.vtkFiltersSources import vtkSphereSource, vtkLineSource
from vtkmodules.vtkFiltersCore import vtkTubeFilter
from vtkmodules.vtkIOLegacy import vtkPolyDataWriter
from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter
from vtkmodules.vtkCommonTransforms import vtkTransform

from collections import defaultdict

sphere = vtkSphereSource()
sphere.SetRadius(5)
sphere.SetCenter(.5, .5, .5)
sphere.SetPhiResolution(100)
sphere.SetThetaResolution(100)
sphere.Update()

pd = sphere.GetOutput()

center = pd.GetCenter()

line = vtkLineSource()
line.SetPoint1(0, 0, 0)
line.SetPoint2(0, 0, 4)

tube = vtkTubeFilter()
tube.SetInputConnection(line.GetOutputPort())
tube.SetRadius(1)
tube.SetNumberOfSides(50)
tube.CappingOn()

transform = vtkTransform()
transform.PostMultiply()
transform.Translate(center[0], center[1], center[2]+3)

tf = vtkTransformPolyDataFilter()
tf.SetInputConnection(tube.GetOutputPort())
tf.SetTransform(transform)

moves = [(-.5, 0), (-.5, 0), (-.5, 0), (-.5, 0), (-.5, 0), (-.5, 0), (-.5, 0), (-.5, 0), (-.5, 0)]

# moves = [(.5, 0), (.5, 0), (.5, 0), (.5, 0),
#     (0, -.5), (0, -.5), (0, -.5), (0, -.5),
#     (-.5, 0), (-.5, 0), (-.5, 0), (-.5, 0),
#     (0, .5)]

# moves = [(0, -.5), (0, -.5), (0, -.5), (0, -.5), (0, -.5), (0, -.5), (0, -.5)]

for i, xy in enumerate(moves):
    transform.Translate(*xy, 0)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputData(0, pd)
    bf.SetInputConnection(1, tf.GetOutputPort())
    bf.SetOperModeToDifference()

    writer = vtkPolyDataWriter()
    writer.SetFileName(f'sphere{i}.vtk')
    writer.SetInputData(pd)
    writer.Update()

    writer2 = vtkPolyDataWriter()
    writer2.SetFileName(f'tube{i}.vtk')
    writer2.SetInputConnection(tf.GetOutputPort())
    writer2.Update()

    writer3 = vtkPolyDataWriter()
    writer3.SetFileName(f'bool{i}.vtk')
    writer3.SetInputConnection(bf.GetOutputPort())
    writer3.Update()

    pd.Initialize()
    pd.DeepCopy(bf.GetOutput())
