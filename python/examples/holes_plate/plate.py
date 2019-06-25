#!/usr/bin/env python2
# *-* coding: UTF-8 *-*

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

# export LD_LIBRARY_PATH=/home/zippy/VTK6/lib

import sys
sys.path.extend(['/home/zippy/VTK6/lib/python2.7/site-packages',
    '/home/zippy/vtkbool/build'])

import vtkboolPython
import vtk

import math

# n, m = 6, 6
n, m = 4, 4

app = vtk.vtkAppendPolyData()

cyl = vtk.vtkCylinderSource()
# cyl.SetResolution(4)
cyl.SetResolution(8)

for i in range(n):
    for j in range(m):
        tr = vtk.vtkTransform()
        tr.Translate(1.5*i, 0, 1.5*j)

        tf = vtk.vtkTransformPolyDataFilter()
        tf.SetInputConnection(cyl.GetOutputPort())
        tf.SetTransform(tr)

        app.AddInputConnection(tf.GetOutputPort())

box = vtk.vtkCubeSource()
box.SetXLength(n*1.5)
box.SetZLength(m*1.5)
box.SetCenter((n*1.5)/2-.75, 0, (m*1.5)/2-.75)

bf = vtkboolPython.vtkPolyDataBooleanFilter();
bf.SetInputConnection(0, box.GetOutputPort())
bf.SetInputConnection(1, app.GetOutputPort())
bf.SetOperModeToDifference()

w = vtk.vtkPolyDataWriter()
w.SetInputConnection(bf.GetOutputPort())
w.SetFileName('plate.vtk')
w.Update()
