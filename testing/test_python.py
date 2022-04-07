#!/usr/bin/env python
# *-* coding: UTF-8 *-*

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

# export LD_LIBRARY_PATH=/home/zippy/VTK9/lib

import sys
sys.path.extend(['/home/zippy/VTK9/lib/python3.10/site-packages',
    '/home/zippy/vtkbool/build/lib/python3.10/site-packages/vtkbool'])

import vtk
import vtkBool

cubeA = vtk.vtkCubeSource()
cubeB = vtk.vtkCubeSource()

bf = vtkBool.vtkPolyDataBooleanFilter()
bf.SetInputConnection(0, cubeA.GetOutputPort())
bf.SetInputConnection(1, cubeB.GetOutputPort())

w = vtk.vtkPolyDataWriter()
w.SetInputConnection(bf.GetOutputPort(1))
w.SetFileName('test.vtk')
w.Update()
