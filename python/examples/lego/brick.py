#!/usr/bin/env python
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

# export LD_LIBRARY_PATH=/home/zippy/VTK8/lib

import sys
sys.path.extend(['/home/zippy/VTK8/lib/python3.6/site-packages',
    '/home/zippy/vtkbool/build'])

import vtk
import vtkboolPython

import math

def create_cylinder(r=.5, h=1., res=30):
    pts = vtk.vtkPoints()
    pts.SetNumberOfPoints(10*res)

    pd = vtk.vtkPolyData()
    pd.Allocate(3*res, 1)

    ang = 2*math.pi/res

    ind = 0

    for sign in [-1, 1]:
        y = sign*h/2
        for i in range(res):
            c0, s0 = math.cos(i*ang), math.sin(i*ang)
            c1, s1 = math.cos((i+1)*ang), math.sin((i+1)*ang)

            x0, z0 = round(r*c0, 6), round(r*s0, 6)
            x1, z1 = round(r*c1, 6), round(r*s1, 6)

            pts.InsertPoint(ind, 0, y, 0)
            pts.InsertPoint(ind+1, x0, y, z0)
            pts.InsertPoint(ind+2, x1, y, z1)

            tri = vtk.vtkIdList()
            [ tri.InsertNextId(ind+j) for j in range(3) ]

            t = pd.InsertNextCell(vtk.VTK_TRIANGLE, tri)

            ind += 3

            if sign == 1:
                pd.ReverseCell(t)

                pts.InsertPoint(ind, x0, -h/2, z0)
                pts.InsertPoint(ind+1, x0, y, z0)

                pts.InsertPoint(ind+2, x1, y, z1)
                pts.InsertPoint(ind+3, x1, -h/2, z1)

                poly = vtk.vtkIdList()
                [ poly.InsertNextId(ind+j) for j in range(4) ]

                pd.InsertNextCell(vtk.VTK_POLYGON, poly)

                ind += 4

    pd.SetPoints(pts)

    prod = vtk.vtkTrivialProducer()
    prod.SetOutput(pd)

    return prod

def add_tube(prev, x, y):
    cylA = create_cylinder(r=3.155, h=8.6)

    cube = vtk.vtkCubeSource()
    cube.SetXLength(.75)
    cube.SetYLength(6.6)
    cube.SetZLength(13.6)
    cube.SetCenter(0, 1, 0)

    bA = vtkboolPython.vtkPolyDataBooleanFilter()
    bA.SetInputConnection(cylA.GetOutputPort())
    bA.SetInputConnection(1, cube.GetOutputPort())
    bA.SetOperModeToUnion()
    bA.DecPolysOff()

    cylB = create_cylinder(r=2.405, h=8.6)

    bB = vtkboolPython.vtkPolyDataBooleanFilter()
    bB.SetInputConnection(bA.GetOutputPort())
    bB.SetInputConnection(1, cylB.GetOutputPort())
    bB.SetOperModeToDifference()
    bB.DecPolysOff()

    tr = vtk.vtkTransform()
    tr.Translate(x, y, 4.30001)
    tr.RotateX(90)

    tf = vtk.vtkTransformPolyDataFilter()
    tf.SetTransform(tr)
    tf.SetInputConnection(bB.GetOutputPort())

    bC = vtkboolPython.vtkPolyDataBooleanFilter()
    bC.SetInputConnection(prev.GetOutputPort())
    bC.SetInputConnection(1, tf.GetOutputPort())
    bC.SetOperModeToUnion()
    bC.DecPolysOff()

    return bC

def add_stud(prev, x, y):
    cylA = create_cylinder(r=2.5, h=2.)

    trA = vtk.vtkTransform()
    trA.Translate(x, y, 10.6)
    trA.RotateX(90)

    tfA = vtk.vtkTransformPolyDataFilter()
    tfA.SetTransform(trA)
    tfA.SetInputConnection(cylA.GetOutputPort())

    bA = vtkboolPython.vtkPolyDataBooleanFilter()
    bA.SetInputConnection(prev.GetOutputPort())
    bA.SetInputConnection(1, tfA.GetOutputPort())
    bA.SetOperModeToUnion()
    bA.DecPolysOff()

    cylB = create_cylinder(r=1.5, h=1.)

    trB = vtk.vtkTransform()
    trB.Translate(x, y, 9.1)
    trB.RotateX(90)

    tfB = vtk.vtkTransformPolyDataFilter()
    tfB.SetTransform(trB)
    tfB.SetInputConnection(cylB.GetOutputPort())

    bB = vtkboolPython.vtkPolyDataBooleanFilter()
    bB.SetInputConnection(bA.GetOutputPort())
    bB.SetInputConnection(1, tfB.GetOutputPort())
    bB.SetOperModeToDifference()
    bB.DecPolysOff()

    return bB

def add_bound(prev, x, y, phi):
    cube = vtk.vtkCubeSource()
    cube.SetXLength(.75)
    cube.SetYLength(.3)
    cube.SetZLength(8.6)

    tr = vtk.vtkTransform()
    tr.Translate(x, y, 4.3)
    tr.Push()
    tr.RotateZ(phi)
    tr.Translate(0, 2.65, 0)

    tf = vtk.vtkTransformPolyDataFilter()
    tf.SetTransform(tr)
    tf.SetInputConnection(cube.GetOutputPort())

    b = vtkboolPython.vtkPolyDataBooleanFilter()
    b.SetInputConnection(prev.GetOutputPort())
    b.SetInputConnection(1, tf.GetOutputPort())
    b.SetOperModeToUnion()
    b.DecPolysOff()

    return b


cubeA = vtk.vtkCubeSource()
cubeA.SetBounds(-16, 16, -8, 8, 0, 9.6)

triA = vtk.vtkTriangleFilter()
triA.SetInputConnection(cubeA.GetOutputPort())

subA = vtk.vtkLinearSubdivisionFilter()
subA.SetInputConnection(triA.GetOutputPort())
subA.SetNumberOfSubdivisions(4)

cubeB = vtk.vtkCubeSource()
cubeB.SetBounds(-14.8, 14.8, -6.8, 6.8, 0, 8.6)

triB = vtk.vtkTriangleFilter()
triB.SetInputConnection(cubeB.GetOutputPort())

subB = vtk.vtkLinearSubdivisionFilter()
subB.SetInputConnection(triB.GetOutputPort())
subB.SetNumberOfSubdivisions(4)

boolA = vtkboolPython.vtkPolyDataBooleanFilter()
boolA.SetInputConnection(subA.GetOutputPort())
boolA.SetInputConnection(1, subB.GetOutputPort())
boolA.SetOperModeToDifference()
boolA.DecPolysOff()

tubeA = add_tube(boolA, 0, 0)
tubeB = add_tube(tubeA, -8, 0)
tubeC = add_tube(tubeB, 8, 0)

studA = add_stud(tubeC, -12, -4)
studB = add_stud(studA, -12, 4)
studC = add_stud(studB, -4, -4)
studD = add_stud(studC, -4, 4)
studE = add_stud(studD, 4, -4)
studF = add_stud(studE, 4, 4)
studG = add_stud(studF, 12, -4)
studH = add_stud(studG, 12, 4)

boundA = add_bound(studH, -12, -4, 180)
boundB = add_bound(boundA, -12, 4, 0)
boundC = add_bound(boundB, -4, -4, 180)
boundD = add_bound(boundC, -4, 4, 0)
boundE = add_bound(boundD, 4, -4, 180)
boundF = add_bound(boundE, 4, 4, 0)
boundG = add_bound(boundF, 12, -4, 180)
boundH = add_bound(boundG, 12, 4, 0)

boundI = add_bound(boundH, -12, -4, 90)
boundJ = add_bound(boundI, -12, 4, 90)

boundK = add_bound(boundJ, 12, -4, -90)
boundL = add_bound(boundK, 12, 4, -90)

clean = vtk.vtkCleanPolyData()
clean.SetInputConnection(boundL.GetOutputPort())

writer = vtk.vtkPolyDataWriter()
writer.SetFileName('brick.vtk')
writer.SetInputConnection(clean.GetOutputPort())
writer.Update()
