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
sys.path.extend(['/home/zippy/VTK8/lib/python3.6/site-packages'])

import vtk
import os

from frieze import extrude

def Merge1 ():
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName('test4.vtk')

    reader2 = vtk.vtkPolyDataReader()
    reader2.SetFileName('test14.vtk')

    app = vtk.vtkAppendPolyData()
    app.AddInputConnection(reader.GetOutputPort())

    tr = vtk.vtkTransform()
    tr.RotateZ(90)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(reader2.GetOutputPort())

    app.AddInputConnection(tp.GetOutputPort())

    return app

def Merge2 ():
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName('test5.vtk')
    reader.Update()

    return reader

def Merge3 ():
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName('test0.vtk')

    reader2 = vtk.vtkPolyDataReader()
    reader2.SetFileName('test1.vtk')

    reader3 = vtk.vtkPolyDataReader()
    reader3.SetFileName('test2.vtk')

    reader4 = vtk.vtkPolyDataReader()
    reader4.SetFileName('test3.vtk')

    app = vtk.vtkAppendPolyData()
    app.AddInputConnection(reader.GetOutputPort())

    tr = vtk.vtkTransform()
    tr.RotateZ(90)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(reader2.GetOutputPort())

    app.AddInputConnection(tp.GetOutputPort())

    tr2 = vtk.vtkTransform()
    tr2.Translate(43.333333, 3.25, 0)

    tp2 = vtk.vtkTransformPolyDataFilter()
    tp2.SetTransform(tr2)
    tp2.SetInputConnection(reader3.GetOutputPort())

    app.AddInputConnection(tp2.GetOutputPort())

    tr3 = vtk.vtkTransform()
    tr3.Translate(43.333333, 3.25, 0)
    tr3.RotateZ(90)

    tp3 = vtk.vtkTransformPolyDataFilter()
    tp3.SetTransform(tr3)
    tp3.SetInputConnection(reader4.GetOutputPort())

    app.AddInputConnection(tp3.GetOutputPort())

    return app

def Merge4 ():
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName('test6.vtk')

    reader2 = vtk.vtkPolyDataReader()
    reader2.SetFileName('test7.vtk')

    reader3 = vtk.vtkPolyDataReader()
    reader3.SetFileName('test8.vtk')

    app = vtk.vtkAppendPolyData()
    app.AddInputConnection(reader.GetOutputPort())

    tr = vtk.vtkTransform()
    tr.RotateZ(270)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(reader2.GetOutputPort())

    app.AddInputConnection(tp.GetOutputPort())

    tr2 = vtk.vtkTransform()
    tr2.Translate(-43.333333, 3.25, 0)

    tp2 = vtk.vtkTransformPolyDataFilter()
    tp2.SetTransform(tr2)
    tp2.SetInputConnection(reader3.GetOutputPort())

    app.AddInputConnection(tp2.GetOutputPort())

    tr3 = vtk.vtkTransform()
    tr3.RotateX(180)

    tp3 = vtk.vtkTransformPolyDataFilter()
    tp3.SetTransform(tr3)
    tp3.SetInputConnection(app.GetOutputPort())

    return tp3

def Merge5 ():
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName('test10.vtk')

    reader2 = vtk.vtkPolyDataReader()
    reader2.SetFileName('test12.vtk')

    reader3 = vtk.vtkPolyDataReader()
    reader3.SetFileName('test13.vtk')

    app = vtk.vtkAppendPolyData()

    app.AddInputConnection(reader.GetOutputPort())

    tr = vtk.vtkTransform()
    tr.Translate(-50.916666, -1.083333, 0)
    tr.RotateZ(90)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(reader2.GetOutputPort())

    app.AddInputConnection(tp.GetOutputPort())

    tr2 = vtk.vtkTransform()
    tr2.Translate(-50.916666, -1.083333, 0)

    tp2 = vtk.vtkTransformPolyDataFilter()
    tp2.SetTransform(tr2)
    tp2.SetInputConnection(reader3.GetOutputPort())

    app.AddInputConnection(tp2.GetOutputPort())

    return app

def Merge6 ():
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName('test11.vtk')

    reader2 = vtk.vtkPolyDataReader()
    reader2.SetFileName('test9.vtk')

    app = vtk.vtkAppendPolyData()
    app.AddInputConnection(reader.GetOutputPort())

    tr = vtk.vtkTransform()
    tr.RotateZ(45)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(reader2.GetOutputPort())

    app.AddInputConnection(tp.GetOutputPort())

    return app

s = 3.041666
t = 1

data = [{ 'obj': Merge1(), 'dy': 0 },
    { 'obj': Merge2(), 'dy': s+t },
    { 'obj': Merge3(), 'dy': 2*(s+t) },
    { 'obj': Merge4(), 'dy': 3*(s+t) },
    { 'obj': Merge5(), 'dy': 4*(s+t)+2.166666 },
    { 'obj': Merge6(), 'dy': 5*(s+t)+3.25 }]

for d in data:
    d['obj'].Update()
    d['ext'] = d['obj'].GetOutput().GetBounds()

merg = vtk.vtkAppendPolyData()

for d in data:
    tr = vtk.vtkTransform()
    tr.Translate(-d['ext'][0], -d['ext'][2]+d['dy'], 0)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(d['obj'].GetOutputPort())

    merg.AddInputConnection(tp.GetOutputPort())

def Bridge (x, z = 1.316666):
    app = vtk.vtkAppendPolyData()

    plane = vtk.vtkPlane()
    plane.SetOrigin(x, 0, 0)
    plane.SetNormal(1, 0, 0)

    cutter = vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputConnection(merg.GetOutputPort())

    cf = vtk.vtkConnectivityFilter()
    cf.SetInputConnection(cutter.GetOutputPort())
    cf.SetExtractionModeToAllRegions()

    cf.Update()

    regs = cf.GetNumberOfExtractedRegions()

    cf.SetExtractionModeToSpecifiedRegions()

    ms = []

    for i in range(regs):
        cf.InitializeSpecifiedRegionList()
        cf.AddSpecifiedRegion(i)

        gf = vtk.vtkGeometryFilter()
        gf.SetInputConnection(cf.GetOutputPort())

        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(gf.GetOutputPort())
        clean.Update()

        pd = clean.GetOutput()

        ys = []

        for j in range(pd.GetNumberOfPoints()):
            pt = pd.GetPoint(j)

            if pt[2] > z-1e-5:
                ys.append(pt[1])

        if len(ys) > 1:
            ys.sort()

            m = ys[-1]+(ys[0]-ys[-1])/2

            cyl = vtk.vtkCylinderSource()
            cyl.SetCenter(0, 0, m)
            cyl.SetResolution(12)
            cyl.SetRadius(.4)

            app.AddInputConnection(cyl.GetOutputPort())

            ms.append(m)

    ov = 2

    ex = extrude([(0, 0), (2, 0), (2, .75), (0, .75)], 0, ms[-1]-ms[0]+2*ov)

    tr = vtk.vtkTransform()
    tr.Translate(-1, -1.25, ms[0]-ov)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(ex.GetOutputPort())

    app.AddInputConnection(tp.GetOutputPort())

    tr2 = vtk.vtkTransform()
    tr2.Translate(x, 0, z+.5)
    tr2.RotateX(270)

    tp2 = vtk.vtkTransformPolyDataFilter()
    tp2.SetTransform(tr2)
    tp2.SetInputConnection(app.GetOutputPort())

    return tp2

merg2 = vtk.vtkAppendPolyData()
merg2.AddInputConnection(merg.GetOutputPort())

mx = max([ d['ext'][1]-d['ext'][0] for d in data ])/8

for i in range(8):
    br = Bridge(mx/2+i*mx)
    merg2.AddInputConnection(br.GetOutputPort())

writer = vtk.vtkPolyDataWriter()
writer.SetInputConnection(merg2.GetOutputPort())
writer.SetFileName('frieze.vtk')
writer.Update()

if os.path.exists('stl'):

    tri = vtk.vtkTriangleFilter()
    tri.SetInputConnection(merg2.GetOutputPort())

    stl = vtk.vtkSTLWriter()
    stl.SetInputConnection(tri.GetOutputPort())
    stl.SetFileName('stl/frieze.stl')
    stl.Update()
