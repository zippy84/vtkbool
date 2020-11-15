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

# export LD_LIBRARY_PATH=/home/zippy/VTK8/lib

import sys
sys.path.extend(['/home/zippy/VTK8/lib/python3.6/site-packages',
    '/home/zippy/vtkbool/build'])

import vtkboolPython
import vtk

import os
import re

from collections import namedtuple

Bnds = namedtuple('Bounds', 'x1 x2 y1 y2 z1 z2')

class Alignment:
    def __init__(self, first):
        self.aligned = [first]

    def add_right(self, flt, ref):
        flt.Update()
        ref.Update()

        bnds_a = Bnds(*flt.GetOutput().GetBounds())
        bnds_b = Bnds(*ref.GetOutput().GetBounds())

        tr = vtk.vtkTransform()
        tr.Translate(bnds_b.x2-bnds_a.x1-1, bnds_b.y1-bnds_a.y1, 0)

        tp = vtk.vtkTransformPolyDataFilter()
        tp.SetTransform(tr)
        tp.SetInputConnection(flt.GetOutputPort())

        try:
            self.aligned.remove(flt)
        except ValueError:
            pass

        self.aligned.append(tp)

        return tp

    def add_top(self, flt, ref):
        flt.Update()
        ref.Update()

        bnds_a = Bnds(*flt.GetOutput().GetBounds())
        bnds_b = Bnds(*ref.GetOutput().GetBounds())

        tr = vtk.vtkTransform()
        tr.Translate(bnds_b.x1-bnds_a.x1, bnds_b.y2-bnds_a.y1-1, 0)

        tp = vtk.vtkTransformPolyDataFilter()
        tp.SetTransform(tr)
        tp.SetInputConnection(flt.GetOutputPort())

        try:
            self.aligned.remove(flt)
        except ValueError:
            pass

        self.aligned.append(tp)

        return tp

    def write(self, *names):
        flts = [self.aligned[0]]

        for i, f in enumerate(self.aligned[1:]):
            bf = vtkboolPython.vtkPolyDataBooleanFilter()
            bf.SetInputConnection(flts[-1].GetOutputPort())
            bf.SetInputConnection(1, f.GetOutputPort())

            flts.append(bf)

        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(flts[-1].GetOutputPort())

        for name in names:
            if name.endswith('stl'):
                w = vtk.vtkSTLWriter()
            else:
                w = vtk.vtkPolyDataWriter()

            w.SetInputConnection(clean.GetOutputPort())
            w.SetFileName(name)
            w.Update()

def extrude(pts, h, z=0):
    cell = vtk.vtkIdList()

    _pts = vtk.vtkPoints()
    _pts.SetDataTypeToDouble()

    [ (_pts.InsertNextPoint(pt[0], pt[1], z), cell.InsertNextId(i)) for i, pt in enumerate(pts) ]

    pd = vtk.vtkPolyData()
    pd.Allocate(1, 1)
    pd.SetPoints(_pts)
    pd.InsertNextCell(vtk.VTK_POLYGON, cell)

    prod = vtk.vtkTrivialProducer()
    prod.SetOutput(pd)

    extr = vtk.vtkLinearExtrusionFilter()
    extr.SetInputConnection(prod.GetOutputPort())
    extr.SetVector(0, 0, h)

    pn = vtk.vtkPolyDataNormals()
    pn.SetInputConnection(extr.GetOutputPort())
    pn.AutoOrientNormalsOn()

    return pn

def add_frame(flt, holds = []):
    flt.Update()
    bnds = Bnds(*flt.GetOutput().GetBounds())

    pts1 = [(bnds.x1-2, bnds.y1-2, 0), (bnds.x2+2, bnds.y1-2, 0), (bnds.x2+2, bnds.y2+2, 0), (bnds.x1-2, bnds.y2+2, 0)]
    pts2 = [(bnds.x1-1, bnds.y1-1, 0), (bnds.x2+1, bnds.y1-1, 0), (bnds.x2+1, bnds.y2+1, 0), (bnds.x1-1, bnds.y2+1, 0)]

    extr1 = extrude(pts1, 1)
    extr2 = extrude(pts2, 1)

    bf = vtkboolPython.vtkPolyDataBooleanFilter()
    bf.SetInputConnection(extr1.GetOutputPort())
    bf.SetInputConnection(1, extr2.GetOutputPort())
    bf.SetOperModeToDifference()

    app = vtk.vtkAppendPolyData()
    app.AddInputConnection(flt.GetOutputPort())
    app.AddInputConnection(bf.GetOutputPort())

    conf = { 'left': ((.5, 0), (0, -90), lambda x, y: abs(bnds.x1-x)),
        'right': ((-.5, 0), (0, 90), lambda x, y: abs(bnds.x2-x)),
        'top': ((0, -.5), (-90, 0), lambda x, y: abs(bnds.y2-y)),
        'bottom': ((0, .5), (90, 0), lambda x, y: abs(bnds.y1-y)) }

    if len(holds) > 0:
        app2 = vtk.vtkAppendPolyData()

        for x, y, dir_ in holds:
            mv, rot, fct = conf[dir_]
            pts = [(.5+mv[0], -.5+mv[1], 0), (.5+mv[0], .5+mv[1], 0), (-.5+mv[0], .5+mv[1], 0), (-.5+mv[0], -.5+mv[1], 0)]

            extr = extrude(pts, fct(x, y)+1)

            tr = vtk.vtkTransform()
            tr.Translate(x, y, 0)
            tr.RotateX(rot[0])
            tr.RotateY(rot[1])

            tp = vtk.vtkTransformPolyDataFilter()
            tp.SetTransform(tr)
            tp.SetInputConnection(extr.GetOutputPort())

            app2.AddInputConnection(tp.GetOutputPort())

        bf2 = vtkboolPython.vtkPolyDataBooleanFilter()
        bf2.SetInputConnection(app.GetOutputPort())
        bf2.SetInputConnection(1, app2.GetOutputPort())

        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(bf2.GetOutputPort())

        return clean

        # app2.AddInputConnection(app.GetOutputPort())
        # return app2

    return app

def repeat(flt, count):
    flt.Update()

    bnds = Bnds(*flt.GetOutput().GetBounds())

    app = vtk.vtkAppendPolyData()

    for i in range(count):
        tr = vtk.vtkTransform()
        tr.Translate(i*(bnds.x2-bnds.x1+1), 0, 0)

        tp = vtk.vtkTransformPolyDataFilter()
        tp.SetTransform(tr)
        tp.SetInputConnection(flt.GetOutputPort())

        app.AddInputConnection(tp.GetOutputPort())

    return app

def read1():
    r = vtk.vtkSTLReader()
    r.SetFileName('Schuerze1.stl')
    fr = add_frame(r, [(-1.5, 10.25, 'left'), (-1.5, 28, 'left'), (-1.5, 45.75, 'left')])
    return fr

def read2():
    r = vtk.vtkSTLReader()
    r.SetFileName('Schuerze2.stl')
    fr = add_frame(r, [(-1.5, 23.1913, 'left'), (-1.5, 43.7588, 'left')])
    return fr

def read3():
    r = vtk.vtkSTLReader()
    r.SetFileName('Schuerze3.stl')
    fr = add_frame(r, [(-1.5, -11.2838, 'left'), (-1.5, -31.8512, 'left')])
    return fr

def split4(_=False):
    r = vtk.vtkSTLReader()
    r.SetFileName('Schuerze4.stl')
    r.Update()

    bnds = Bnds(*r.GetOutput().GetBounds())

    y = bnds.y1/2

    pts = [[bnds.x1-1, y, bnds.z1-1], [bnds.x2+1, y, bnds.z1-1], [bnds.x2+1, y, bnds.z2+1], [bnds.x1-1, y, bnds.z2+1]]

    if _:
        pts.reverse()

    cell = vtk.vtkIdList()

    vtk_pts = vtk.vtkPoints()
    vtk_pts.SetDataTypeToDouble()

    [ (vtk_pts.InsertNextPoint(pt), cell.InsertNextId(i)) for i, pt in enumerate(pts) ]

    pd = vtk.vtkPolyData()
    pd.Allocate(1, 1)
    pd.SetPoints(vtk_pts)
    pd.InsertNextCell(vtk.VTK_POLYGON, cell)

    prod = vtk.vtkTrivialProducer()
    prod.SetOutput(pd)

    bf = vtkboolPython.vtkPolyDataBooleanFilter()
    bf.SetInputConnection(r.GetOutputPort())
    bf.SetInputConnection(1, prod.GetOutputPort())
    bf.SetOperModeToDifference()

    clean = vtk.vtkCleanPolyData()
    clean.SetInputConnection(bf.GetOutputPort())

    pn = vtk.vtkPolyDataNormals()
    pn.SetInputConnection(clean.GetOutputPort())
    pn.AutoOrientNormalsOn()

    return pn

def read4():
    s = split4()
    fr = add_frame(s, [(-1.5, -80.0285, 'left'), (-1.5, -102.6081, 'left'), (-1.5, -125.1877, 'left')])
    return fr

def read5():
    s = split4(True)
    fr = add_frame(s, [(-1.5, -12.2898, 'left'), (-1.5, -34.8694, 'left'), (-1.5, -57.449, 'left')])
    return fr

def create_sill(l):
    pr1 = [(0, 0), (1.6666, 0), (1.6666, 1.875), (0, .91666)]
    pr2 = [(.6666, 0), (1.6666, 0), (1.6666, 1.875), (.6666, 1.875)]

    extr = extrude(pr1, l+2*.8333)

    extr2 = extrude(pr2, .8333)
    extr3 = extrude(pr2, .8333, l+.8333)

    bf1 = vtkboolPython.vtkPolyDataBooleanFilter()
    bf1.SetInputConnection(extr.GetOutputPort())
    bf1.SetInputConnection(1, extr2.GetOutputPort())

    bf2 = vtkboolPython.vtkPolyDataBooleanFilter()
    bf2.SetInputConnection(bf1.GetOutputPort())
    bf2.SetInputConnection(1, extr3.GetOutputPort())

    clean = vtk.vtkCleanPolyData()
    clean.SetInputConnection(bf2.GetOutputPort())

    pn = vtk.vtkPolyDataNormals()
    pn.SetInputConnection(clean.GetOutputPort())
    pn.AutoOrientNormalsOn()

    tr = vtk.vtkTransform()
    tr.Translate(1.6666, 0, 0)
    tr.RotateZ(180)
    tr.RotateX(90)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(pn.GetOutputPort())

    return tp

def create_sill2():
    pr1 = [(0, 0), (1.6666, 0), (1.6666, 1.875), (0, .91666)]
    pr2 = [(.6666, 0), (1.6666, 0), (1.6666, 1.875), (.6666, 1.875)]

    extr = extrude(pr1, 11.9075+2*.8333)

    extr2 = extrude(pr2, .8333)
    extr3 = extrude(pr2, .8333, 11.9075+.8333)
    extr4 = extrude(pr2, 3.25/3*2, 11.9075/2+.8333-3.25/3)

    bf1 = vtkboolPython.vtkPolyDataBooleanFilter()
    bf1.SetInputConnection(extr.GetOutputPort())
    bf1.SetInputConnection(1, extr2.GetOutputPort())

    bf2 = vtkboolPython.vtkPolyDataBooleanFilter()
    bf2.SetInputConnection(bf1.GetOutputPort())
    bf2.SetInputConnection(1, extr3.GetOutputPort())

    bf3 = vtkboolPython.vtkPolyDataBooleanFilter()
    bf3.SetInputConnection(bf2.GetOutputPort())
    bf3.SetInputConnection(1, extr4.GetOutputPort())

    clean = vtk.vtkCleanPolyData()
    clean.SetInputConnection(bf3.GetOutputPort())

    pn = vtk.vtkPolyDataNormals()
    pn.SetInputConnection(clean.GetOutputPort())
    pn.AutoOrientNormalsOn()

    tr = vtk.vtkTransform()
    tr.Translate(1.6666, 0, 0)
    tr.RotateZ(180)
    tr.RotateX(90)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(pn.GetOutputPort())

    return tp

def create_sill3():
    pr1 = [(0, 0), (1.75, 0), (1.75, 5./6*1.875), (0, 1.0833)]
    pr2 = [(.75, 0), (1.75, 0), (1.75, 5./6*1.875), (.75, 5./6*1.875)]

    extr = extrude(pr1, 22.7325+2*.54125)

    extr2 = extrude(pr2, .54125)
    extr3 = extrude(pr2, .54125, 22.7325+.54125)

    bf1 = vtkboolPython.vtkPolyDataBooleanFilter()
    bf1.SetInputConnection(extr.GetOutputPort())
    bf1.SetInputConnection(1, extr2.GetOutputPort())

    bf2 = vtkboolPython.vtkPolyDataBooleanFilter()
    bf2.SetInputConnection(bf1.GetOutputPort())
    bf2.SetInputConnection(1, extr3.GetOutputPort())

    clean = vtk.vtkCleanPolyData()
    clean.SetInputConnection(bf2.GetOutputPort())

    pn = vtk.vtkPolyDataNormals()
    pn.SetInputConnection(clean.GetOutputPort())
    pn.AutoOrientNormalsOn()

    tr = vtk.vtkTransform()
    tr.Translate(1.75, 0, 0)
    tr.RotateZ(180)
    tr.RotateX(90)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(pn.GetOutputPort())

    return tp

def create_sill4(l):
    pr1 = [(0, 0), (1.75, 0), (1.75, 5./6*1.875), (0, 1.0833)]
    pr2 = [(.75, 0), (1.75, 0), (1.75, 5./6*1.875), (.75, 5./6*1.875)]

    extr1 = extrude(pr1, l+.75+.54125)
    extr2 = extrude(pr2, .54125)

    bf = vtkboolPython.vtkPolyDataBooleanFilter()
    bf.SetInputConnection(extr1.GetOutputPort())
    bf.SetInputConnection(1, extr2.GetOutputPort())
    bf.Update()

    bnds = Bnds(*bf.GetOutput().GetBounds())

    pts = [(bnds.x1-2, bnds.y1-2, 0), (bnds.x2+2, bnds.y1-2, 0), (bnds.x2+2, bnds.y2+2, 0), (bnds.x1-2, bnds.y2+2, 0)]

    cell = vtk.vtkIdList()

    _pts = vtk.vtkPoints()
    _pts.SetDataTypeToDouble()

    [ (_pts.InsertNextPoint(pt), cell.InsertNextId(i)) for i, pt in enumerate(pts[::-1]) ]

    pd = vtk.vtkPolyData()
    pd.Allocate(1, 1)
    pd.SetPoints(_pts)
    pd.InsertNextCell(vtk.VTK_POLYGON, cell)

    prod = vtk.vtkTrivialProducer()
    prod.SetOutput(pd)

    tr = vtk.vtkTransform()
    tr.Translate(0, 0, bnds.z2)
    tr.RotateY(45)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(prod.GetOutputPort())

    bf2 = vtkboolPython.vtkPolyDataBooleanFilter()
    bf2.SetInputConnection(bf.GetOutputPort())
    bf2.SetInputConnection(1, tp.GetOutputPort())
    bf2.SetOperModeToDifference()

    clean = vtk.vtkCleanPolyData()
    clean.SetInputConnection(bf2.GetOutputPort())

    pn = vtk.vtkPolyDataNormals()
    pn.SetInputConnection(clean.GetOutputPort())
    pn.AutoOrientNormalsOn()

    tr2 = vtk.vtkTransform()
    tr2.Translate(1.75, 0, 0)
    tr2.RotateZ(180)
    tr2.RotateX(90)

    tp2 = vtk.vtkTransformPolyDataFilter()
    tp2.SetTransform(tr2)
    tp2.SetInputConnection(pn.GetOutputPort())

    return tp2

def mirror(flt):
    tr = vtk.vtkTransform()
    tr.Scale(-1, 1, 1)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(flt.GetOutputPort())

    rs = vtk.vtkReverseSense()
    rs.SetInputConnection(tp.GetOutputPort())

    return rs

def rotate(flt):
    tr = vtk.vtkTransform()
    tr.RotateZ(90)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(flt.GetOutputPort())

    return tp

if __name__ == '__main__':
    stl1 = read1()
    stl2 = read2()
    stl3 = read3()
    stl4 = read4()
    stl5 = read5()

    A = Alignment(stl1)

    a2 = A.add_right(stl2, stl1)
    a3 = A.add_right(stl3, a2)
    a4 = A.add_right(stl4, a3)
    a5 = A.add_right(stl5, a4)

    holds1 = sum([ [(.5+i*2.6666, 0, 'bottom'), (.5+i*2.6666, 11.9075+2*.8333, 'top')] for i in range(4) ], [])
    holds2 = sum([ [(.5+i*2.6666, 0, 'bottom'), (.5+i*2.6666, 9.7425+2*.8333, 'top')] for i in range(7) ], [])
    holds3 = sum([ [(.5+i*2.6666, 0, 'bottom'), (.5+i*2.6666, 8.66+2*.8333, 'top')] for i in range(4) ], [])
    holds4 = sum([ [(.5+i*2.6666, 0, 'bottom'), (.5+i*2.6666, 11.9075+2*.8333, 'top')] for i in range(2) ], [])

    a6 = A.add_right(rotate(add_frame(repeat(create_sill(11.9075), 4), holds1)), a5)
    a7 = A.add_top(rotate(add_frame(repeat(create_sill(9.7425), 7), holds2)), a6)
    a8 = A.add_top(rotate(add_frame(repeat(create_sill(8.66), 4), holds3)), a7)
    a9 = A.add_top(rotate(add_frame(repeat(create_sill2(), 2), holds4)), a8)

    # anbau

    holds5 = [(.5, 0, 'bottom'), (.5, 22.7325+2*.54125, 'top')]
    holds6= [(.5, 0, 'bottom')]

    p1 = add_frame(create_sill3(), holds5)

    a10 = A.add_right(p1, a6)

    p2 = add_frame(create_sill4(5.95375+1.0825), holds6)
    p3 = add_frame(create_sill4(11.9075+1.0825), holds6)

    a11 = A.add_right(p2, a10)
    a12 = A.add_top(mirror(p3), a11)

    a13 = A.add_right(mirror(p2), a11)
    a14 = A.add_top(p3, a13)

    A.write('all.vtk', 'all.stl')

    os.makedirs('frames', exist_ok=True)

    for k, v in dict(locals()).items():
        if k == 'stl1' or re.match(r'a\d+', k):
            w = vtk.vtkPolyDataWriter()
            w.SetInputConnection(v.GetOutputPort())
            w.SetFileName(f'frames/{k}.vtk')
            w.Update()
