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

import math
import os
import re

from teile import Alignment, extrude, repeat, add_frame

def merge1():
    r = vtk.vtkPolyDataReader()
    r.SetFileName('einzeln/test4.vtk')
    fr = add_frame(r, [(-9.858333, 2.5, 'top'), (-29.575, 2.5, 'top'), (-49.291666, 2.5, 'top'), (-69.008333, 2.5, 'top'), (-88.725, 2.5, 'top')])
    return fr

def merge2():
    r = vtk.vtkPolyDataReader()
    r.SetFileName('einzeln/test5.vtk')
    fr = add_frame(r, [(-9.208333, 2.5, 'top'), (-27.625, 2.5, 'top'), (-46.041666, 2.5, 'top'), (-64.458333, 2.5, 'top'), (-82.875, 2.5, 'top')])
    return fr

def merge3():
    r = vtk.vtkPolyDataReader()
    r.SetFileName('einzeln/test0.vtk')

    r2 = vtk.vtkPolyDataReader()
    r2.SetFileName('einzeln/test1.vtk')

    tr = vtk.vtkTransform()
    tr.RotateZ(90)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(r2.GetOutputPort())

    bf = vtkboolPython.vtkPolyDataBooleanFilter()
    bf.SetInputConnection(r.GetOutputPort())
    bf.SetInputConnection(1, tp.GetOutputPort())
    bf.DecPolysOff()

    fr = add_frame(bf, [(-10.8333, 2.5, 'top'), (-32.5, 2.5, 'top')])

    return fr

def merge4():
    r = vtk.vtkPolyDataReader()
    r.SetFileName('einzeln/test6.vtk')

    r2 = vtk.vtkPolyDataReader()
    r2.SetFileName('einzeln/test7.vtk')

    tr = vtk.vtkTransform()
    tr.RotateZ(270)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(r2.GetOutputPort())

    bf = vtkboolPython.vtkPolyDataBooleanFilter()
    bf.SetInputConnection(r.GetOutputPort())
    bf.SetInputConnection(1, tp.GetOutputPort())
    bf.DecPolysOff()

    fr = add_frame(bf, [(10.8333, 2.5, 'top'), (32.5, 2.5, 'top')])

    return fr

def merge5():
    r = vtk.vtkPolyDataReader()
    r.SetFileName('einzeln/test2.vtk')
    fr = add_frame(r, [(-10.8333, 2.5, 'top'), (-32.5, 2.5, 'top')])
    return fr

def merge6():
    r = vtk.vtkPolyDataReader()
    r.SetFileName('einzeln/test8.vtk')
    fr = add_frame(r, [(17.333, 2.5, 'top'), (34.666, 2.5, 'top')])
    return fr

def merge7():
    r = vtk.vtkPolyDataReader()
    r.SetFileName('einzeln/test14.vtk')
    fr = add_frame(r, [(6.2291666, 2.5, 'top')])
    return fr

def merge8():
    r = vtk.vtkPolyDataReader()
    r.SetFileName('einzeln/test12.vtk')

    r2 = vtk.vtkPolyDataReader()
    r2.SetFileName('einzeln/test13.vtk')

    tr = vtk.vtkTransform()
    tr.RotateZ(90)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(r.GetOutputPort())

    bf = vtkboolPython.vtkPolyDataBooleanFilter()
    bf.SetInputConnection(r2.GetOutputPort())
    bf.SetInputConnection(1, tp.GetOutputPort())
    bf.DecPolysOff()

    fr = add_frame(bf, [(-6.2291666, 2.5, 'top')])

    return fr

def merge9():
    r = vtk.vtkPolyDataReader()
    r.SetFileName('einzeln/test10.vtk')
    fr = add_frame(r, [(-8.486111, 2.5, 'top'), (-25.458333, 2.5, 'top'), (-42.430555, 2.5, 'top')])
    return fr

def merge10():
    r = vtk.vtkPolyDataReader()
    r.SetFileName('einzeln/test11.vtk')
    fr = add_frame(r, [(-9.1924, 2.5, 'top')])
    return fr

def merge11():
    r = vtk.vtkPolyDataReader()
    r.SetFileName('einzeln/test9.vtk')
    fr = add_frame(r, [(4.291666, 2.5, 'top')])
    return fr

def merge12():
    r = vtk.vtkPolyDataReader()
    r.SetFileName('einzeln/test3.vtk')
    fr = add_frame(r, [(5.6875, 2.5, 'top')])
    return fr

def create_rest():
    pr = [(0, 0), (1.5, 0), (1.5, 1.875/3), (1.75, 1.875/3), (1.75, 2*1.875/3), (2, 2*1.875/3), (2, 1.875), (0, 1.875)]
    extr = extrude(pr, 1.625)

    tr = vtk.vtkTransform()
    tr.Translate(1.875/2, 0, 0)
    tr.RotateZ(90)

    tp = vtk.vtkTransformPolyDataFilter()
    tp.SetTransform(tr)
    tp.SetInputConnection(extr.GetOutputPort())

    return tp

if __name__ == '__main__':
    p1 = merge1()
    p2 = merge2()
    p3 = merge3()
    p4 = merge4()
    p5 = merge5()
    p6 = merge6()
    p7 = merge7()
    p8 = merge8()
    p9 = merge9()
    p10 = merge10()
    p11 = merge11()
    p12 = merge12()

    A = Alignment(p1)

    a2 = A.add_top(p2, p1)
    a3 = A.add_top(p3, a2)
    a4 = A.add_right(p4, a3)
    a5 = A.add_top(p5, a3)
    a6 = A.add_right(p6, a5)
    a7 = A.add_top(p7, a5)
    a8 = A.add_right(p8, a7)
    a9 = A.add_right(p9, a8)
    a10 = A.add_top(p10, a7)
    a11 = A.add_right(p11, a10)
    a12 = A.add_right(p12, a11)

    holds = [ (i*2.875, 0, 'bottom') for i in range(12) ]

    a13 = A.add_top(add_frame(repeat(create_rest(), 12), holds), a10)

    A.write('band.vtk', 'band.stl')

    os.makedirs('frames2', exist_ok=True)

    for k, v in dict(locals()).items():
        if k == 'p1' or re.match(r'a\d+', k):
            w = vtk.vtkPolyDataWriter()
            w.SetInputConnection(v.GetOutputPort())
            w.SetFileName(f'frames2/{k}.vtk')
            w.Update()
