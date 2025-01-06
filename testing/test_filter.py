#!/usr/bin/env python
# *-* coding: UTF-8 *-*

# Copyright 2012-2025 Ronald Römer
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
sys.path.extend(['/home/zippy/vtkbool/build/lib/python3.11/site-packages/vtkbool'])

from collections import defaultdict
from operator import itemgetter
import math
import pytest
from functools import reduce

from vtkmodules.vtkCommonCore import vtkIdList, vtkIdTypeArray, vtkPoints
from vtkmodules.vtkCommonDataModel import vtkPolyData, VTK_POLYGON
from vtkmodules.vtkFiltersSources import vtkCubeSource, vtkCylinderSource, vtkSphereSource, vtkPolyLineSource
from vtkmodules.vtkCommonDataModel import vtkKdTreePointLocator
from vtkmodules.vtkFiltersCore import vtkAppendPolyData, vtkCleanPolyData
from vtkmodules.vtkIOLegacy import vtkPolyDataWriter, vtkPolyDataReader
from vtkmodules.vtkCommonExecutionModel import vtkTrivialProducer
from vtkmodules.vtkFiltersModeling import vtkRotationalExtrusionFilter
from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter
from vtkmodules.vtkCommonTransforms import vtkTransform
from vtkmodules.vtkCommonExecutionModel import vtkAlgorithm

from vtkBool import vtkPolyDataBooleanFilter, vtkPolyDataContactFilter

from generate_frieze import extrude

def check_result(bf, expected_regs=None):
    lines = bf.GetOutput(2)

    lines.BuildLinks()

    pdA = bf.GetOutput(0)
    pdB = bf.GetOutput(1)

    pdA.BuildLinks()
    pdB.BuildLinks()

    contsA = lines.GetCellData().GetArray('cA')
    contsB = lines.GetCellData().GetArray('cB')

    assert isinstance(contsA, vtkIdTypeArray)
    assert isinstance(contsB, vtkIdTypeArray)

    regionsA = pdA.GetPointData().GetArray('RegionId')
    regionsB = pdB.GetPointData().GetArray('RegionId')

    assert isinstance(regionsA, vtkIdTypeArray)
    assert isinstance(regionsB, vtkIdTypeArray)

    if isinstance(expected_regs, list):
        valuesA = set(regionsA.GetValue(i) for i in range(regionsA.GetNumberOfValues()))
        valuesB = set(regionsB.GetValue(i) for i in range(regionsB.GetNumberOfValues()))

        assert len(valuesA) == expected_regs[0]
        assert len(valuesB) == expected_regs[1]

    for name, pd in [('pdA', pdA), ('pdB', pdB)]:
        print(f'checking {name}')

        loc = vtkKdTreePointLocator()
        loc.SetDataSet(pd)
        loc.BuildLocator()

        it = lines.GetLines().NewIterator()

        while not it.IsDoneWithTraversal():

            line = it.GetCurrentCell()

            # print('line', it.GetCurrentCellId())

            idA = line.GetId(0)
            idB = line.GetId(1)

            linkedA = vtkIdList()
            linkedB = vtkIdList()

            lines.GetPointCells(idA, linkedA)
            lines.GetPointCells(idB, linkedB)

            _linkedA = linkedA.GetNumberOfIds()
            _linkedB = linkedB.GetNumberOfIds()

            ptA = lines.GetPoint(idA)
            ptB = lines.GetPoint(idB)

            # print(ptA)
            # print(ptB)

            ptsA = vtkIdList()
            ptsB = vtkIdList()

            loc.FindPointsWithinRadius(1e-5, ptA, ptsA)
            loc.FindPointsWithinRadius(1e-5, ptB, ptsB)

            neigs = defaultdict(list)

            for pts in [ptsA, ptsB]:
                polys = vtkIdList()

                for i in range(pts.GetNumberOfIds()):
                    pd.GetPointCells(pts.GetId(i), polys)

                    for j in range(polys.GetNumberOfIds()):
                        neigs[polys.GetId(j)].append(pts.GetId(i))

            direct_neigs = {}

            for poly_id, point_ids in neigs.items():
                if len(point_ids) > 1:
                    equal_points = defaultdict(list)

                    for point_id in point_ids:
                        x, y, z = pd.GetPoint(point_id)

                        equal_points[f'{x:.5f},{y:.5f},{z:.5f}'].append(point_id)

                    groups = list(equal_points.values())

                    if len(groups) == 2:
                        poly = vtkIdList()
                        pd.GetCellPoints(poly_id, poly)

                        next_inds = [ i+1 for i in range(poly.GetNumberOfIds()-1) ] + [0]

                        groups_ = [ [ (id_, poly.IsId(id_)) for id_ in group ] for group in groups ]

                        for group in groups_:
                            if len(group) > 1:
                                group.sort(key=itemgetter(1))

                                # mehr als 2 punkte an gleicher stelle?
                                assert len(group) < 3

                                it_a = map(itemgetter(1), group + group[:1])
                                it_b = iter(it_a)
                                next(it_b)

                                for a, b in zip(it_a, it_b):
                                    # benachbart?
                                    c = next_inds[a]
                                    assert c != b

                        group_a, group_b = groups_

                        possible_edges = [ [i, j] for i in group_a for j in group_b ]

                        edges = []

                        for a, b in possible_edges:
                            id_i, i = a
                            id_j, j = b

                            next_i = next_inds[i]
                            next_j = next_inds[j]

                            if next_i == j \
                                or next_j == i:

                                # reihenfolge ist dann nicht mehr wichtig
                                edges.append([id_i, id_j])

                        # print(edges)

                        direct_neigs[poly_id] = edges

            # print(direct_neigs)

            match len(direct_neigs):
                case 2:
                    a, b = direct_neigs.values()

                    assert len(a) == 1
                    assert len(b) == 1

                    ids = set(a[0] + b[0])

                    if _linkedA == _linkedB:
                        assert len(ids) == 4

                    elif _linkedA == 2:
                        end_ids = set()

                        for i in range(_linkedA):
                            _line = vtkIdList()
                            lines.GetCellPoints(linkedA.GetId(i), _line)

                            end_ids.add(_line.GetId(0))
                            end_ids.add(_line.GetId(1))

                        end_ids.remove(idA)

                        if len(end_ids) == 2:
                            assert len(ids) == 4
                        else:
                            assert len(ids) == 3

                    elif _linkedB == 2:
                        end_ids = set()

                        for i in range(_linkedB):
                            _line = vtkIdList()
                            lines.GetCellPoints(linkedB.GetId(i), _line)

                            end_ids.add(_line.GetId(0))
                            end_ids.add(_line.GetId(1))

                        end_ids.remove(idB)

                        if len(end_ids) == 2:
                            assert len(ids) == 4
                        else:
                            assert len(ids) == 3

                    else:
                        # bspw. bei _linkedA == 4, _linkedB == 6
                        assert len(ids) == 4

                case 1:
                    a, *_ = direct_neigs.values()

                    assert len(a) == 2

                    ids = set(sum(a, []))

                    if _linkedA == 2 \
                        or _linkedB == 2:

                        assert len(ids) == 3

                    else:
                        assert len(ids) == 4

                case 0:
                    assert False

            it.GoToNextCell()

        print('-> ok')

def write_result(bf, d):
    for i, name in enumerate(['pdA', 'pdB', 'lines']):
        writer = vtkPolyDataWriter()
        writer.SetFileName(d / f'{name}.vtk')
        writer.SetInputConnection(bf.GetOutputPort(i))
        writer.Update()

def test_simple(tmp_path):
    cube = vtkCubeSource()
    cube.SetYLength(.5)

    cyl = vtkCylinderSource()
    cyl.SetResolution(32)
    cyl.SetHeight(.5)
    cyl.SetCenter(0, .5, 0)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, cyl.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [2, 2])

def test_simple_2(tmp_path):
    cube = vtkCubeSource()
    cube.SetYLength(.5)

    cyl = vtkCylinderSource()
    cyl.SetResolution(32)
    cyl.SetHeight(.5)
    cyl.SetCenter(0, .25, 0)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, cyl.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [2, 2])

def test_simple_3(tmp_path):
    cube = vtkCubeSource()
    cube.SetYLength(.5)

    cyl = vtkCylinderSource()
    cyl.SetResolution(32)
    cyl.SetHeight(.5)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, cyl.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [6, 6])

def test_simple_4(tmp_path):
    cube = vtkCubeSource()

    cyl = vtkPolyData()
    cyl.Allocate(1)

    pts = vtkPoints()

    bottom = vtkIdList()
    top = vtkIdList()

    phi = math.pi/16

    for i in range(32):
        x0 = .5*math.cos(i*phi)
        z0 = .5*math.sin(i*phi)

        x1 = .5*math.cos((i+1)*phi)
        z1 = .5*math.sin((i+1)*phi)

        top.InsertNextId(pts.InsertNextPoint(x0, .75, z0))
        bottom.InsertNextId(pts.InsertNextPoint(x0, -.25, z0))

        for j in range(4):
            side = vtkIdList()

            side.InsertNextId(pts.InsertNextPoint(x0, -.25+j/4, z0))
            side.InsertNextId(pts.InsertNextPoint(x0, -.25+(j+1)/4, z0))
            side.InsertNextId(pts.InsertNextPoint(x1, -.25+(j+1)/4, z1))
            side.InsertNextId(pts.InsertNextPoint(x1, -.25+j/4, z1))

            cyl.InsertNextCell(VTK_POLYGON, side)

    cyl.SetPoints(pts)

    cyl.ReverseCell(cyl.InsertNextCell(VTK_POLYGON, top))
    cyl.InsertNextCell(VTK_POLYGON, bottom)

    prod = vtkTrivialProducer()
    prod.SetOutput(cyl)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, prod.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [2, 2])

@pytest.mark.xfail
def test_simple_5():
    cube = vtkCubeSource()

    cyl = vtkCylinderSource()
    cyl.SetResolution(32)
    cyl.SetHeight(.75)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, cyl.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    check_result(bf)

def test_same(tmp_path):
    sphere = vtkSphereSource()

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, sphere.GetOutputPort())
    bf.SetInputConnection(1, sphere.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [56, 56])

@pytest.mark.xfail
def test_strips():
    cubeA = vtkCubeSource()

    cubeB = vtkCubeSource()
    cubeB.SetCenter(.3, .3, 0)

    cubeC = vtkCubeSource()
    cubeC.SetCenter(.15, 1.15, 0)
    cubeC.SetXLength(2)
    cubeC.SetYLength(2)

    app = vtkAppendPolyData()
    app.AddInputConnection(cubeA.GetOutputPort())
    app.AddInputConnection(cubeB.GetOutputPort())

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, app.GetOutputPort())
    bf.SetInputConnection(1, cubeC.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    check_result(bf)

def test_touch(tmp_path):
    cube = vtkCubeSource()

    cylinder = vtkCylinderSource()
    cylinder.SetResolution(32)
    cylinder.SetRadius(.25)
    cylinder.SetCenter(.25, 0, 0)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, cylinder.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [3, 3])

def test_nearly_congruent(tmp_path):
    phi = math.radians(.0005)

    z = math.tan(phi)*.5

    ptsA = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]
    polysA = [[0, 1, 2, 3]]

    ptsB = [[0, 0, 0], [0, .5, 0], [1, .5, 0], [1, 0, 0], [0, 1, z], [1, 1, z]]
    polysB = [
        [0, 1, 2, 3],
        [2, 1, 4, 5]
    ]

    _ptsA = vtkPoints()
    [ _ptsA.InsertNextPoint(*pt) for pt in ptsA ]

    _ptsB = vtkPoints()
    [ _ptsB.InsertNextPoint(*pt) for pt in ptsB ]

    pdA = vtkPolyData()
    pdA.Allocate(1)
    pdA.SetPoints(_ptsA)

    for poly in polysA:
        cell = vtkIdList()
        [ cell.InsertNextId(i) for i in poly ]
        pdA.InsertNextCell(VTK_POLYGON, cell)

    pdB = vtkPolyData()
    pdB.Allocate(1)
    pdB.SetPoints(_ptsB)

    for poly in polysB:
        cell = vtkIdList()
        [ cell.InsertNextId(i) for i in poly ]
        pdB.InsertNextCell(VTK_POLYGON, cell)

    writerA = vtkPolyDataWriter()
    writerA.SetInputData(pdA)
    writerA.SetFileName(tmp_path / 'pdA.vtk')
    writerA.Update()

    writerB = vtkPolyDataWriter()
    writerB.SetInputData(pdB)
    writerB.SetFileName(tmp_path / 'pdB.vtk')
    writerB.Update()

    cf = vtkPolyDataContactFilter()
    cf.SetInputData(0, pdA)
    cf.SetInputData(1, pdB)

    cf.Update()

    lines = cf.GetOutput()

    assert lines.GetNumberOfCells() == 1

    it = lines.GetLines().NewIterator()
    line = it.GetCellAtId(0)

    a = lines.GetPoint(line.GetId(0))
    b = lines.GetPoint(line.GetId(1))

    print(a, b)

    writerC = vtkPolyDataWriter()
    writerC.SetInputData(lines)
    writerC.SetFileName(tmp_path / 'lines.vtk')
    writerC.Update()

def test_merger(tmp_path):
    cube = vtkCubeSource()

    cyl = vtkCylinderSource()
    cyl.SetRadius(.25)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, cyl.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [3, 3])

def test_merger_2(tmp_path):
    cube = vtkCubeSource()

    pts = [
        [-.4, 0],
        [0, -.4],
        [.4, 0],
        [0, .4],
        [-.25, 0],
        [-.15, -.1],
        [-.05, 0],
        [-.15, .1],
        [.05, 0],
        [.15, -.1],
        [.25, 0],
        [.15, .1]
    ]

    polys = [
        [0, 2, 4, 20, 18, 16, 12, 10, 8],
        [4, 6, 0, 8, 14, 12, 16, 22, 20],
        [5, 3, 1, 9, 11, 13, 17, 19, 21],
        [1, 7, 5, 21, 23, 17, 13, 15, 9],
        [0, 1, 3, 2],
        [2, 3, 5, 4],
        [4, 5, 7, 6],
        [6, 7, 1, 0],
        [8, 10, 11, 9],
        [10, 12, 13, 11],
        [12, 14, 15, 13],
        [14, 8, 9, 15],
        [16, 18, 19, 17],
        [18, 20, 21, 19],
        [20, 22, 23, 21],
        [22, 16, 17, 23]
    ]

    _pts = vtkPoints()

    for pt in pts:
        _pts.InsertNextPoint(*pt, .5)
        _pts.InsertNextPoint(*pt, -.5)

    pd = vtkPolyData()
    pd.Allocate(1)
    pd.SetPoints(_pts)

    for poly in polys:
        cell = vtkIdList()
        [ cell.InsertNextId(i) for i in poly ]
        pd.InsertNextCell(VTK_POLYGON, cell)

    prod = vtkTrivialProducer()
    prod.SetOutput(pd)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, prod.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [7, 5])

def test_quads(tmp_path):
    sphereA = vtkSphereSource()
    sphereA.LatLongTessellationOn()
    sphereA.SetCenter(-.25, 0, 0)

    sphereB = vtkSphereSource()
    sphereB.LatLongTessellationOn()
    sphereB.SetCenter(.25, 0, 0)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, sphereA.GetOutputPort())
    bf.SetInputConnection(1, sphereB.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [2, 2])

def test_triangle_strips(tmp_path):
    line = vtkPolyLineSource()
    line.SetNumberOfPoints(19)

    phi = math.pi/18

    for i in range(19):
        line.SetPoint(i, 0, math.cos(i*phi), math.sin(i*phi))

    extr = vtkRotationalExtrusionFilter()
    extr.SetInputConnection(line.GetOutputPort())
    extr.SetResolution(18)
    extr.SetRotationAxis(0, 1, 0)

    sphere = vtkSphereSource()
    sphere.SetRadius(1)

    tra = vtkTransform()
    tra.RotateX(90)

    tf = vtkTransformPolyDataFilter()
    tf.SetInputConnection(sphere.GetOutputPort())
    tf.SetTransform(tra)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, extr.GetOutputPort())
    bf.SetInputConnection(1, tf.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [49, 49])

def test_special(tmp_path):
    cubeA = vtkCubeSource()
    cubeA.SetBounds(-5, 5, -10, 0, 0, 10)

    pts = [
        [0, 0],
        [1, 0],
        [2, -1],
        [3, 0],
        [4, -1],
        [5, 0],
        [5, 10],
        [-5, 10],
        [-5, 0],
        [-4, 1],
        [-3, 0],
        [-2, 1],
        [-1, 0]
    ]

    cubeB = extrude(pts, 10)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cubeA.GetOutputPort())
    bf.SetInputConnection(1, cubeB.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    write_result(bf, tmp_path)
    # check_result(bf, [10, 8])

    assert bf.GetOutput(2).GetNumberOfCells() == 32

@pytest.mark.xfail
def test_non_manifolds():
    reader = vtkPolyDataReader()
    reader.SetFileName('data/test.vtk')

    cyl = vtkCylinderSource()
    cyl.SetCenter(.75, 0, 0)
    cyl.SetRadius(.75)
    cyl.SetResolution(18)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, reader.GetOutputPort())
    bf.SetInputConnection(1, cyl.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    check_result(bf)

def test_branched(tmp_path):
    reader = vtkPolyDataReader()
    reader.SetFileName('data/branched.vtk')

    cube = vtkCubeSource()
    cube.SetBounds(-.3283653157, .3283653157, -.5, .5, -.375, .375)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, reader.GetOutputPort())
    bf.SetOperModeToNone()
    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [4, 4])

def test_branched_2(tmp_path):
    reader = vtkPolyDataReader()
    reader.SetFileName('data/branched.vtk')

    cube = vtkCubeSource()
    cube.SetBounds(-.3283653157, .3283653157, -.5, .5, -.75, .75)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, reader.GetOutputPort())
    bf.SetOperModeToNone()
    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [8, 8])

@pytest.mark.xfail
def test_branched_3():
    reader = vtkPolyDataReader()
    reader.SetFileName('data/branched3.vtk')

    cube = vtkCubeSource()
    cube.SetBounds(-.3283653157, .3283653157, -.5, .5, -.75, .75)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, reader.GetOutputPort())
    bf.SetOperModeToNone()
    bf.Update()

    check_result(bf)

def test_branched_4(tmp_path):
    reader = vtkPolyDataReader()
    reader.SetFileName('data/branched4.vtk')

    cube = vtkCubeSource()
    cube.SetBounds(-.3283653157, .3283653157, -.5, .5, -.375, .375)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, reader.GetOutputPort())
    bf.SetOperModeToNone()
    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [6, 7])

def test_branched_5(tmp_path):
    reader = vtkPolyDataReader()
    reader.SetFileName('data/branched.vtk')

    traA = vtkTransform()
    traA.Translate(0, -.5, 0)

    tfA = vtkTransformPolyDataFilter()
    tfA.SetInputConnection(reader.GetOutputPort())
    tfA.SetTransform(traA)

    traB = vtkTransform()
    traB.Translate(0, .5, 0)

    tfB = vtkTransformPolyDataFilter()
    tfB.SetInputConnection(reader.GetOutputPort())
    tfB.SetTransform(traB)

    app = vtkAppendPolyData()
    app.AddInputConnection(tfA.GetOutputPort())
    app.AddInputConnection(tfB.GetOutputPort())

    cube = vtkCubeSource()
    cube.SetBounds(-.3283653157, .3283653157, -1, 1, -.375, .375)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, app.GetOutputPort())
    bf.SetOperModeToNone()
    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [7, 8])

def test_branched_6(tmp_path):
    reader = vtkPolyDataReader()
    reader.SetFileName('data/branched6.vtk')

    cube = vtkCubeSource()
    cube.SetBounds(-.3283653157, .3283653157, -1, 1, -.33371031284, .33371031284)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, reader.GetOutputPort())
    bf.SetOperModeToNone()
    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [6, 5])

@pytest.mark.xfail
def test_bad_shaped():
    cube = vtkCubeSource()
    cube.SetBounds(-2.5, 2.5, 0, 5, 0, 5)

    z = 4.999

    pts = [
        [2.5, -2.5, 0],
        [2.5, -2.5, 5],
        [0, -2.5, z],
        [-2.5, -2.5, 5],
        [-2.5, -2.5, 0],

        [-2.5, 2.5, 0],
        [-2.5, 2.5, 5],
        [0, 2.5, z],
        [2.5, 2.5, 5],
        [2.5, 2.5, 0],
    ]

    polys = [
        [0, 1, 2, 3, 4],
        [5, 6, 7, 8, 9],
        [9, 8, 1, 0],
        [4, 3, 6, 5],
        [9, 0, 4, 5],
        [1, 8, 7, 6, 3, 2]
    ]

    _pts = vtkPoints()

    for pt in pts:
        _pts.InsertNextPoint(*pt)

    pd = vtkPolyData()
    pd.Allocate(1)
    pd.SetPoints(_pts)

    for poly in polys:
        cell = vtkIdList()
        [ cell.InsertNextId(i) for i in poly ]
        pd.InsertNextCell(VTK_POLYGON, cell)

    prod = vtkTrivialProducer()
    prod.SetOutput(pd)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cube.GetOutputPort())
    bf.SetInputConnection(1, prod.GetOutputPort())
    bf.SetOperModeToNone()
    bf.Update()

    check_result(bf)

@pytest.mark.xfail
def test_self_intersecting_polys():
    cube = vtkCubeSource()
    cube.SetCenter(0, 0, 1.375) # 1, 1.25, 1.375

    pts = [
        [.5, .5, 0],
        [-.5, .5, 0],
        [-.5, -.5, 0],
        [.5 , -.5, 0],

        [.5, 0, .5], # 4
        [0, .5, .5],
        [-.5, 0, .5],
        [0, -.5, .5],

        [.5, 0, .75], # 8
        [0, .5, .75],
        [-.5, 0, .75],
        [0, -.5, .75],

        [.5, 0, 1], # 12
        [.5, .5, 1],
        [0, .5, 1],
        [-.5, .5, 1],
        [-.5, 0, 1], # 16
        [-.5, -.5, 1],
        [0, -.5, 1],
        [.5, -.5, 1]
    ]

    polys = [
        [3, 2, 1, 0],
        [12, 13, 14, 15, 16, 17, 18, 19],
        [3, 19, 18, 11, 7, 11, 18, 17, 2],
        [0, 13, 12, 8, 4, 8, 12, 19, 3],
        [1, 15, 14, 9, 5, 9, 14, 13, 0],
        [2, 17, 16, 10, 6, 10, 16, 15, 1]
    ]

    _pts = vtkPoints()

    for pt in pts:
        _pts.InsertNextPoint(*pt)

    pd = vtkPolyData()
    pd.Allocate(1)
    pd.SetPoints(_pts)

    for poly in polys:
        cell = vtkIdList()
        [ cell.InsertNextId(i) for i in poly ]
        pd.InsertNextCell(VTK_POLYGON, cell)

    prod = vtkTrivialProducer()
    prod.SetOutput(pd)

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, prod.GetOutputPort())
    bf.SetInputConnection(1, cube.GetOutputPort())
    bf.SetOperModeToNone()
    bf.Update()

    check_result(bf)

def test_equal_capt_pts(tmp_path):
    cyl = vtkCylinderSource()
    cyl.SetHeight(2)
    cyl.SetResolution(12)

    z = .0000025

    tra = vtkTransform()
    tra.RotateZ(90)
    tra.Translate(0, 0, z)

    tf = vtkTransformPolyDataFilter()
    tf.SetTransform(tra)
    tf.SetInputConnection(cyl.GetOutputPort())

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cyl.GetOutputPort())
    bf.SetInputConnection(1, tf.GetOutputPort())
    bf.SetOperModeToNone()
    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [4, 4])

def test_equal_capt_pts_2(tmp_path):
    reader = vtkPolyDataReader()
    reader.SetFileName('data/cross.vtk')

    z = .000001

    tra = vtkTransform()
    tra.RotateZ(45)
    tra.Translate(0, 0, z)

    tf = vtkTransformPolyDataFilter()
    tf.SetTransform(tra)
    tf.SetInputConnection(reader.GetOutputPort())

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, reader.GetOutputPort())
    bf.SetInputConnection(1, tf.GetOutputPort())
    bf.SetOperModeToNone()
    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [24, 24])

def test_equal_capt_pts_3(tmp_path):

    # from vtkmodules.util.execution_model import select_ports

    # cyl = vtkCylinderSource()
    # cyl.SetHeight(2)
    # cyl.SetResolution(12)

    # cyls = []

    # for x in [-.5, -.25, 0, .25, .5]:
    #     tra = vtkTransform()
    #     tra.Translate(x, 0, 0)

    #     cyls.append(cyl >> vtkTransformPolyDataFilter(transform=tra))

    # def fct(a, b):
    #     bf = vtkPolyDataBooleanFilter()

    #     a >> bf
    #     b >> select_ports(1, bf)

    #     return bf

    # union = reduce(fct, cyls)

    # (union >> vtkCleanPolyData(output_points_precision=vtkAlgorithm.DOUBLE_PRECISION) >> vtkPolyDataWriter(file_name='data/complex.vtk')).update()

    cyl = vtkCylinderSource()
    cyl.SetHeight(2)
    cyl.SetResolution(12)

    reader = vtkPolyDataReader()
    reader.SetFileName('data/complex.vtk')

    z = .0000025

    tra = vtkTransform()
    tra.RotateZ(90)
    tra.Translate(0, 0, z)

    tf = vtkTransformPolyDataFilter()
    tf.SetTransform(tra)
    tf.SetInputConnection(reader.GetOutputPort())

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, cyl.GetOutputPort())
    bf.SetInputConnection(1, tf.GetOutputPort())
    bf.SetOperModeToNone()
    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf, [14, 14])
