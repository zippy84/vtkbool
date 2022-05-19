#!/usr/bin/env python
# *-* coding: UTF-8 *-*

# Copyright 2012-2022 Ronald Römer
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

from collections import defaultdict
from operator import itemgetter
import math
import pytest

from vtkmodules.vtkCommonCore import vtkIdList, vtkIdTypeArray, vtkPoints
from vtkmodules.vtkCommonDataModel import vtkPolyData, VTK_POLYGON
from vtkmodules.vtkFiltersSources import vtkCubeSource, vtkCylinderSource, vtkSphereSource
from vtkmodules.vtkCommonDataModel import vtkKdTreePointLocator
from vtkmodules.vtkFiltersCore import vtkAppendPolyData
from vtkmodules.vtkIOLegacy import vtkPolyDataWriter
from vtkmodules.vtkCommonExecutionModel import vtkTrivialProducer

from vtkBool import vtkPolyDataBooleanFilter, vtkPolyDataContactFilter

def find_pts(loc, pt, pts):
    loc.FindPointsWithinRadius(1e-5, pt, pts)

def get_next(id_, poly):
    return 0 if id_+1 == poly.GetNumberOfIds() else id_+1

def check_result(bf):
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

    for pd in [pdA, pdB]:
        print('checking pd')

        loc = vtkKdTreePointLocator()
        loc.SetDataSet(pd)
        loc.BuildLocator()

        it = lines.GetLines().NewIterator()

        while not it.IsDoneWithTraversal():

            line = it.GetCurrentCell()

            print('line', it.GetCurrentCellId())

            idA = line.GetId(0)
            idB = line.GetId(1)

            linkedA = vtkIdList()
            linkedB = vtkIdList()

            lines.GetPointCells(idA, linkedA)
            lines.GetPointCells(idB, linkedB)

            _linkedA = linkedA.GetNumberOfIds()
            _linkedB = linkedB.GetNumberOfIds()

            pA = lines.GetPoint(idA)
            pB = lines.GetPoint(idB)

            # print(pA)
            # print(pB)

            ptsA = vtkIdList()
            ptsB = vtkIdList()

            find_pts(loc, pA, ptsA)
            find_pts(loc, pB, ptsB)

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
                                    c = get_next(a, poly)
                                    assert c != b

                        group_a, group_b = groups_

                        possible_edges = [ [i, j] for i in group_a for j in group_b ]

                        edges = []

                        for a, b in possible_edges:
                            id_i, i = a
                            id_j, j = b

                            next_i = get_next(i, poly)
                            next_j = get_next(j, poly)

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
    check_result(bf)

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
    check_result(bf)

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
    check_result(bf)

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
    check_result(bf)

def test_same(tmp_path):
    sphere = vtkSphereSource()

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, sphere.GetOutputPort())
    bf.SetInputConnection(1, sphere.GetOutputPort())
    bf.SetOperModeToNone()

    bf.Update()

    write_result(bf, tmp_path)
    check_result(bf)

# def test_disjoin(tmp_path):
#     pts = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1], [.5, .5, 0], [.5, .5, 1]]
#     polys = [
#         [0, 1, 5, 4],
#         [1, 2, 6, 5],
#         [2, 3, 7, 6],
#         [3, 0, 4, 7],
#         [0, 8, 1],
#         [1, 8, 2],
#         [2, 8, 3],
#         [3, 8, 0],
#         [4, 5, 9],
#         [5, 6, 9],
#         [6, 7, 9],
#         [7, 4, 9]
#     ]

#     _pts = vtkPoints()
#     [ _pts.InsertNextPoint(*pt) for pt in pts ]

#     pd = vtkPolyData()
#     pd.Allocate(1)
#     pd.SetPoints(_pts)

#     for poly in polys:
#         cell = vtkIdList()
#         [ cell.InsertNextId(i) for i in poly ]
#         pd.InsertNextCell(VTK_POLYGON, cell)

#     cube = vtkCubeSource()
#     cube.SetCenter(1, .5, .5)

#     bf = vtkPolyDataBooleanFilter()
#     bf.SetInputData(0, pd)
#     bf.SetInputConnection(1, cube.GetOutputPort())
#     bf.SetOperModeToNone()

#     bf.Update()

#     write_result(bf, tmp_path)
#     check_result(bf)

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
    check_result(bf)

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
    # check_result(bf)

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
    # check_result(bf)
