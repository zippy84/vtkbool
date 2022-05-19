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
import math

from collections import namedtuple, defaultdict
from operator import attrgetter, itemgetter

from vtkmodules.vtkCommonCore import vtkIdList, vtkPoints
from vtkmodules.vtkCommonDataModel import vtkPolyData, VTK_POLYGON, VTK_LINE
from vtkmodules.vtkIOLegacy import vtkPolyDataWriter
from vtkmodules.vtkCommonDataModel import vtkKdTree
from vtkmodules.vtkFiltersFlowPaths import vtkModifiedBSPTree

Conn = namedtuple('Conn', 'd,i,j')

def write_vtk(pts, polys, name):
    _pts = vtkPoints()

    pd = vtkPolyData()
    pd.Allocate(1)

    for poly in polys:
        cell = vtkIdList()

        for i in poly:
            cell.InsertNextId(_pts.InsertNextPoint(pts[i]))

        pd.InsertNextCell(VTK_POLYGON, cell)

    pd.SetPoints(_pts)

    w = vtkPolyDataWriter()
    w.SetInputData(pd)
    w.SetFileName(f'{name}.vtk')
    w.Update()

def write_vtk_2(pd, name):
    w = vtkPolyDataWriter()
    w.SetInputData(pd)
    w.SetFileName(f'{name}.vtk')
    w.Update()

def draw(radius, step, x, y, rotate=0):
    pts = []

    phi = 2*math.pi/step

    for i in range(step):
        pts.append([radius*math.cos(i*phi), radius*math.sin(i*phi), 0])

    phi = math.radians(rotate)

    for pt in pts:
        pt[:2] = [math.cos(phi)*pt[0]-math.sin(phi)*pt[1], math.sin(phi)*pt[0]+math.cos(phi)*pt[1],]

        pt[0] += x
        pt[1] += y

    return pts

_t = 6

_polys = [
    draw(8, 18, 0, 0),
    draw(.25, _t, 3.5, 0),
    draw(.25, _t, -3.5, 0),

    draw(.5, _t, 1, 0),
    draw(.5, _t, -1, 0),
    draw(.5, _t, 0, 1, 30),
    draw(.5, _t, 0, -1, 30)
]

_pts = vtkPoints()

pts = sum(_polys, [])
polys = [ [ _pts.InsertNextPoint(p) for p in poly ] for poly in _polys ]

write_vtk(pts, polys, 'pd')

kd_tree = vtkKdTree()
kd_tree.OmitZPartitioning()
kd_tree.BuildLocatorFromPoints(_pts)

sources = []

lines = vtkPolyData()
lines.Allocate(1)
lines.SetPoints(_pts)

for i, poly in enumerate(polys):
    for p in poly:
        sources.append(i)

    _poly = poly[1:]
    _poly.append(poly[0])

    for a, b in zip(poly, _poly):
        line = vtkIdList()
        line.InsertNextId(a)
        line.InsertNextId(b)
        lines.InsertNextCell(VTK_LINE, line)

tree = vtkModifiedBSPTree()
tree.SetDataSet(lines)
tree.BuildLocator()

n = 0

conns_per_poly = defaultdict(list)

def find_conns(ids=[]):
    global n

    if n > len(pts):
        return False

    n += 5

    new_conns = defaultdict(list)

    for _id, poly in enumerate(polys):
        if ids \
            and _id not in ids:

            continue

        for i in poly:
            # print(i)

            found_pts = vtkIdList()
            kd_tree.FindClosestNPoints(n, pts[i], found_pts)

            for _j in range(found_pts.GetNumberOfIds()):
                j = found_pts.GetId(_j)

                src_a = sources[i]
                src_b = sources[j]

                pt_a = pts[i]
                pt_b = pts[j]

                if src_a != src_b:
                    # print(i, j)

                    _pts = vtkPoints()
                    _cells = vtkIdList()

                    good = True

                    if tree.IntersectWithLine(pt_a, pt_b, 1e-5, _pts, _cells) == 1:
                        for k in range(_cells.GetNumberOfIds()):
                            line = vtkIdList()
                            lines.GetCellPoints(_cells.GetId(k), line)

                            id0 = line.GetId(0)
                            id1 = line.GetId(1)

                            if id0 != i and id1 != i \
                                and id0 != j and id1 != j:
                                # print((i, j), '->', (id0, id1))

                                good = False

                                break

                    if good:
                        v = [pt_a[0]-pt_b[0], pt_a[1]-pt_b[1]]
                        d = v[0]*v[0]+v[1]*v[1]

                        # print((src_a, src_b))

                        new_conns[(src_a, src_b)].append(Conn(d, i, j))
                        new_conns[(src_b, src_a)].append(Conn(d, j, i))

    for i, conns in conns_per_poly.items():
        if not ids \
            or i in ids:

            conns.clear()

    for k, v in new_conns.items():
        a, b = k
        conns_per_poly[a].extend(v)

    for conns in conns_per_poly.values():
        conns.sort(key=attrgetter('d'))

    return True

find_conns()

connected = defaultdict(list, {0: []})

restricted = [] # keine der conns darf im gleichen punkt beginnen

lines_2 = vtkPolyData()
lines_2.Allocate(1)
lines_2.SetPoints(_pts)

tree_2 = vtkModifiedBSPTree()
tree_2.SetDataSet(lines_2)

while len(connected) < len(polys):
    found_one = False

    for i, poly in enumerate(polys):
        if i in connected:
            continue

        for conn in conns_per_poly[i]:
            if sources[conn.j] in connected \
                and conn.j not in restricted \
                and tree_2.IntersectWithLine(pts[conn.i], pts[conn.j], 1e-5, None, None) == 0:

                print(i, conn)
                connected[sources[conn.i]].append(conn)
                restricted.extend([conn.i, conn.j])

                # das andere poly auch aktualisieren
                connected[sources[conn.j]].append(Conn(conn.d, conn.j, conn.i))

                found_one = True

                line = vtkIdList()
                line.InsertNextId(conn.i)
                line.InsertNextId(conn.j)

                lines_2.InsertNextCell(VTK_LINE, line)

                tree_2.Modified()

                break

    if not found_one:
        if not find_conns():
            break

if len(connected) != len(polys):
    sys.exit('merging failed')

write_vtk_2(lines_2, 'lines_before')

print(connected)

chains = {}

for i in connected:
    chain = [i]

    while chain[-1] != 0:
        chain.append(sources[connected[chain[-1]][0].j])

    chains[i] = chain

print('chains', chains)

solved = set([0])

for i, poly in enumerate(polys):
    if i > 0 \
        and connected[i]:

        conns = connected[i]

        while len(conns) < 2:
            first = conns[0]

            prios = {}

            for conn in conns_per_poly[i]:
                src = sources[conn.j]

                if src in prios:
                    continue

                if conn.i != first.i \
                    and conn.j != first.j \
                    and conn.i not in restricted \
                    and conn.j not in restricted \
                    and tree_2.IntersectWithLine(pts[conn.i], pts[conn.j], 1e-5, None, None) == 0:

                    chain_a = chains[i]
                    chain_b = chains[src]

                    # gemeinsame eltern

                    shared = set(chain_a) & set(chain_b)

                    # print(chain_a, chain_b, '->', shared)

                    # gemeinsame eltern müssen sich alle in solved befinden

                    if shared <= solved:
                        new_chain_a = [ p for p in chain_a if p not in solved ]
                        new_chain_b = [ p for p in chain_b if p not in solved ]

                        solvable = set(new_chain_a + new_chain_b)

                        prios[src] = (conn, solvable, len(solvable), -conn.d)

            if prios:
                _prios = sorted(prios.values(), key=itemgetter(2, 3))

                print(f'\033[46mprios {_prios}\033[0m')

                conn = _prios[-1][0]
                conns.append(conn)

                connected[sources[conn.j]].append(Conn(conn.d, conn.j, conn.i))

                restricted.extend([conn.i, conn.j])

                line = vtkIdList()
                line.InsertNextId(conn.i)
                line.InsertNextId(conn.j)

                lines_2.InsertNextCell(VTK_LINE, line)

                tree_2.Modified()

                solved.update(_prios[-1][1])

                print(f'\033[41m{solved}\033[0m')

                break

            else:
                if not find_conns([i]):
                    sys.exit('merging failed')

print(connected)

write_vtk_2(lines_2, 'lines_after')

print(restricted)

# stage 1

used_conns = []

merged = polys[0]

for i, conns in connected.items():
    if i == 0:
        continue

    f = conns[0]

    print(f)

    assert f.j in merged

    index_a = merged.index(f.j)

    rotated_a = merged[index_a:] + merged[:index_a]
    rotated_a.append(merged[index_a])

    print(rotated_a)

    poly_b = polys[sources[f.i]]
    index_b = poly_b.index(f.i)

    rotated_b = poly_b[index_b:] + poly_b[:index_b]
    rotated_b.append(poly_b[index_b])

    print(rotated_b)

    merged = rotated_a + list(reversed(rotated_b))

    print(merged)

    used_conns.append(f)

write_vtk(pts, [merged], 'merged_stage1')

# stage 2

print(used_conns)

splitted = [merged]

left_conns = set()

for i, conns in connected.items():
    if i == 0:
        continue

    for conn in conns[1:]:
        if not any(u.i == conn.j and u.j == conn.i for u in used_conns):
            print(conn)

            if conn.i < conn.j:
                left_conns.add(Conn(0, conn.i, conn.j))
            else:
                left_conns.add(Conn(0, conn.j, conn.i))

print(left_conns)

for conn in left_conns:
    print(conn)
    for i, poly in enumerate(splitted):
        try:
            index_a = poly.index(conn.i) # kann nur in einem polygon existieren
        except ValueError:
            continue

        # einbauen

        rotated = poly[index_a:] + poly[:index_a]

        print(rotated)

        assert conn.j in rotated

        index_b = rotated.index(conn.j)

        print(index_a, index_b)

        poly_a = rotated[:index_b]
        poly_b = rotated[index_b:]

        poly_a.append(rotated[index_b])
        poly_b.append(rotated[0])

        print(poly_a)
        print(poly_b)

        splitted[i:i+1] = [poly_a, poly_b]

        break

splitted.extend(polys[1:])

write_vtk(pts, splitted, 'merged_stage2')

# probe

for poly in splitted:
    unique_pts = set()

    for i in poly:
        p = f'{pts[i][0]:.5f}, {pts[i][1]:.5f}'

        if p in unique_pts:
            sys.exit('merging failed')
        else:
            unique_pts.add(p)

# assert len(splitted) == 10
