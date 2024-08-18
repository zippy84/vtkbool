#!/usr/bin/env python
# *-* coding: UTF-8 *-*

import math
from collections import defaultdict, namedtuple
from operator import itemgetter

from vtkmodules.vtkCommonCore import vtkIdList, vtkIdTypeArray, vtkPoints, vtkMath, reference
from vtkmodules.vtkCommonDataModel import vtkPolyData, VTK_POLYGON, VTK_QUAD, VTK_VERTEX, VTK_TRIANGLE, vtkKdTreePointLocator
from vtkmodules.vtkFiltersSources import vtkCylinderSource
from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter
from vtkmodules.vtkCommonTransforms import vtkTransform
from vtkmodules.vtkIOLegacy import vtkPolyDataWriter, vtkPolyDataReader
from vtkmodules.vtkFiltersFlowPaths import vtkModifiedBSPTree
from vtkmodules.vtkFiltersCore import vtkCleanPolyData

import sys
sys.path.extend(['/home/zippy/vtkbool/build/lib/python3.12/site-packages/vtkbool'])

from vtkBool import vtkPolyDataBooleanFilter

from test import SnapPoint, SnapEdge, SnapPoly, write, compute_normal, Base, in_poly, _proj_line, proj_line, get_points, clean, Prepare

class Align:
    def __init__(self, mesh_a, mesh_b):
        self.mesh_a = clean(mesh_a)
        self.mesh_b = clean(mesh_b)

        write(self.mesh_a, 'cleaned_a_2.vtk')
        write(self.mesh_b, 'cleaned_b_2.vtk')

        self.congr_ids_a = set()
        self.congr_ids_b = set()

        self.align_congruent_points()

        deb_a = None
        deb_b = None

        # deb_a = 17696
        # deb_b = 29901

        self.lines_a = defaultdict(set)
        self.lines_b = defaultdict(set)

        self.find_inters(self.mesh_a, self.mesh_b, self.congr_ids_a, self.lines_a, self.lines_b, 'a', deb_a)
        self.find_inters(self.mesh_b, self.mesh_a, self.congr_ids_b, self.lines_b, self.lines_a, 'b', deb_b)

        self.lines_a = { edge: { ind: [sys.float_info.max] for ind in dat } for edge, dat in self.lines_a.items() }
        self.lines_b = { edge: { ind: [sys.float_info.max] for ind in dat } for edge, dat in self.lines_b.items() }

        while self.adjust_snaps(self.mesh_a, self.mesh_b, self.lines_a, 'a') > 0 \
            or self.adjust_snaps(self.mesh_b, self.mesh_a, self.lines_b, 'b') > 0:

            pass

        self.final_check(self.lines_a, 'a')
        self.final_check(self.lines_b, 'b')


    def final_check(self, connected, name):
        print(f'final_check({name})')

        data = []

        for edge, dat in connected.items():
            for ind, dists in dat.items():

                if len(dists) > 1 and dists[-1] > 1e-10:
                    data.append((edge, ind, math.sqrt(dists[-1])))

        data.sort(key=itemgetter(2))

        for dat in data:
            print(*dat)


    def adjust_snaps(self, mesh, other_mesh, connected, name):
        i = 0

        for edge, dat in connected.items():
            # edge auf die gesnapt wurde

            for ind, dists in dat.items():
                # ind ist der punkt aus mesh, der auf edge liegen sollte

                pt = mesh.GetPoint(ind)

                proj, d, l, t = _proj_line(other_mesh, edge, pt)

                if d > 1e-10:
                     # konvergenz zur 0 ist erstrebenswert

                    if d < dists[-1]:
                        Prepare.move_pt(mesh, ind, proj)

                        dists.append(d)

                        i += 1

        return i



    def align_congruent_points(self):
        tree = vtkKdTreePointLocator()

        tree.SetDataSet(self.mesh_b)
        tree.BuildLocator()

        for i in range(self.mesh_a.GetNumberOfPoints()):
            p = self.mesh_a.GetPoint(i)

            d = reference(0)

            j = tree.FindClosestPointWithinRadius(1e-5, p, d)

            if j != -1:
                Prepare.move_pt(self.mesh_a, i, self.mesh_b.GetPoint(j))

                self.congr_ids_a.add(i)
                self.congr_ids_b.add(j)

    def find_inters(self, mesh, other_mesh, omit_ids, connected, other_connected, name, deb_id = None):
        if deb_id is not None:
            print(f'find_iters({name}, {deb_id})')

        cell_lines = defaultdict(list)

        lines = set()

        itr = mesh.GetPolys().NewIterator()

        while not itr.IsDoneWithTraversal():
            curr = itr.GetCurrentCell()

            ids = [ curr.GetId(i) for i in range(curr.GetNumberOfIds()) ]

            ids = ids + ids[:1]

            for line in zip(ids, ids[1:]):
                _line = tuple(sorted(line))

                cell_lines[itr.GetCurrentCellId()].append(_line)

                lines.add(_line)

            itr.GoToNextCell()

        _lines = []

        if deb_id is not None:
            _lines = cell_lines[deb_id]

            print(_lines)

        tree = vtkModifiedBSPTree()
        tree.SetDataSet(other_mesh)
        tree.BuildLocator()

        pd = vtkPolyData()
        pd.Allocate(1)

        pd_pts = vtkPoints()

        inters = defaultdict(list)

        other_inters = defaultdict(list)

        pts = vtkPoints()
        cells = vtkIdList()

        for line in lines:
            if line in _lines:
                print('line', line)

            a, b = line

            p_a = mesh.GetPoint(a)
            p_b = mesh.GetPoint(b)

            v = [0, 0, 0]
            vtkMath.Subtract(p_b, p_a, v)

            l = vtkMath.Normalize(v)

            q_b = [0, 0, 0]
            vtkMath.Add(p_b, v, q_b)

            q_a = [0, 0, 0]
            vtkMath.Subtract(p_a, v, q_a)

            if tree.IntersectWithLine(q_a, q_b, 1e-5, pts, cells) == 0:
                continue

            for i in range(pts.GetNumberOfPoints()):
                cell_id = cells.GetId(i)

                s = pts.GetPoint(i)

                # projektion s auf line, da schnittpunkt in abh. von cell nicht exakt sein kann

                s, d, *rest = _proj_line(mesh, line, s)

                if line in _lines:
                    print('->', cell_id, d)

                ids, cell_pts = get_points(other_mesh, cell_id)

                d_a = vtkMath.Distance2BetweenPoints(s, p_a)
                d_b = vtkMath.Distance2BetweenPoints(s, p_b)

                # abstand innerhalb 1e-5?

                if d_a < 1e-10 or d_b < 1e-10:

                    # lage innerhalb der cell?

                    snap = next(( SnapPoint(cell_id, s, id_, pt, line) for id_, pt in zip(ids, cell_pts) if vtkMath.Distance2BetweenPoints(s, pt) < 1e-10 ), None)

                    if not snap:
                        try:
                            edge, proj = proj_line(other_mesh, ids, s)

                            snap = SnapEdge(cell_id, s, edge, proj, line)
                        except TypeError:
                            pass

                    if snap:
                        if d_a < 1e-10 and a not in omit_ids:
                            inters[a].append(snap)

                        elif d_b < 1e-10 and b not in omit_ids:
                            inters[b].append(snap)

                        vert = vtkIdList()
                        vert.InsertNextId(pd_pts.InsertNextPoint(s))

                        pd.InsertNextCell(VTK_VERTEX, vert)

                else:

                    snap = next(( SnapPoint(cell_id, s, id_, pt, line) for id_, pt in zip(ids, cell_pts) if vtkMath.Distance2BetweenPoints(s, pt) < 1e-10 ), None)

                    if snap:
                        if snap.s != snap.pt:
                            other_inters[snap.id].append(snap)

                    else:

                        projs = [ (id_, pt, vtkMath.Distance2BetweenPoints(s, pt), *_proj_line(mesh, line, pt)) for id_, pt in zip(ids, cell_pts) ]

                        dists = { id_: (math.sqrt(d2), d) for id_, pt, d2, proj, d, l, t in projs if d < 1e-5 and t > 0 and t < l }

                        if snap is None:
                            if line in _lines:
                                print(cell_id, dists)

                        snap = next(( SnapPoint(cell_id, proj, id_, pt, line) for id_, pt, d2, proj, d, l, t in projs if d < 1e-5 and t > 0 and t < l ), None)

                        if snap and snap.s != snap.pt:
                            other_inters[snap.id].append(snap)




        pd.SetPoints(pd_pts)

        write(pd, f'verts_{name}_2.vtk')


        for ind, snaps in inters.items():

            if all( isinstance(snap, SnapEdge) for snap in snaps ):

                assert len(set( frozenset(snap.edge) for snap in snaps )) == 1

                snap = snaps[0]

                if snap.line in _lines:
                    print('1 ->', ind, snap)

                Prepare.move_pt(mesh, ind, snap.proj)

                connected[snap.edge].add(ind)

            elif all( isinstance(snap, SnapPoint) for snap in snaps ):
                assert len(set( snap.id for snap in snaps )) == 1

                snap = snaps[0]

                if snap.line in _lines:
                    print('2 ->', ind, snap)

                Prepare.move_pt(mesh, ind, snap.pt)

            else:
                point_snaps = [ snap for snap in snaps if isinstance(snap, SnapPoint) ]

                edge_snaps = [ snap for snap in snaps if isinstance(snap, SnapEdge) ]

                if point_snaps:
                    snap = point_snaps[0]

                    if snap.line in _lines:
                        print('4 ->', ind, snap)

                    Prepare.move_pt(mesh, ind, snap.pt)

                elif edge_snaps:
                    snap = edge_snaps[0]

                    if snap.line in _lines:
                        print('5 ->', ind, snap)

                    Prepare.move_pt(mesh, ind, snap.proj)

        mesh.RemoveDeletedCells()

        for ind, snaps in other_inters.items():
            snap, = snaps[:1]

            if snap.line in _lines:
                print('6 ->', ind, snap)

            Prepare.move_pt(other_mesh, ind, snap.s)

            other_connected[snap.line].add(snap.id)

            # cell_id,s,edge,proj,line

        write(mesh, f'new_pd_{name}_2.vtk')


if __name__ == '__main__':
    reader_a = vtkPolyDataReader()
    reader_a.SetFileName('fibulaA.vtk')

    reader_b = vtkPolyDataReader()
    reader_b.SetFileName('fibulaB.vtk')

    reader_a.Update()
    reader_b.Update()

    pd_a = reader_a.GetOutput()
    pd_b = reader_b.GetOutput()

    align = Align(pd_a, pd_b)

    # verify

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputData(0, align.mesh_a)
    bf.SetInputData(1, align.mesh_b)

    writer = vtkPolyDataWriter()
    writer.SetFileName('union_2.vtk')
    writer.SetInputConnection(bf.GetOutputPort())
    writer.Update()
