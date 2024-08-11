#!/usr/bin/env python
# *-* coding: UTF-8 *-*

import math
from collections import defaultdict, namedtuple
from operator import itemgetter

from vtkmodules.vtkCommonCore import vtkIdList, vtkIdTypeArray, vtkPoints, vtkMath, mutable
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

        self.find_inters(self.mesh_a, self.mesh_b, self.congr_ids_a, 'a') # , 18340)
        self.find_inters(self.mesh_b, self.mesh_a, self.congr_ids_b, 'b') # , 1659)

    def align_congruent_points(self):
        tree = vtkKdTreePointLocator()

        tree.SetDataSet(self.mesh_b)
        tree.BuildLocator()

        for i in range(self.mesh_a.GetNumberOfPoints()):
            p = self.mesh_a.GetPoint(i)

            d = mutable(0)

            j = tree.FindClosestPointWithinRadius(1e-5, p, d)

            if j != -1:
                Prepare.move_pt(self.mesh_a, i, self.mesh_b.GetPoint(j))

                self.congr_ids_a.add(i)
                self.congr_ids_b.add(j)

    def find_inters(self, mesh, other_mesh, omit_ids, file_name, deb_id = None):
        if deb_id is not None:
            print(f'find_iters({file_name})')

        lines = set()

        itr = mesh.GetPolys().NewIterator()

        while not itr.IsDoneWithTraversal():
            curr = itr.GetCurrentCell()

            if deb_id is not None and deb_id != itr.GetCurrentCellId():
                pass

            else:

                ids = [ curr.GetId(i) for i in range(curr.GetNumberOfIds()) ]

                ids = ids + ids[:1]

                for line in zip(ids, ids[1:]):
                    if line[::-1] not in lines:
                        lines.add(tuple(line))

            itr.GoToNextCell()

        if deb_id is not None:
            print(lines)

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
            # print(line)

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

                s, d = _proj_line(mesh, line, s)

                if deb_id is not None:
                    print('->', cell_id, d)

                d_a = vtkMath.Distance2BetweenPoints(s, p_a)
                d_b = vtkMath.Distance2BetweenPoints(s, p_b)

                # abstand innerhalb 1e-5?

                if d_a < 1e-10 or d_b < 1e-10:

                    # lage innerhalb der cell?

                    ids, cell_pts = get_points(other_mesh, cell_id)

                    snap = next(( SnapPoint(cell_id, s, id_, pt, line) for id_, pt in zip(ids, cell_pts) if vtkMath.Distance2BetweenPoints(s, pt) < 1e-10 ), None)

                    if not snap:
                        try:
                            edge, proj = proj_line(other_mesh, ids, s)

                            snap = SnapEdge(cell_id, s, edge, proj, line)
                        except TypeError:
                            pass

                    if not snap:
                        normal = compute_normal(cell_pts)

                        base = Base(cell_pts, normal)

                        if in_poly(base.poly, base.tr_forward(s)):
                            snap = SnapPoly(cell_id, s, line)

                    if snap:
                        if d_a < 1e-10 and a not in omit_ids:
                            inters[a].append(snap)

                        elif d_b < 1e-10 and b not in omit_ids:
                            inters[b].append(snap)

                        vert = vtkIdList()
                        vert.InsertNextId(pd_pts.InsertNextPoint(s))

                        pd.InsertNextCell(VTK_VERTEX, vert)

                else:
                    ids, cell_pts = get_points(other_mesh, cell_id)

                    snap = next(( SnapPoint(cell_id, s, id_, pt, line) for id_, pt in zip(ids, cell_pts) if vtkMath.Distance2BetweenPoints(s, pt) < 1e-10 ), None)

                    if snap and snap.s != snap.pt:
                        other_inters[snap.id].append(snap)



        pd.SetPoints(pd_pts)

        write(pd, f'verts_{file_name}_2.vtk')

        for ind, snaps in inters.items():

            if all( isinstance(snap, SnapEdge) for snap in snaps ):

                assert len(set( frozenset(snap.edge) for snap in snaps )) == 1

                if deb_id is not None:
                    print(ind, snaps[0])

                Prepare.move_pt(mesh, ind, snaps[0].proj)

            elif all( isinstance(snap, SnapPoint) for snap in snaps ):
                assert len(set( snap.id for snap in snaps )) == 1

                if deb_id is not None:
                    print(ind, snaps[0])

                Prepare.move_pt(mesh, ind, snaps[0].pt)

            elif all( isinstance(snap, SnapPoly) for snap in snaps ):
                if deb_id is not None:
                    print(ind, snaps[0])

                Prepare.move_pt(mesh, ind, snaps[0].s)

            else:
                point_snaps = [ snap for snap in snaps if isinstance(snap, SnapPoint) ]

                edge_snaps = [ snap for snap in snaps if isinstance(snap, SnapEdge) ]

                if point_snaps:
                    if deb_id is not None:
                        print(ind, point_snaps[0])

                    Prepare.move_pt(mesh, ind, point_snaps[0].pt)

                elif edge_snaps:
                    if deb_id is not None:
                        print(ind, edge_snaps[0])

                    Prepare.move_pt(mesh, ind, edge_snaps[0].proj)

        mesh.RemoveDeletedCells()

        for ind, snaps in other_inters.items():
            snap, = snaps[:1]

            if deb_id is not None:
                print(ind, snap)

            Prepare.move_pt(other_mesh, ind, snap.s)

        write(mesh, f'new_pd_{file_name}_2.vtk')


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
