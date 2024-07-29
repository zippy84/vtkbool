#!/usr/bin/env python
# *-* coding: UTF-8 *-*

import math
from collections import defaultdict, namedtuple
from operator import itemgetter

from vtkmodules.vtkCommonCore import vtkIdList, vtkIdTypeArray, vtkPoints, vtkMath
from vtkmodules.vtkCommonDataModel import vtkPolyData, VTK_POLYGON, VTK_VERTEX, VTK_TRIANGLE
from vtkmodules.vtkFiltersSources import vtkCylinderSource
from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter
from vtkmodules.vtkCommonTransforms import vtkTransform
from vtkmodules.vtkIOLegacy import vtkPolyDataWriter, vtkPolyDataReader
from vtkmodules.vtkFiltersFlowPaths import vtkModifiedBSPTree
from vtkmodules.vtkFiltersCore import vtkCleanPolyData

from test import SnapPoint, SnapEdge, SnapPoly, write, compute_normal, Base, in_poly, proj_line, get_points, clean

class Align:
    def __init__(self, mesh_a, mesh_b):
        self.mesh_a = clean(mesh_a)
        self.mesh_b = clean(mesh_b)

        self.find_inters(self.mesh_a, self.mesh_b, 'verts_a_2.vtk')
        self.find_inters(self.mesh_b, self.mesh_a, 'verts_b_2.vtk')

    def find_inters(self, mesh, other_mesh, file_name):

        lines = set()

        itr = mesh.GetPolys().NewIterator()

        while not itr.IsDoneWithTraversal():
            curr = itr.GetCurrentCell()

            ids = [ curr.GetId(i) for i in range(curr.GetNumberOfIds()) ]

            ids = ids + ids[:1]

            for line in zip(ids, ids[1:]):
                if line[::-1] not in lines:
                    lines.add(tuple(line))

            itr.GoToNextCell()

        # print(lines)

        tree = vtkModifiedBSPTree()
        tree.SetDataSet(other_mesh)
        tree.BuildLocator()

        pd = vtkPolyData()
        pd.Allocate(1)

        pd_pts = vtkPoints()

        inters = defaultdict(list)

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

                pt = pts.GetPoint(i)

                d_a = vtkMath.Distance2BetweenPoints(pt, p_a)
                d_b = vtkMath.Distance2BetweenPoints(pt, p_b)

                # abstand innerhalb (1e-6, 1e-5)?

                if (d_a < 1e-10 and d_a > 1e-12) or (d_b < 1e-10 and d_b > 1e-12):

                    # lage innerhalb der cell?

                    snap = None

                    ids, cell_pts = get_points(other_mesh, cell_id)

                    for id_, pt_ in zip(ids, cell_pts):
                        if vtkMath.Distance2BetweenPoints(pt, pt_) < 1e-10:
                            snap = SnapPoint(cell_id, pt, id_, pt_, line)

                            break

                    if not snap:
                        try:
                            edge, proj = proj_line(other_mesh, ids, pt)

                            snap = SnapEdge(cell_id, pt, edge, proj, line)
                        except TypeError:
                            pass

                    if not snap:
                        normal = compute_normal(cell_pts)

                        base = Base(cell_pts, normal)

                        if in_poly(base.poly, base.tr_forward(pt)):
                            snap = SnapPoly(cell_id, pt, line)

                    if snap:
                        # print(snap)

                        if d_a < 1e-10 and d_a > 1e-12:
                            inters[a].append(snap)

                        elif d_b < 1e-10 and d_b > 1e-12:
                            inters[b].append(snap)

                        vert = vtkIdList()
                        vert.InsertNextId(pd_pts.InsertNextPoint(pt))

                        pd.InsertNextCell(VTK_VERTEX, vert)


        pd.SetPoints(pd_pts)

        write(pd, file_name)

        for ind, snaps in inters.items():

            if all( isinstance(snap, SnapEdge) for snap in snaps ):

                assert len(set( frozenset(snap.edge) for snap in snaps )) == 1

                self.move_pt(mesh, ind, snaps[0].proj)

            elif all( isinstance(snap, SnapPoint) for snap in snaps ):
                assert len(set( snap.id for snap in snaps )) == 1

                self.move_pt(mesh, ind, snaps[0].pt)

            elif all( isinstance(snap, SnapPoly) for snap in snaps ):
                self.move_pt(mesh, ind, snaps[0].s)

            else:
                point_snaps = [ snap for snap in snaps if isinstance(snap, SnapPoint) ]

                edge_snaps = [ snap for snap in snaps if isinstance(snap, SnapEdge) ]

                if point_snaps:
                    self.move_pt(mesh, ind, point_snaps[0].pt)

                elif edge_snaps:
                    self.move_pt(mesh, ind, edge_snaps[0].proj)


    def move_pt(self, mesh, ind, dest_pt):
        src_pt = mesh.GetPoint(ind)

        d = vtkMath.Distance2BetweenPoints(src_pt, dest_pt)

        print(ind, '->', dest_pt, d)

        cells = vtkIdList()

        mesh.GetPointCells(ind, cells)

        types = [ mesh.GetCellType(cells.GetId(i)) for i in range(cells.GetNumberOfIds()) ]

        if all(t == VTK_TRIANGLE for t in types):
            pass


if __name__ == '__main__':
    reader_a = vtkPolyDataReader()
    reader_a.SetFileName('fibulaA.vtk')

    reader_b = vtkPolyDataReader()
    reader_b.SetFileName('fibulaB.vtk')

    reader_a.Update()
    reader_b.Update()

    pd_a = reader_a.GetOutput()
    pd_b = reader_b.GetOutput()

    Align(pd_a, pd_b)
