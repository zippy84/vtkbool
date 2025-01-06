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

def get_plane(pts):
    normal = compute_normal(pts)

    d = vtkMath.Dot(normal, pts[0])

    return normal, d

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

        deb_a = 17920 # 17007
        deb_b = 1339

        self.edges_a = defaultdict(set)
        self.edges_b = defaultdict(set)

        self.polys_a = defaultdict(set)
        self.polys_b = defaultdict(set)

        self.find_inters(self.mesh_a, self.mesh_b, self.congr_ids_a, self.edges_a, self.edges_b, self.polys_a, 'a', deb_a)
        self.find_inters(self.mesh_b, self.mesh_a, self.congr_ids_b, self.edges_b, self.edges_a, self.polys_b, 'b', deb_b)

        deps_a = self.collect_dependencies(self.edges_b, 'a')
        deps_b = self.collect_dependencies(self.edges_a, 'b')

        self.edges_a = { edge: { ind: [sys.float_info.max] for ind in inds } for edge, inds in self.edges_a.items() }
        self.edges_b = { edge: { ind: [sys.float_info.max] for ind in inds } for edge, inds in self.edges_b.items() }

        while self.adjust_snaps(self.mesh_a, self.mesh_b, self.edges_a, deps_a, 'a') > 0 \
            or self.adjust_snaps(self.mesh_b, self.mesh_a, self.edges_b, deps_b, 'b') > 0:

            pass

        self.check_edge_projs(self.edges_a, 'a')
        self.check_edge_projs(self.edges_b, 'b')

        self.check_poly_projs(self.mesh_a, self.mesh_b, self.polys_a, deps_a, 'a')
        self.check_poly_projs(self.mesh_b, self.mesh_a, self.polys_b, deps_b, 'b')

    def collect_dependencies(self, edges, name):
        print(f'collect_dependencies({name})')

        deps = defaultdict(lambda: defaultdict(set))

        for edge, inds in edges.items():
            edge_ = tuple(sorted(edge))

            a, b = edge_

            for ind in inds:
                deps[a][edge_].add(ind)
                deps[b][edge_].add(ind)

        for k, v in deps.items():
            print(k, '->', sum(map(list, v.values()), []))

        return deps


    def check_poly_projs(self, mesh, other_mesh, polys, deps, name):
        print(f'check_poly_projs({name})')

        data = []

        for cell_id, inds in polys.items():
            ids, cell_pts = get_points(other_mesh, cell_id)

            plane_normal, plane_d = get_plane(cell_pts)

            for ind in inds:
                pt = mesh.GetPoint(ind)

                d = vtkMath.Dot(plane_normal, pt)-plane_d

                data.append((cell_id, ind, abs(d)))

                if ind in deps:
                    print(ind, '-->', dict(deps[ind]))

        data.sort(key=itemgetter(2))

        for dat in data:
            print(*dat)


    def check_edge_projs(self, edges, name):
        print(f'check_edge_projs({name})')

        data = []

        for edge, dat in edges.items():
            for ind, dists in dat.items():

                if len(dists) > 1 and dists[-1] > 1e-10:
                    data.append((edge, ind, math.sqrt(dists[-1])))

        data.sort(key=itemgetter(2))

        for dat in data:
            print(*dat)


    def adjust_snaps(self, mesh, other_mesh, edges, deps, name):
        i = 0

        for edge, dat in edges.items():
            # edge auf die gesnapt wurde

            for ind, dists in dat.items():
                # ind ist der punkt aus mesh, der auf edge liegen sollte

                pt = mesh.GetPoint(ind)

                if ind in deps:
                    print(ind, '-->', dict(deps[ind]))

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

    def find_inters(self, mesh, other_mesh, omit_ids, edges, other_edges, polys, name, deb_id = None):
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

                if line in _lines:
                    print(d_a, d_b)

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

                        if not snap:
                            plane_normal, plane_d = get_plane(cell_pts)

                            plane_d_a = vtkMath.Dot(plane_normal, p_a)-plane_d
                            plane_d_b = vtkMath.Dot(plane_normal, p_b)-plane_d

                            if line in _lines:
                                print(plane_d_a, plane_d_b)

                            base = Base(cell_pts, plane_normal)

                            if abs(plane_d_a) < 1e-5:
                                v = plane_normal[:]

                                vtkMath.MultiplyScalar(v, plane_d_a)

                                w = [0, 0, 0]

                                vtkMath.Subtract(p_a, v, w)

                                if in_poly(base.poly, base.tr_forward(w)):
                                    snap = SnapPoly(cell_id, w, line)

                                    if line in _lines:
                                        print(a, snap)

                                    inters[a].append(snap)

                            if abs(plane_d_b) < 1e-5:
                                v = plane_normal[:]

                                vtkMath.MultiplyScalar(v, plane_d_b)

                                w = [0, 0, 0]

                                vtkMath.Subtract(p_b, v, w)

                                if in_poly(base.poly, base.tr_forward(w)):
                                    snap = SnapPoly(cell_id, w, line)

                                    if line in _lines:
                                        print(b, snap)

                                    inters[b].append(snap)




        pd.SetPoints(pd_pts)

        write(pd, f'verts_{name}_2.vtk')


        for ind, snaps in inters.items():

            if all( isinstance(snap, SnapPoint) for snap in snaps ):
                assert len(set( snap.id for snap in snaps )) == 1

                snap = snaps[0]

                if snap.line in _lines:
                    print('1 ->', ind, snap)

                Prepare.move_pt(mesh, ind, snap.pt)

            elif all( isinstance(snap, SnapEdge) for snap in snaps ):

                assert len(set( frozenset(snap.edge) for snap in snaps )) == 1

                snap = snaps[0]

                if snap.line in _lines:
                    print('2 ->', ind, snap)

                Prepare.move_pt(mesh, ind, snap.proj)

                edges[snap.edge].add(ind)

            elif all( isinstance(snap, SnapPoly) for snap in snaps ):
                snap = snaps[0]

                if snap.line in _lines:
                    print('3 ->', ind, snap)

                Prepare.move_pt(mesh, ind, snap.s)

                polys[snap.cell_id].add(ind)

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

            other_edges[snap.line].add(snap.id)


        write(mesh, f'new_pd_{name}_2.vtk')


if __name__ == '__main__':
    reader_a = vtkPolyDataReader()
    reader_a.SetFileName('fibulaA.vtk')

    reader_b = vtkPolyDataReader()
    reader_b.SetFileName('fibulaB.vtk')

    # reader_a = vtkPolyDataReader()
    # reader_a.SetFileName('../_bugs/issue77/a.vtk')

    # reader_b = vtkPolyDataReader()
    # reader_b.SetFileName('../_bugs/issue77/b.vtk')

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
