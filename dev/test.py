#!/usr/bin/env python
# *-* coding: UTF-8 *-*

import math
from collections import defaultdict, namedtuple
from operator import itemgetter
from itertools import batched

from vtkmodules.vtkCommonCore import vtkIdList, vtkIdTypeArray, vtkPoints, vtkMath
from vtkmodules.vtkCommonDataModel import vtkPolyData, vtkPolygon, VTK_POLYGON, VTK_QUAD, VTK_TRIANGLE, VTK_VERTEX
from vtkmodules.vtkFiltersSources import vtkCylinderSource
from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter
from vtkmodules.vtkCommonTransforms import vtkTransform
from vtkmodules.vtkIOLegacy import vtkPolyDataWriter, vtkPolyDataReader
from vtkmodules.vtkFiltersFlowPaths import vtkModifiedBSPTree
from vtkmodules.vtkFiltersCore import vtkCleanPolyData
from vtkmodules.vtkCommonExecutionModel import vtkAlgorithm

import sys
sys.path.extend(['/home/zippy/vtkbool/build/lib/python3.12/site-packages/vtkbool'])

from vtkBool import vtkPolyDataBooleanFilter

SnapPoint = namedtuple('SnapPoint', 'cell_id,s,id,pt,line')
SnapEdge = namedtuple('SnapEdge', 'cell_id,s,edge,proj,line')
SnapPoly = namedtuple('SnapPoly', 'cell_id,s,line')

Point = namedtuple('Point', 'id,pt')

def write(pd, file_name):
    writer = vtkPolyDataWriter()
    writer.SetInputData(pd)
    writer.SetFileName(file_name)
    writer.Update()

def compute_normal(pts):
    n = [0, 0, 0]

    pts_ = pts + pts[:1]

    for pA, pB in zip(pts_, pts_[1:]):
        n[0] += (pA[1]-pB[1])*(pA[2]+pB[2])
        n[1] += (pA[2]-pB[2])*(pA[0]+pB[0])
        n[2] += (pA[0]-pB[0])*(pA[1]+pB[1])

    vtkMath.Normalize(n)

    return n

class Base:
    def __init__(self, pts, n):
        self.pts = pts

        a, b = pts[:2]

        ei = [b[0]-a[0], b[1]-a[1], b[2]-a[2]]

        vtkMath.Normalize(ei)

        ej = [n[1]*ei[2]-n[2]*ei[1], -n[0]*ei[2]+n[2]*ei[0], n[0]*ei[1]-n[1]*ei[0]]

        vtkMath.Normalize(ej)

        d = n[0]*a[0]+n[1]*a[1]+n[2]*a[2]

        self.ei = ei
        self.ej = ej
        self.d = d

        self.poly = [ self.tr_forward(pt) for pt in self.pts ]

    def tr_forward(self, pt):
        return [self.ei[0]*pt[0]+self.ei[1]*pt[1]+self.ei[2]*pt[2], self.ej[0]*pt[0]+self.ej[1]*pt[1]+self.ej[2]*pt[2]]

    def tr_backward(self, pt):
        return [pt[0]*self.ei[0]+pt[1]*self.ej[0]+self.d*self.n[0],
            pt[0]*self.ei[1]+pt[1]*self.ej[1]+self.d*self.n[1],
            pt[0]*self.ei[2]+pt[1]*self.ej[2]+self.d*self.n[2]]

def in_poly(poly, pt):
    poly_ = poly + poly[:1]

    _in = False

    for a, b in zip(poly_, poly_[1:]):
        if (a[0] <= pt[0] or b[0] <= pt[0]) and ((a[1] < pt[1] and b[1] >= pt[1]) or (b[1] < pt[1] and a[1] >= pt[1])):
            if a[0]+(pt[1]-a[1])*(b[0]-a[0])/(b[1]-a[1]) < pt[0]:
                _in = not _in

    return _in

def _proj_line(mesh, edge, pt):
    a, b = edge

    p_a = mesh.GetPoint(a)
    p_b = mesh.GetPoint(b)

    v = [0, 0, 0]

    vtkMath.Subtract(p_b, p_a, v)
    vtkMath.Normalize(v)

    w = [0, 0, 0]

    vtkMath.Subtract(pt, p_a, w)

    l = vtkMath.Normalize(w)

    prod = min(max(vtkMath.Dot(v, w), -1), 1)

    d = math.sin(math.acos(prod))*l

    proj = [0, 0, 0]

    v_ = v[:]

    vtkMath.MultiplyScalar(v_, math.sqrt(l*l-d*d))

    vtkMath.Add(p_a, v_, proj)

    return proj, d

def proj_line(mesh, ids, pt):
    ids = ids + ids[:1]

    for edge in zip(ids, ids[1:]):
        proj, d = _proj_line(mesh, edge, pt)

        if d > 0 and d < 1e-5:
            return edge, proj


def get_points(pd, cell_id):
    cell = vtkIdList()
    pd.GetCellPoints(cell_id, cell)

    ids = [ cell.GetId(i) for i in range(cell.GetNumberOfIds()) ]

    return ids, [ pd.GetPoint(id_) for id_ in ids ]

def clean(pd):
    clean = vtkCleanPolyData()
    clean.SetOutputPointsPrecision(vtkAlgorithm.DOUBLE_PRECISION)
    clean.SetInputData(pd)
    clean.Update()

    return clean.GetOutput()

class Prepare:
    def __init__(self, mesh_a, mesh_b):
        self.mesh_a = clean(mesh_a)
        self.mesh_b = clean(mesh_b)

        write(self.mesh_a, 'cleaned_a.vtk')
        write(self.mesh_b, 'cleaned_b.vtk')

        self.find(self.mesh_a, self.mesh_b, 'a')
        self.find(self.mesh_b, self.mesh_a, 'b')

    def find(self, mesh, other_mesh, file_name):

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

        pts = vtkPoints()
        cells = vtkIdList()

        all_edge_snaps = defaultdict(list)

        for line in lines:
            # print(line)

            a, b = line

            p_a = mesh.GetPoint(a)
            p_b = mesh.GetPoint(b)

            if tree.IntersectWithLine(p_a, p_b, 1e-5, pts, cells) == 0:
                continue

            point_snaps = defaultdict(list)

            edge_snaps = defaultdict(list)

            for i in range(pts.GetNumberOfPoints()):
                pt = pts.GetPoint(i)

                cell_id = cells.GetId(i)

                ids, cell_pts = get_points(other_mesh, cell_id)

                normal = compute_normal(cell_pts)

                base = Base(cell_pts, normal)

                if not in_poly(base.poly, base.tr_forward(pt)):
                    continue

                vert = vtkIdList()
                vert.InsertNextId(pd_pts.InsertNextPoint(pt))

                pd.InsertNextCell(VTK_VERTEX, vert)

                snap = next(( SnapPoint(cell_id, pt, id_, pt_, line) for id_, pt_ in zip(ids, cell_pts) if vtkMath.Distance2BetweenPoints(pt, pt_) < 1e-10 ), None)

                if snap:
                    point_snaps[snap.id].append(snap)

                else:
                    try:
                        edge, proj = proj_line(other_mesh, ids, pt)

                        k = ','.join([ f'{x:z.5f}' for x in proj ])

                        snap = SnapEdge(cell_id, pt, edge, proj, line)

                        edge_snaps[k].append(snap)

                    except TypeError:
                        pass

            point_snaps = [ snaps for snaps in point_snaps.values() if len(snaps) > 1 ]

            print(point_snaps)

            for snaps in point_snaps:
                assert len(snaps) == 2 # sind hier nicht auch mehr als 2 möglich?

                a, b = snaps

                assert a.id == b.id

                proj, d = _proj_line(mesh, a.line, a.pt)

                print(a.id, '->', proj, d)

                Prepare.move_pt(other_mesh, a.id, proj)


            edge_snaps = [ snaps for snaps in edge_snaps.values() if len(snaps) > 1 ]

            print(edge_snaps)

            for snaps in edge_snaps:
                assert len(snaps) == 2

                a, b = snaps

                assert frozenset(a.edge) == frozenset(b.edge)

                proj, d = _proj_line(mesh, a.line, a.proj)

                all_edge_snaps[frozenset(a.edge)].append((a, proj, d))


        pd.SetPoints(pd_pts)

        write(pd, f'verts_{file_name}.vtk')

        _cells = defaultdict(dict)

        other_mesh.BuildLinks()

        other_pts = other_mesh.GetPoints()

        for projs in all_edge_snaps.values():
            first, *_ = projs[0]

            print(first)

            a, b = first.edge

            p = other_mesh.GetPoint(a)

            projs[:] = [ (snap, dest_proj, d, vtkMath.Distance2BetweenPoints(p, snap.proj)) for snap, dest_proj, d in projs ]

            projs.sort(key=itemgetter(3))

            neigs = vtkIdList()

            other_mesh.GetCellEdgeNeighbors(first.cell_id, *first.edge, neigs)

            assert neigs.GetNumberOfIds() == 1

            neig_id = neigs.GetId(0)

            print(first.cell_id, neig_id)

            new_pts = [ Point(other_pts.InsertNextPoint(p[1]), p[1]) for p in projs ]

            print(new_pts)

            _cells[first.cell_id][first.edge] = new_pts

            _cells[neig_id][first.edge[::-1]] = new_pts[::-1]

        print(all_edge_snaps)

        print(_cells)

        for cell_id, edges in _cells.items():
            Prepare.tringulate_cell(other_mesh, cell_id, edges)

        other_mesh.RemoveDeletedCells()

        write(other_mesh, f'new_pd_{file_name}.vtk')


    @staticmethod
    def tringulate_cell(mesh, cell_id, edges = None):

        ids, pts = get_points(mesh, cell_id)

        normal = compute_normal(pts)

        base = Base(pts, normal)

        cell_pts = [ Point(id_, pt) for id_, pt in zip(ids, pts) ]

        new_cell = []

        if edges:
            cell_pts = cell_pts + cell_pts[:1]

            for a, b in zip(cell_pts, cell_pts[1:]):
                new_cell.append(a)

                edge = a.id, b.id

                if edge in edges:
                    new_cell.extend(edges[edge])

        else:
            new_cell = cell_pts[:]

        print([ p.id for p in new_cell ])

        poly = vtkPolygon()

        poly.GetPointIds().SetNumberOfIds(len(new_cell))
        poly.GetPoints().SetNumberOfPoints(len(new_cell))

        for i, p in enumerate(new_cell):
            poly.GetPointIds().SetId(i, i)

            pt = base.tr_forward(p.pt) + [0]

            poly.GetPoints().SetPoint(i, pt)

        triangles = vtkIdList()

        assert poly.Triangulate(triangles) != 0

        print(triangles.GetNumberOfIds())

        triangle_ids = [ triangles.GetId(i) for i in range(triangles.GetNumberOfIds()) ]

        mesh.DeleteCell(cell_id)

        for triangle in batched(triangle_ids, 3):
            cell = vtkIdList()

            [ cell.InsertNextId(new_cell[i].id) for i in triangle ]

            mesh.InsertNextCell(VTK_TRIANGLE, cell)


    @staticmethod
    def move_pt(mesh, ind, dest_pt):

        cells = vtkIdList()

        mesh.GetPointCells(ind, cells)

        for i in range(cells.GetNumberOfIds()):
            cell_id = cells.GetId(i)

            t = mesh.GetCellType(cell_id)

            if t == VTK_POLYGON or t == VTK_QUAD:
                Prepare.tringulate_cell(mesh, cell_id)

        mesh.GetPoints().SetPoint(ind, dest_pt)

        assert mesh.GetPoints().GetPoint(ind) == tuple(dest_pt)


def create_complex():
    cyl = vtkCylinderSource()
    cyl.SetHeight(2)
    cyl.SetResolution(12)

    tra_a = vtkTransform()
    tra_a.Translate(.25, 0, 0)

    tf_a = vtkTransformPolyDataFilter()
    tf_a.SetTransform(tra_a)
    tf_a.SetInputConnection(cyl.GetOutputPort())

    tra_b = vtkTransform()
    tra_b.Translate(-.25, 0, 0)

    tf_b = vtkTransformPolyDataFilter()
    tf_b.SetTransform(tra_b)
    tf_b.SetInputConnection(cyl.GetOutputPort())

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputConnection(0, tf_a.GetOutputPort())
    bf.SetInputConnection(1, tf_b.GetOutputPort())

    bf.Update()

    cleaned = clean(bf.GetOutput())

    writer = vtkPolyDataWriter()
    writer.SetFileName('complex.vtk')
    writer.SetInputData(cleaned)
    writer.Update()

def test():
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

    cyl.Update()
    tf.Update()

    pd_a = cyl.GetOutput()
    pd_b = tf.GetOutput()

    return Prepare(pd_a, pd_b)

def test2():
    cyl = vtkCylinderSource()
    cyl.SetHeight(2)
    cyl.SetResolution(12)

    reader = vtkPolyDataReader()
    reader.SetFileName('complex.vtk')

    z = .0000025

    tra = vtkTransform()
    tra.RotateZ(90)
    tra.Translate(0, 0, z)

    tf = vtkTransformPolyDataFilter()
    tf.SetTransform(tra)
    tf.SetInputConnection(reader.GetOutputPort())

    cyl.Update()
    tf.Update()

    pd_a = cyl.GetOutput()
    pd_b = tf.GetOutput()

    return Prepare(pd_a, pd_b)

def test3():
    reader = vtkPolyDataReader()
    reader.SetFileName('../testing/data/cross.vtk')

    z = .000001

    tra = vtkTransform()
    tra.RotateZ(45)
    tra.Translate(0, 0, z)

    tf = vtkTransformPolyDataFilter()
    tf.SetTransform(tra)
    tf.SetInputConnection(reader.GetOutputPort())

    tf.Update()

    pd_a = reader.GetOutput()
    pd_b = tf.GetOutput()

    return Prepare(pd_a, pd_b)


if __name__ == '__main__':
    r = sys.argv[1]

    if r == '0':
        prepare = test()
    elif r == '1':
        create_complex()
        prepare = test2()
    elif r == '2':
        prepare = test3()

    # verify

    bf = vtkPolyDataBooleanFilter()
    bf.SetInputData(0, prepare.mesh_a)
    bf.SetInputData(1, prepare.mesh_b)

    writer = vtkPolyDataWriter()
    writer.SetFileName('union.vtk')
    writer.SetInputConnection(bf.GetOutputPort())
    writer.Update()
