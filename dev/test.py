#!/usr/bin/env python
# *-* coding: UTF-8 *-*

import math
from collections import defaultdict, namedtuple
from operator import itemgetter

from vtkmodules.vtkCommonCore import vtkIdList, vtkIdTypeArray, vtkPoints, vtkMath
from vtkmodules.vtkCommonDataModel import vtkPolyData, VTK_POLYGON, VTK_VERTEX
from vtkmodules.vtkFiltersSources import vtkCylinderSource
from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter
from vtkmodules.vtkCommonTransforms import vtkTransform
from vtkmodules.vtkIOLegacy import vtkPolyDataWriter, vtkPolyDataReader
from vtkmodules.vtkFiltersFlowPaths import vtkModifiedBSPTree
from vtkmodules.vtkFiltersCore import vtkCleanPolyData

SnapPoint = namedtuple('SnapPoint', 'cell_id,s,id,pt,line')
SnapEdge = namedtuple('SnapEdge', 'cell_id,s,edge,proj,line')
SnapPoly = namedtuple('SnapPoly', 'cell_id,s,line')

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

        if d < 1e-5:
            return edge, proj


def get_points(pd, cell_id):
    cell = vtkIdList()
    pd.GetCellPoints(cell_id, cell)

    ids = [ cell.GetId(id_) for id_ in range(cell.GetNumberOfIds()) ]

    return ids, [ pd.GetPoint(id_) for id_ in ids ]

def clean(pd):
    clean = vtkCleanPolyData()
    clean.SetInputData(pd)
    clean.Update()

    return clean.GetOutput()

class Prepare:
    def __init__(self, mesh_a, mesh_b):
        self.mesh_a = clean(mesh_a)
        self.mesh_b = clean(mesh_b)

        self.find(self.mesh_a, self.mesh_b, 'verts_a.vtk')
        self.find(self.mesh_b, self.mesh_a, 'verts_b.vtk')

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

                for id_, pt_ in zip(ids, cell_pts):
                    if vtkMath.Distance2BetweenPoints(pt, pt_) < 1e-10:
                        k = ','.join([ f'{x:z.5f}' for x in pt_ ])

                        point_snaps[k].append(SnapPoint(cell_id, pt, id_, pt_, line))

                        break

                try:
                    edge, proj = proj_line(other_mesh, ids, pt)

                    k = ','.join([ f'{x:z.5f}' for x in proj ])

                    edge_snaps[k].append(SnapEdge(cell_id, pt, edge, proj, line))

                except TypeError:
                    pass

            point_snaps = { k: v for k, v in point_snaps.items() if len(v) > 1 }

            print(point_snaps)

            for snaps in point_snaps.values():

                assert len(set( snap.id for snap in snaps )) == 1

                proj, d = _proj_line(mesh, snaps[0].line, snaps[0].pt)

                print(snaps[0].id, '->', proj, d)


            edge_snaps = { k: v for k, v in edge_snaps.items() if len(v) > 1 }

            print(edge_snaps)

            for snaps in edge_snaps.values():

                assert len(set( frozenset(snap.edge) for snap in snaps )) == 1

                proj, d = _proj_line(mesh, snaps[0].line, snaps[0].proj)

                all_edge_snaps[frozenset(snaps[0].edge)].append((snaps[0].proj, proj, d))


        pd.SetPoints(pd_pts)

        write(pd, file_name)

        for edge, projs in all_edge_snaps.items():
            a, b = edge

            p = other_mesh.GetPoint(a)

            projs[:] = [ (src_proj, dest_proj, d, vtkMath.Distance2BetweenPoints(p, src_proj)) for src_proj, dest_proj, d in projs ]

            projs.sort(key=itemgetter(2))

        print(all_edge_snaps)

if __name__ == '__main__':
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

    Prepare(pd_a, pd_b)
