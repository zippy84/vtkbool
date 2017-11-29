#!/usr/bin/env python
# *-* coding: UTF-8 *-*

# export LD_LIBRARY_PATH=/home/zippy/VTK7/lib

import sys
import json
from collections import deque

sys.path.extend(['/home/zippy/VTK6/lib/python2.7/site-packages'])
import vtk

sys.path.append('../../vp/dev')
from tools import intersect2, get_angle, is_cw

E = 1e-5

# funktioniert nicht mit holes.json, da die holes selber holes enthalten

with open('holes2.json', 'r') as f:
    polys = json.load(f)['holes']

    for i, p in enumerate(polys):
        num = len(p)

        for j in range(1, num):
            p[j][0] += p[j-1][0]
            p[j][1] += p[j-1][1]

        assert not is_cw(p)


    pts = vtk.vtkPoints()

    src = []

    for i, p in enumerate(polys):
        map(pts.InsertNextPoint, [ pt + [0] for pt in p ])
        l = len(p)
        src.extend([i]*l)

    print pts

    num_pts = pts.GetNumberOfPoints()

    tree = vtk.vtkKdTree()
    tree.OmitZPartitioning()
    tree.BuildLocatorFromPoints(pts)
    print tree

    print src

    edges = []

    o = 0

    os = []
    defs = {}

    for i, p in enumerate(polys):
        l = len(p)
        for j in range(l):
            edges.append([o+j, o+(j+1)%l])

        os.append(l)

        defs[i] = [ o+j for j in range(l) ]

        o += l

    print edges

    results = []

    for p in polys:
        results.append({})

    ids = range(num_pts)
    curr = 1

    num_polys = len(polys)

    print num_pts, num_polys, len(edges)

    cons = []

    while len(ids) > 0:
        for pt_id in ids:
            print 'pt_id', pt_id

            print pts.GetPoint(pt_id)

            found = vtk.vtkIdList()
            tree.FindClosestNPoints(curr*5, pts.GetPoint(pt_id), found)

            #print found

            valids = []

            for i in range((curr-1)*5, found.GetNumberOfIds()):
                ind = found.GetId(i)
                #print ind, src[ind]

                # darf nicht zum gleichen poly gehören
                if src[ind] != src[pt_id]:
                    valids.append(ind)


            print 'valids', valids

            invalids = []

            for v in valids:
                #print '%i --' % v

                for edge in edges:
                    if pt_id in edge \
                        or v in edge:

                        continue

                    ps = map(pts.GetPoint, [pt_id, v, edge[0], edge[1]])
                    d = intersect2(*ps)

                    if d:
                        invalids.append(v)
                        break

            print invalids

            for v in invalids:
                valids.remove(v)

            for v in valids:
                #print v, src[v]

                con = tuple([pt_id, v])

                pt_a = pts.GetPoint(pt_id)
                pt_b = pts.GetPoint(v)

                vec = [0, 0, 0]
                vtk.vtkMath.Subtract(pt_a, pt_b, vec)

                d = vtk.vtkMath.Norm(vec)

                p_a = src[pt_id]
                p_b = src[v]

                results[p_a][con] = d
                results[p_b][tuple(reversed(con))] = d


        del ids[:]
        del cons[:]

        viewed = [0]

        while len(viewed) < num_polys:
            f = []

            g = None

            for n in viewed:
                for con, d in results[n].items():
                    if src[con[1]] not in viewed:
                        if not g or d < g['d']:
                            g = {'d': d, 'con': con}

            if g:
                viewed.append(src[g['con'][1]])
                cons.append(g['con'])

            else:
                for i, p in enumerate(polys):
                    if i not in viewed:
                        _o = sum(os[:i])
                        for j in range(os[i]):
                            ids.append(_o+j)

                print ids

                curr += 1

                break

    print curr

    print cons

    def _fct (poly, a, b, s):
        pt_b = pts.GetPoint(b)
        num = len(poly)

        for i, id_ in enumerate(poly):
            if id_ == a:
                pt_a = pts.GetPoint(id_)

                v = [0, 0, 0]
                vtk.vtkMath.Subtract(pt_b, pt_a, v)
                vtk.vtkMath.Normalize(v)

                i_a = poly[(i+1)%num]
                i_b = poly[(i+num-1)%num]

                p_a = pts.GetPoint(i_a)
                p_b = pts.GetPoint(i_b)

                w_a = [0, 0, 0]
                vtk.vtkMath.Subtract(p_a, pt_a, w_a)
                vtk.vtkMath.Normalize(w_a)

                w_b = [0, 0, 0]
                vtk.vtkMath.Subtract(p_b, pt_a, w_b)
                vtk.vtkMath.Normalize(w_b)

                ref = get_angle(w_a, w_b)
                phi = get_angle(w_a, v)

                if phi > ref or (s == 0 and phi < ref):
                    break

        return i

    repls = dict([ (i, i) for i in range(num_polys) ])

    print repls
    print defs

    num = num_polys

    for con in cons:
        a, b = con
        p_a = repls[src[a]]
        p_b = repls[src[b]]

        d_a = defs[p_a]
        d_b = defs[p_b]

        print p_a, p_b

        #i_a = d_a.index(a)
        #i_b = d_b.index(b)

        i_a = _fct(d_a, a, b, p_a)
        i_b = _fct(d_b, b, a, p_b)

        q_a = deque(d_a)
        q_a.rotate(-i_a)

        q_b = deque(d_b)
        q_b.rotate(-i_b)

        #print q_a
        #print q_b

        l_a = list(q_a)
        l_b = list(q_b)

        assert l_a[0] == a
        assert l_b[0] == b

        m = []

        if p_a == 0:
            m.append(l_a[0])
            m.extend(reversed(l_a))
        else:
            m.extend(l_a)
            m.append(l_a[0])

        if p_b == 0:
            m.append(l_b[0])
            m.extend(reversed(l_b))
        else:
            m.extend(l_b)
            m.append(l_b[0])

        #print m

        defs[num] = m

        repls[src[a]] = num
        repls[src[b]] = num

        num += 1

    print defs
    print repls

    print num

    svg = []

    for idx in defs[num-1]:
        pt = [0, 0, 0]
        pts.GetPoint(idx, pt)
        svg.append('%f,%f' % (pt[0], pt[1]))

    print 'M' + ' L'.join(svg) + 'Z'
