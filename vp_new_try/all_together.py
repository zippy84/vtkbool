#!/usr/bin/env python
# *-* coding: UTF-8 *-*

from vis_poly import vis_poly, to_abs_path, cross, is_near
from rm_trivials import rm_trivials
import json
from copy import deepcopy

E = 1e-5

def is_on_seg (a, b, pt):
    global E

    if is_near(a, pt) or is_near(b, pt):
        return False

    if (pt[0] < a[0] and pt[0] < b[0]) \
        or (pt[0] > a[0] and pt[0] > b[0]):
        return False

    if (pt[1] < a[1] and pt[1] < b[1]) \
        or (pt[1] > a[1] and pt[1] > b[1]):
        return False

    return abs(cross(a, b, pt)) < E

def add_internals (pts, poly):
    num = len(poly)

    res = []

    for i in range(num):
        p_a = poly[i]
        p_b = poly[(i+1)%num]

        res.append(p_a)

        for j, pt in enumerate(pts):
            if is_on_seg(p_a['pt'], p_b['pt'], pt):
                res.append({ 'pt': pt, 'idx': j })

    return res

def vis_poly_wrapper (poly, idx):
    pts = deepcopy(poly)
    res = rm_trivials(pts, idx)
    p = vis_poly(res)

    return add_internals(pts, p)

with open('complex.json', 'r') as f:
    polys = json.load(f)['polys']

    for i, poly in enumerate(polys):
        num = len(poly)

        for j in range(1, num):
            poly[j][0] += poly[j-1][0]
            poly[j][1] += poly[j-1][1]

        all_res = {}

        for j in range(num):
            all_res[j] = to_abs_path([ p['pt'] for p in vis_poly_wrapper(poly, j) ])

        with open('data_files/data_%i.js' % i, 'w') as out:
            out.write('var polys = %s;' % json.dumps(all_res))
