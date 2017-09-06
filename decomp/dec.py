#!/usr/bin/env python
# *-* coding: UTF-8 *-*

import sys
import json

from collections import deque

sys.path.append('../vp_new_try')

from vis_poly import vis_poly_wrapper
from rm_trivials import rm_internals
from tools import cross

def decompose (pts):
    poly = deque([ { 'pt': pt, 'idx': i } for i, pt in enumerate(pts) ])

    rm_internals(poly)

    num = len(poly)

    def is_refl (a, b, c):
        return cross(poly[a]['pt'], poly[b]['pt'], poly[c]['pt']) > 0

    for i in range(num):
        j = (i+num-1)%num
        k = (i+1)%num

        poly[i]['refl'] = is_refl(i, j, k)

    first = next(( i for i, p in enumerate(poly) if p['refl'] ))

    poly.rotate(-first)

    print poly[0]

    poly_ = [ p['pt'] for p in poly ]

    pairs = set()

    for i, p in enumerate(poly):
        if p['refl']:
            vp = filter(lambda v: v['idx'] is not None, vis_poly_wrapper(poly_, i))
            # idx bezieht sich auf poly_

            for v in vp[1:]:
                pairs.add(tuple(sorted([vp[0]['idx'], v['idx']])))

    print pairs

with open('../vp_new_try/complex.json', 'r') as f:
    polys = json.load(f)['polys']

    for i, poly in enumerate(polys):
        if i != 1:
            continue

        num = len(poly)

        for j in range(1, num):
            poly[j][0] += poly[j-1][0]
            poly[j][1] += poly[j-1][1]

        decompose(poly)
