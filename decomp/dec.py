#!/usr/bin/env python
# *-* coding: UTF-8 *-*

import sys
import json

sys.path.append('../vp_new_try')

from vis_poly import vis_poly_wrapper
from tools import cross

def decompose (pts):
    pass

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
