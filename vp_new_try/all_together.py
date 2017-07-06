#!/usr/bin/env python
# *-* coding: UTF-8 *-*

from vis_poly import vis_poly, to_abs_path
from rm_trivials import rm_trivials
import json
from copy import deepcopy

with open('complex.json', 'r') as f:
    polys = json.load(f)['polys']

    for i, poly in enumerate(polys):
        num = len(poly)

        for j in range(1, num):
            poly[j][0] += poly[j-1][0]
            poly[j][1] += poly[j-1][1]

        all_res = {}

        for j in range(num):
            pts = deepcopy(poly)

            #all_res[j] = to_abs_path(rm_trivials(pts, j))

            res = rm_trivials(pts, j)
            all_res[j] = to_abs_path(vis_poly(res))

        with open('data_files/data_%i.js' % i, 'w') as out:
            out.write('var polys = %s;' % json.dumps(all_res))
