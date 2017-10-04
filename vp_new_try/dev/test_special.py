#!/usr/bin/env python
# *-* coding: UTF-8 *-*

import json
from copy import deepcopy
from tools import *
from rm_trivials import rm_internals, rm_trivials, add_internals

with open('special.json', 'r') as f:
    polys = json.load(f)['polys']

    for i, poly in enumerate(polys):
        if i != 1:
            continue

        num = len(poly)

        for j in range(1, num):
            poly[j][0] += poly[j-1][0]
            poly[j][1] += poly[j-1][1]

        assert is_cw(poly)

        all_ = {}

        for j in range(num):
            # dimlose taschen
            #if j in [3, 9]:#
            if j in [2, 3, 9, 10]:
                pass#continue

            if j != 5:
                pass#continue

            poly_ = deepcopy(poly)
            res, pts = rm_trivials(poly_, j)

            all_[j] = to_abs_path([ r['pt'] for r in res ])

            with open('data_files/special_%i.js' % i, 'w') as out:
                out.write('var polys = %s;' % json.dumps(all_))

        print to_abs_path(poly)
