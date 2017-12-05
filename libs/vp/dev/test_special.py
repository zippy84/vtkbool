#!/usr/bin/env python
# *-* coding: UTF-8 *-*

# Copyright 2012-2018 Ronald Römer
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import json
from copy import deepcopy
from tools import *
from rm_trivials import rm_internals, rm_trivials, add_internals

with open('special.json', 'r') as f:
    polys = json.load(f)['polys']

    for i, poly in enumerate(polys):
        if i != 0:
            pass#continue

        num = len(poly)

        for j in range(1, num):
            poly[j][0] += poly[j-1][0]
            poly[j][1] += poly[j-1][1]

        assert is_cw(poly)

        all_ = {}

        for j in range(num):
            poly_ = deepcopy(poly)
            res, pts = rm_trivials(poly_, j)

            res2 = add_internals(pts, res)

            all_[j] = to_abs_path([ r['pt'] for r in res2 ])

            with open('data_files/special_%i.js' % i, 'w') as out:
                out.write('var pts = {!r}; var polys = {!s};'.format(to_abs_path(poly), json.dumps(all_)))
