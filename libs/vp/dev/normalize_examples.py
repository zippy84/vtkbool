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
import os

size = [200, 300]
pad = 40

paths = []
data = []

for i, path in enumerate(paths):
    pts = [ [float(c) for c in coords.split(',')] for coords in path.split(' ') ]

    num = len(pts)

    for j in range(1, num):
        pts[j][0] += pts[j-1][0]
        pts[j][1] += pts[j-1][1]

    xs = [ pt[0] for pt in pts ]
    ys = [ pt[1] for pt in pts ]

    min_x = min(xs)
    min_y = min(ys)

    max_x = max(xs)
    max_y = max(ys)

    fac_x = (size[0]-pad)/(max_x-min_x)
    fac_y = (size[1]-pad)/(max_y-min_y)

    pts_ = zip([ (x-min_x)*fac_x+pad/2 for x in xs ],
        [ (y-min_y)*fac_y+pad/2 for y in ys ])

    pts_ = [ list(pt) for pt in pts_ ]

    for j in range(1, num):
        pts_[-j][0] = pts_[-j][0]-pts_[-(j+1)][0]
        pts_[-j][1] = pts_[-j][1]-pts_[-(j+1)][1]

    data.append(pts_)

print json.dumps({'paths': data})
