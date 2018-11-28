#!/usr/bin/env python2
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

path = 'M-0.188903,0.681407 L0.681407,0.188903 L0.188903,-0.681407 L-0.681407,-0.188903 L-0.371304,-0.248241 L-0.371553,-0.256357 L-0.372303,-0.280802 L-0.363765,-0.295019 L-0.346752,-0.323345 L-0.327831,-0.341416 L-0.298539,-0.369394 L-0.271473,-0.387517 L-0.235004,-0.411937 L-0.204777,-0.426163 L-0.165821,-0.444498 L-0.13843,-0.452005 L-0.101522,-0.46212 L-0.0820466,-0.46212 L-0.0518952,-0.46212 L-0.0427393,-0.456231 L-0.0244964,-0.444498 L-0.0242473,-0.436382 L-0.0234968,-0.411937 L-0.027061,-0.406002 L-0.0490485,-0.369393 L-0.0499043,-0.368576 L-0.0581019,-0.360747 L-0.0922776,-0.325232 L-0.0960884,-0.321272 L-0.15046,-0.283798 L-0.159115,-0.277833 L-0.22123,-0.24961 L-0.228806,-0.246168 L-0.282204,-0.233928 L-0.294278,-0.230619 L-0.295159,-0.230619 L-0.343905,-0.230619 L-0.347727,-0.233077 L-0.371304,-0.248241 L-0.681407,-0.188903 L-0.188903,0.681407 L-0.133981,0.14034 L-0.143137,0.134451 L-0.16138,0.122718 L-0.161629,0.114602 L-0.162379,0.0901566 L-0.158815,0.0842223 L-0.136828,0.0476135 L-0.135972,0.0467961 L-0.127774,0.0389666 L-0.0935986,0.0034521 L-0.0897879,-0.000507916 L-0.035416,-0.037982 L-0.0267608,-0.0439473 L0.0353537,-0.0721699 L0.0429295,-0.0756121 L0.0963279,-0.0878519 L0.108402,-0.0911609 L0.109283,-0.0911609 L0.158029,-0.0911609 L0.161851,-0.0887028 L0.185427,-0.0735389 L0.185677,-0.0654232 L0.186427,-0.0409778 L0.177888,-0.0267611 L0.160875,0.00156527 L0.141955,0.0196363 L0.112662,0.0476136 L0.0855968,0.0657369 L0.049128,0.0901567 L0.0189005,0.104383 L-0.020055,0.122718 L-0.0474465,0.130225 L-0.0843544,0.14034 L-0.10383,0.14034 L-0.133981,0.14034 Z'
f0 = 1

pts = map(lambda pt: map(lambda c: float(c), pt.split(',')), path[1:-2].split(' L'))

for i in range(1, len(pts)):
    pts[-i][0] = pts[-i][0]-pts[-(i+1)][0]
    pts[-i][1] = pts[-i][1]-pts[-(i+1)][1]

print 'path => m' + ' '.join([ ','.join(map(str, p)) for p in pts ])

polys = []

for f in [1, 10, 100, 1000]:
    pts2 = deepcopy(pts)

    for pt in pts2:
        pt[0] *= f0*f
        pt[1] *= f0*f

    polys.append(pts2)

with open('diagn.json', 'w') as f:
    json.dump({'polys': polys}, f)
