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

from copy import deepcopy
import math

E = 1e-5

def move (a, b, c):
    # damit a, b oder c nicht der ursprung ist
    # det(a,b,c) wäre dann auch 0

    z = [0, 0]

    if is_near(a, z) or is_near(b, z) or is_near(c, z):
        x = min(a[0], b[0], c[0])
        y = min(a[1], b[1], c[1])

        m = [1-x, 1-y]

        a[0] += m[0]; a[1] += m[1]
        b[0] += m[0]; b[1] += m[1]
        c[0] += m[0]; c[1] += m[1]

    assert not(is_near(a, z) or is_near(b, z) or is_near(c, z))

def ld (a, b, c):
    # muss unabh. von der skalierung sein

    global E

    v_a = [b[0]-a[0], b[1]-a[1]]
    v_b = [c[0]-a[0], c[1]-a[1]]

    l_a = normalize(v_a)
    l_b = normalize(v_b)

    f = 10/max(l_a, l_b);

    _a = [f*a[0], f*a[1]]
    _b = [f*b[0], f*b[1]]
    _c = [f*c[0], f*c[1]]

    move(_a, _b, _c)

    e = abs(v_a[0]*v_b[0]+v_a[1]*v_b[1])

    if e > E:
        if l_a > l_b:
            _c[0] = _a[0]+(5/e)*v_b[0]
            _c[1] = _a[1]+(5/e)*v_b[1]
        else:
            _b[0] = _a[0]+(5/e)*v_a[0]
            _b[1] = _a[1]+(5/e)*v_a[1]

    return abs(_a[0]*(_b[1]-_c[1])-_a[1]*(_b[0]-_c[0])+(_b[0]*_c[1]-_c[0]*_b[1]))

def cross (a, b, c):
    # kreuzprodukt der vektoren ab und ac
    return (b[1]-a[1])*(c[0]-a[0])-(b[0]-a[0])*(c[1]-a[1])

def normalize (v):
    try:
        l = (v[0]*v[0]+v[1]*v[1])**-.5
    except ZeroDivisionError:
        l = float('inf')

    v[0] *= l
    v[1] *= l

    return 1/l

def intersect (o, r, pA, pB):
    global E
    oB = [o[0]+r[0], o[1]+r[1]]

    m11 = oB[0]-o[0]
    m12 = pA[0]-pB[0]
    m21 = oB[1]-o[1]
    m22 = pA[1]-pB[1]
    v1 = pA[0]-o[0]
    v2 = pA[1]-o[1]

    det = m11*m22-m12*m21

    if abs(det) < E:
        return

    t1 = (v1*m22-m12*v2)/det
    t2 = (m11*v2-v1*m21)/det

    if t2 > -E and t2 < 1-E:
        return {'s': [pA[0]+t2*(pB[0]-pA[0]), pA[1]+t2*(pB[1]-pA[1])],
            't1': t1,
            't2': t2}

def intersect2 (oA, oB, pA, pB):
    global E

    m11 = oB[0]-oA[0]
    m12 = pA[0]-pB[0]
    m21 = oB[1]-oA[1]
    m22 = pA[1]-pB[1]
    v1 = pA[0]-oA[0]
    v2 = pA[1]-oA[1]

    det = m11*m22-m12*m21

    if abs(det) < E:
        return

    t1 = (v1*m22-m12*v2)/det
    t2 = (m11*v2-v1*m21)/det

    if t1 > -E and t1 < 1+E:
        if t2 > E and t2 < 1+E:
            return {'s': [pA[0]+t2*(pB[0]-pA[0]), pA[1]+t2*(pB[1]-pA[1])],
                't1': t1,
                't2': t2}

def is_frontfaced (r, a, b):
    r_ab = [a[1]-b[1], b[0]-a[0]] # um pi/2 gedreht
    normalize(r_ab)

    _x = r[0]*r_ab[0]+r[1]*r_ab[1]
    return _x < 0

def get_angle (a, b):
    ang = math.atan2(a[0]*b[1]-b[0]*a[1], a[0]*b[0]+a[1]*b[1])
    if ang < 0:
        ang += 2*math.pi
    return ang

def is_near (a, b):
    global E
    return abs(b[0]-a[0]) < E \
        and abs(b[1]-a[1]) < E

def from_path (path):
    pts = [ map(float, p.split(',')) for p in path[2:].split(' ')]

    num = len(pts)

    for i in range(1, num):
        pts[i][0] += pts[i-1][0]
        pts[i][1] += pts[i-1][1]

    return pts

def to_path (pts, x=0, y=0):
    pts_ = deepcopy(pts)

    for i in range(1, len(pts_)):
        pts_[-i][0] = pts_[-i][0]-pts_[-(i+1)][0]
        pts_[-i][1] = pts_[-i][1]-pts_[-(i+1)][1]

    pt = pts_[0]
    pt[0] += x
    pt[1] += y

    return 'm' + ' '.join([ ','.join(map(str, p)) for p in pts_ ])

def to_abs_path (pts):
    return 'M' + ' L'.join([ ','.join(map(str, p)) for p in pts ]) + 'Z'

def get_t (a, b, c):
    global E

    if is_near(b, c):
        return 1

    _a = a[:]
    _b = b[:]
    _c = c[:]

    move(_a, _b, _c)

    m11 = _a[0]
    m12 = _b[0]-_a[0]
    m21 = _a[1]
    m22 = _b[1]-_a[1]
    v1 = _c[0]
    v2 = _c[1]

    det = m11*m22-m12*m21

    assert abs(det) > E

    return (m11*v2-v1*m21)/det

def is_on_seg (a, b, c):
    global E

    if is_near(a, c) or is_near(b, c):
        return False

    if abs(a[0]-b[0]) > E:
        if (c[0] < a[0] and c[0] < b[0]) \
            or (c[0] > a[0] and c[0] > b[0]):
            return False
    else:
        if abs(a[0]-c[0]) > E:
            return False

    if abs(a[1]-b[1]) > E:
        if (c[1] < a[1] and c[1] < b[1]) \
            or (c[1] > a[1] and c[1] > b[1]):
            return False
    else:
        if abs(a[1]-c[1]) > E:
            return False

    return ld(a, b, c) < 1e-3

def is_cw (poly):
    # http://mathworld.wolfram.com/PolygonArea.html
    num = len(poly)
    s = 0

    for i in range(num):
        a = poly[i]
        b = poly[(i+1)%num]
        s += a[0]*b[1]-b[0]*a[1]

    return s < 0

def is_pip (poly, pt):
    for p in poly:
        if is_near(p, pt):
            return False

    num = len(poly)
    in_ = False

    for i in range(num):
        a = poly[i]
        b = poly[(i+1)%num]

        if (a[0] <= pt[0] or b[0] <= pt[0]) \
            and (a[1] < pt[1] and b[1] >= pt[1] \
                or b[1] < pt[1] and a[1] >= pt[1]):

            # schnittpunkt mit bounding box und strahlensatz
            if a[0]+(pt[1]-a[1])*(b[0]-a[0])/(b[1]-a[1]) < pt[0]:
                in_ = not(in_)

    return in_
