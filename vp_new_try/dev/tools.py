# *-* coding: UTF-8 *-*

from copy import deepcopy
import math

def ld (a, b, c):
    return abs(a[0]*(b[1]-c[1])-a[1]*(b[0]-c[0])+(b[0]*c[1]-c[0]*b[1]))

def cross (a, b, c):
    # kreuzprodukt der vektoren ab und ac
    return (b[1]-a[1])*(c[0]-a[0])-(b[0]-a[0])*(c[1]-a[1])

def normalize (v):
    l = (v[0]*v[0]+v[1]*v[1])**-.5
    v[0] *= l
    v[1] *= l

    return 1/l

E = 1e-5

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
    return math.acos(r[0]*r_ab[0]+r[1]*r_ab[1]) > math.pi/2

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
    m11 = a[0]
    m12 = b[0]-a[0]
    m21 = a[1]
    m22 = b[1]-a[1]
    v1 = c[0]
    v2 = c[1]

    det = m11*m22-m12*m21

    return (m11*v2-v1*m21)/det

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
