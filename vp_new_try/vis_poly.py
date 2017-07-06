#!/usr/bin/env python
# *-* coding: UTF-8 *-*

import math
from itertools import tee, izip
from copy import deepcopy
import json
from jinja2 import Environment, Template

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


def vis_poly (pts, ind=0):

    num = len(pts)

    pt_x = pts[ind]

    print 'M' + ' L'.join([ ','.join(map(str, p)) for p in pts ]) + 'Z'

    class Vert (object):
        def __init__(self, pt, ind = None, nxt = None):
            self.r = [pt[0]-pt_x[0], pt[1]-pt_x[1]]
            self.nxt = nxt
            self.pt = pt
            self.ind = ind

            normalize(self.r)

            self.phi = None

            if 'ref' in locals():
                self.phi = get_angle(ref, self.r)

    verts = []

    for i in range(num-1):
        _i = (ind+i+1)%num
        verts.append(Vert(pts[_i], ind=_i))

    ref = [verts[0].pt[0]-pt_x[0], verts[0].pt[1]-pt_x[1]]
    normalize(ref)

    for v in verts:
        v.phi = get_angle(ref, v.r)

        print v.phi*180/math.pi

    #exit()

    for i, v in enumerate(verts[:-1]):
        v.nxt = i+1

    print verts

    vp = [0, 1]

    t, u, v = (0, None, None)

    leftBags = []

    while 1:
        u = verts[t].nxt
        v = verts[u].nxt

        if not v:
            break

        print '>', u, v

        pt_u = verts[u].pt
        pt_v = verts[v].pt

        print 'orig', verts[u].ind, verts[v].ind

        if ld(pt_x, pt_u, pt_v) < 1e-2:
            print 'skipping'

            t = u
            continue

        c_a = cross(pt_x, pt_u, pt_v)
        c_b = cross(verts[t].pt, pt_u, pt_v)

        print 'c_a', c_a
        print 'c_b', c_b

        if c_a < 0:
            print 'vis'

            if vp[-1] != u:
                vp.append(u)

            vp.append(v)

            t = u
        else:
            if c_b > 0:

                w = v
                while 1:
                    a = verts[w]
                    b = verts[a.nxt]

                    pt_a = a.pt
                    pt_b = b.pt

                    print '>>', a.ind, b.ind

                    d = intersect(pt_x, verts[u].r, pt_a, pt_b)

                    if d and is_frontfaced(verts[u].r, pt_a, pt_b):
                        print 'edge (%d, %d)' % (a.ind, b.ind)

                        print d

                        if d['t2'] < E:
                            verts[u].nxt = w
                            vp.append(w)
                            t = u
                            leftBags.append({'f': u, 'g': w, 'phi': verts[u].phi})

                        else:

                            _v = Vert(d['s'], nxt = a.nxt)
                            verts.append(_v)

                            k = len(verts)-1

                            verts[u].nxt = k

                            vp.append(k)

                            t = u

                            leftBags.append({'f': u, 'g': k, 'phi': verts[u].phi})

                        break

                    else:
                        w = a.nxt

            elif c_b < 0:
                # schnitt mit leftBags?
                bag = None
                while leftBags and (leftBags[-1]['phi']-verts[v].phi) > -E:
                    bag = leftBags.pop()

                d = None

                if bag:
                    print 'bag', bag, verts[v].phi

                    d = intersect2(verts[bag['f']].pt, verts[bag['g']].pt, pt_u, pt_v)

                if d:
                    while vp and vp[-1] != bag['f']:
                        print 'popping_1', vp[-1]
                        vp.pop()

                    x = v

                    i = 0

                    while 1:
                        a = verts[x]
                        b = verts[a.nxt]

                        pt_a = a.pt
                        pt_b = b.pt

                        print '>>', a.ind, b.ind

                        d2 = {'s': verts[bag['f']].pt} if is_near(verts[bag['f']].pt, pt_v) else intersect2(verts[bag['f']].pt, d['s'], pt_b, pt_a)

                        if d2 and is_frontfaced(verts[bag['f']].r, pt_a, pt_b) \
                            and (i > 0 or cross(pt_a, pt_u, pt_b) < 0):

                            if is_near(verts[bag['f']].pt, d2['s']):
                                verts[bag['f']].nxt = a.nxt

                                vp.append(a.nxt)

                            else:
                                if d2['t2'] > 1-E:
                                    verts[bag['f']].nxt = x

                                    vp.append(x)

                                    leftBags.append({'f': bag['f'], 'g': x, 'phi': bag['phi']})

                                else:

                                    _v = Vert(d2['s'], nxt = a.nxt)
                                    verts.append(_v)

                                    k = len(verts)-1

                                    verts[bag['f']].nxt = k

                                    vp.append(k)

                                    leftBags.append({'f': bag['f'], 'g': k, 'phi': bag['phi']})

                            t = bag['f']

                            break

                        else:
                            x = a.nxt

                        i += 1

                else:
                    while vp:
                        a, b = vp[-2:]

                        print 'popping_2', vp[-1]

                        del vp[-1]

                        d = intersect(pt_x, verts[v].r, verts[a].pt, verts[b].pt)

                        if d:
                            if d['t2'] < E:
                                c = vp[-2]

                                if ld(pt_x, verts[a].pt, verts[c].pt) < 1e-3:
                                    del vp[-1]
                                    t = vp[-1]

                                else:
                                    t = a

                            else:
                                _v = Vert(d['s'])
                                verts.append(_v)

                                k = len(verts)-1

                                verts[a].nxt = k

                                vp.append(k)

                                t = k

                            break

                    pre_w = v

                    w = verts[v].nxt

                    if ld(pt_x, pt_v, verts[w].pt) < 1e-3:
                        pre_w = w
                        w = verts[w].nxt

                    print v, '->', pre_w

                    pt_w = verts[w].pt
                    pt_pre_w = verts[pre_w].pt

                    c_c = cross(pt_x, pt_v, pt_w)
                    c_d = cross(pt_v, pt_u, pt_w)

                    print 'c_c', c_c
                    print 'c_d', c_d

                    if c_c < 0:

                        if c_d < 0:
                            verts[vp[-1]].nxt = pre_w

                            vp.append(pre_w)

                        else:
                            x = w

                            while 1:
                                a = verts[x]
                                b = verts[a.nxt]

                                pt_a = a.pt
                                pt_b = b.pt

                                print '>>', a.ind, b.ind

                                d = intersect(pt_x, verts[v].r, pt_a, pt_b)

                                if d and not is_frontfaced(verts[v].r, pt_a, pt_b):
                                    print 'x'

                                    if d['t2'] < E:
                                        verts[vp[-1]].nxt = x
                                        vp.append(x)

                                    else:

                                        _v = Vert(d['s'], nxt = a.nxt)
                                        verts.append(_v)
                                        k = len(verts)-1

                                        verts[vp[-1]].nxt = k

                                        vp.append(k)

                                    break
                                else:
                                    x = a.nxt



                    else:
                        verts[vp[-1]].nxt = pre_w

                        vp.append(pre_w)

    #return 'M' + ' L'.join([ ','.join(map(str, verts[v_].pt)) for v_ in vp ]) + 'Z'

    pts_ = [pt_x] + [ verts[v_].pt for v_ in vp ]

    return pts_


if __name__ == '__main__':

    env = Environment()
    env.filters['to_path'] = to_path

    tmpl = env.from_string('''<?xml version="1.0" standalone="no"?>
        <svg xmlns="http://www.w3.org/2000/svg"
            version="1.1"
            width="800"
            height="{{ (rows*300)|int }}">
            <g>
                {% for d in data %}
                <path d="{{ d.vp|to_path(d.x, d.y) }}z" style="fill:#b3b3b3;"/>
                <path d="{{ d.pts|to_path(d.x, d.y) }}" style="stroke:black; stroke-width:1; fill:none;"/>
                <rect x="{{ d.x }}" y="{{ d.y }}" height="300" width="200" style="stroke:black; stroke-width:1; fill:none;"/>
                {% endfor %}
            </g>
        </svg>
        ''')

    with open('polys.json', 'r') as f:

        polys = json.load(f)['polys']

        data = []

        size = [200, 300]

        for i, poly in enumerate(polys):
            row, col = divmod(i, 4)

            pts = deepcopy(poly)

            num = len(pts)

            for j in range(1, num):
                pts[j][0] += pts[j-1][0]
                pts[j][1] += pts[j-1][1]

            vp = vis_poly(pts)

            data.append({ 'pts': pts,
                'vp': vp,
                'x': col*size[0],
                'y': row*size[1] })

        with open('output.svg', 'w') as out:
            out.write(tmpl.render(data=data, rows=math.ceil(len(data)/4.)))
