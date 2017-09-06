#!/usr/bin/env python
# *-* coding: UTF-8 *-*

from itertools import tee, izip, groupby, cycle
from operator import itemgetter
from collections import deque
from enum import Enum
import json
from collections import defaultdict

from tools import *

class Dir(Enum):
    FORWARD = 1
    BACKWARD = 2
    UNDEFINED = 3

class Side(Enum):
    IN = 1
    OUT = 2

E = 1e-5

def align_pts (pts, ind):
    verts = []

    pt_x = pts[ind]

    for i, pt in enumerate(pts):
        if i != ind \
            and not is_near(pt, pt_x):

            r = [pt[0]-pt_x[0], pt[1]-pt_x[1]]

            d = normalize(r)
            phi = get_angle([1, 0], r)

            verts.append({ 'i': i, 'r': r, 'd': d, 'phi': phi })

    phis = defaultdict(list)

    for i, v in enumerate(verts):
        phis['%.4f' % v['phi']].append(i)

    print 'phis', phis

    for phi, ids in phis.items():
        if len(ids) > 1:
            ids.sort(key=lambda id_: verts[id_]['d'])

            for id_ in ids:
                v = verts[id_]

                phi_ = float(phi)-v['phi']

                # dreht den ortsvektor um den kleinen winkel

                w = [v['d']*v['r'][0], v['d']*v['r'][1]]
                r = [w[0]-w[1]*phi_, w[0]*phi_+w[1]]

                pts[v['i']] = [pt_x[0]+r[0], pt_x[1]+r[1]]


def rm_internals (pts, skip):
    print 'pts', pts
    num = len(pts)

    internal = []

    for i in range(num):
        j = (i+1)%num
        k = (i+2)%num

        if ld(pts[i]['pt'], pts[j]['pt'], pts[k]['pt']) < 1e-2:
            internal.append(j)

    try:
        internal.remove(skip)
    except ValueError:
        pass

    print internal

    ids = sorted(internal)

    for i in ids[::-1]:
        del pts[i]

    print pts


def rm_trivials (pts, ind):
    global E

    align_pts(pts, ind)

    num = len(pts)

    iA, iB = (ind+1)%num, \
    (ind+num-1)%num

    print iA, iB

    rA = [pts[iA][0]-pts[ind][0], pts[iA][1]-pts[ind][1]]
    rB = [pts[iB][0]-pts[ind][0], pts[iB][1]-pts[ind][1]]

    ind_ = None

    poly = []

    for i in range(num):
        poly.append({ 'pt': pts[i], 'idx': i, 'src': None, 't': None })

        if i == ind:
            ind_ = len(poly)
            continue

        j = (i+1)%num

        new_pts = []
        # entlang rA

        if ld(pts[ind], pts[iA], pts[i]) < 1e-2:
            t = get_t(pts[ind], pts[iA], pts[i])

            if t > E:
                poly[-1].update({ 'src': 'A', 't': t })
        else:
            try:
                d = intersect(pts[ind], rA, pts[i], pts[j])

                if d['t1'] > E:
                    new_pts.append((d['t2'], { 'pt': d['s'], 'idx': None, 'src': 'A', 't': d['t1'] }))
            except:
                pass

        # und entlang rB

        if ld(pts[ind], pts[iB], pts[i]) < 1e-2:
            t = get_t(pts[ind], pts[iB], pts[i])

            if t > E:
                poly[-1].update({ 'src': 'B', 't': t })
        else:
            try:
                d = intersect(pts[ind], rB, pts[i], pts[j])

                if d['t1'] > E:
                    new_pts.append((d['t2'], { 'pt': d['s'], 'idx': None, 'src': 'B', 't': d['t1'] }))
            except:
                pass

        new_pts.sort(key=itemgetter(0))

        poly.extend([ p[1] for p in new_pts ])

        print new_pts

    for i, p in enumerate(poly):
        p['rm'] = False

    #for i, p in enumerate(poly):
    #    print i, p

    #print 'M' + ' L'.join([ ','.join(map(str, p['pt'])) for p in poly ]) + 'Z'

    #print 'A:', [ (i, p) for i, p in enumerate(poly) if p['src'] == 'A' ]
    #print 'B:', [ (i, p) for i, p in enumerate(poly) if p['src'] == 'B' ]

    num_ = len(poly)

    def remove_pockets (_good, rot, d):
        print rot, d

        def add_pocket (pair):
            k = pair['i']
            pair['pocket'] = [k]

            while True:
                k = (k+1)%num_
                pair['pocket'].append(k)
                if k == pair['j']:
                    break

            if len(pair['pocket']) > 2:
                for p_id in pair['pocket'][1:-1]:
                    pt = poly[p_id]['pt']
                    e = (pt[0]*rot[0]+pt[1]*rot[1])-d

                    if abs(e) > E:
                        pair['side'] = Side.IN if e < 0 else Side.OUT
                        break

        pairs = []

        a, b = tee(_good)
        next(b, None)

        pairs = []

        for pair in izip(a, b):
            p_a, p_b = pair

            pairs.append({ 'i': p_a['i'], 'j': p_b['i'],
                'dir': Dir.FORWARD if p_b['t'] > p_a['t'] else Dir.BACKWARD })

            if abs(p_b['t']-p_a['t']) < E:
                pairs[-1]['dir'] = Dir.UNDEFINED

        #print 'pairs', pairs

        def add_pair(pair):
            pairs.append(pair)
            return len(pairs)-1

        if pairs:
            new_pairs = []

            for p in pairs:
                if p['dir'] == Dir.FORWARD:
                    add_pocket(p)

            grps = []
            o = 0
            for k, g in groupby(pairs, key=itemgetter('dir')):
                n = len(list(g))
                grps.append((k, [ o+i for i in range(n) ]))

                o += n

            print ind, grps

            for i, g in enumerate(grps):
                if g[0] == Dir.UNDEFINED:
                    if i > 0 and grps[i-1][0] == Dir.BACKWARD:
                        grps[i] = (Dir.BACKWARD, g[1])
                    else:
                        add_pocket(pairs[g[1][0]])
                        grps[i] = (Dir.FORWARD, g[1])


            if len(grps) > 0 and grps[0][0] == Dir.FORWARD:
                while True:
                    # fasst erneut zusammen
                    grps = [ (f, sum([ h[1] for h in g ], []) ) for f, g in groupby(grps, key=itemgetter(0)) ]

                    ids = range(len(grps))

                    print '>>', grps

                    try:
                        i, g = next(( (i, g) for i, g in izip(ids, grps) if g[0] == Dir.BACKWARD ))
                        print 'i', i
                    except StopIteration:
                        break
                    else:

                        _c = True

                        last_pair = pairs[g[1][-1]]

                        for j, p in enumerate(grps[i-1][1][::-1]):
                            p_ = pairs[p]

                            a, b, c = poly[p_['i']]['t'], \
                                poly[p_['j']]['t'], \
                                poly[last_pair['j']]['t']

                            # liegt c zw (a, b)?

                            q = c-a > E and b-c > E and p_['side'] == Side.OUT

                            print j, (a, b, c), q

                            if q:
                                print '(1) new pair', (p_['i'], last_pair['j'])

                                del grps[i-1][1][-j-1:]
                                grps[i-1][1].append(add_pair({ 'i': p_['i'], 'j': last_pair['j'] }))
                                add_pocket(pairs[grps[i-1][1][-1]])

                                _c = False

                                break

                        if _c and (i+1) in ids:
                            # es gibt einen danach

                            print '(2) new pair', \
                                pairs[grps[i-1][1][-1]]['j'], \
                                pairs[grps[i+1][1][0]]['j']

                            grps[i-1][1].append(add_pair({ 'i': pairs[grps[i-1][1][-1]]['j'], 'j': pairs[grps[i+1][1][0]]['j'] }))

                            add_pocket(pairs[grps[i-1][1][-1]])

                            del grps[i+1][1][0]

                        del grps[i]

                assert len(grps) == 1

                print grps[0][1]

                new_pairs.extend([ pairs[g] for g in grps[0][1] ])

            print 'new_pairs', new_pairs

            for p in new_pairs:
                assert 'pocket' in p

                if 'side' in p and p['side'] == Side.OUT:
                    for i in p['pocket'][1:-1]:
                        poly[i]['rm'] = True

    def A():
        print 'A()'
        for i, p in enumerate(poly):
            p['i'] = i

        print 'M' + ' L'.join([ ','.join(map(str, p['pt'])) for p in poly ]) + 'Z'

        good = [ p for p in poly if p['src'] == 'A' ]

        first = next(( i for i, g in enumerate(good) if g['idx'] == iA ))
        _good = deque(good)
        _good.rotate(-first)

        print [ (g['i'], g['t']) for g in _good ]

        rot = [-rA[1], rA[0]]
        normalize(rot)

        d = rot[0]*pts[ind][0]+rot[1]*pts[ind][1]

        remove_pockets(_good, rot, d)

    A()

    def B():
        print 'B()'
        poly.reverse()

        for i, p in enumerate(poly):
            p['i'] = i

        print 'M' + ' L'.join([ ','.join(map(str, p['pt'])) for p in poly ]) + 'Z'

        good = [ p for p in poly if p['src'] == 'B' ]

        first = next(( i for i, g in enumerate(good) if g['idx'] == iB ))
        _good = deque(good)
        _good.rotate(-first)

        print [ (g['i'], g['t']) for g in _good ]

        rot = [rB[1], -rB[0]]
        normalize(rot)

        d = rot[0]*pts[ind][0]+rot[1]*pts[ind][1]

        remove_pockets(_good, rot, d)

    B()

    res = deque([ p for p in poly if not p['rm'] ])

    ni = next((i for i, p in enumerate(res) if p['idx'] == ind))

    res.rotate(-ni)

    rm_internals(res, 0)

    assert res[0]['idx'] == ind

    return [ { 'pt': p['pt'], 'idx': p['idx'] } for p in res ]


def add_internals (pts, poly):
    num = len(poly)

    res = []

    for i in range(num):
        p_a = poly[i]
        p_b = poly[(i+1)%num]

        res.append(p_a)

        for j, pt in enumerate(pts):
            if is_on_seg(p_a['pt'], p_b['pt'], pt):
                res.append({ 'pt': pt, 'idx': j })

    return res
