#!/usr/bin/env python
# *-* coding: UTF-8 *-*

from itertools import tee, izip, groupby, cycle
from operator import itemgetter
from collections import deque
from enum import Enum
import json
from collections import defaultdict
from copy import deepcopy

import os

from tools import *

class Dir(Enum):
    FORWARD = 1
    BACKWARD = 2
    UNDEFINED = 3

class Side(Enum):
    IN = 1
    OUT = 2

E = 1e-5

def align_pts (poly, ind):
    verts = []

    pt_x = poly[ind]['pt']

    for i, p in enumerate(poly):
        if i != ind \
            and not is_near(p['pt'], pt_x) \
            and 'int' not in p:

            r = [p['pt'][0]-pt_x[0], p['pt'][1]-pt_x[1]]

            d = normalize(r)
            phi = get_angle([1, 0], r)

            verts.append({ 'i': i, 'r': r, 'd': d, 'phi': phi })

    phis = defaultdict(list)

    for i, v in enumerate(verts):
        phis['%.3f' % v['phi']].append(i)

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

                poly[v['i']]['pt'] = [pt_x[0]+r[0], pt_x[1]+r[1]]


    for p in poly:
        if 'int' in p:
            id_a, id_b = p['int']['edge']

            pt_a = poly[id_a]['pt']
            pt_b = poly[id_b]['pt']

            v = [pt_b[0]-pt_a[0], pt_b[1]-pt_a[1]]

            new_pt = [pt_a[0]+p['int']['t']*v[0],
                pt_a[1]+p['int']['t']*v[1]]

            p['pt'] = new_pt


def rm_internals (poly, skip):
    num = len(poly)

    for i in range(num):
        j = (i+1)%num
        k = (i+2)%num

        if j != skip \
            and not is_near(poly[i]['pt'], poly[k]['pt']) \
            and ld(poly[i]['pt'], poly[j]['pt'], poly[k]['pt']) < 1e-2:

            poly[j]['int'] = {}

    ids = [ i for i, p in enumerate(poly) if 'int' not in p ]
    #print ids

    num_ = len(ids)

    for i in range(num_):
        j = (i+1)%num_

        id_a, id_b = ids[i], ids[j]

        pt_a = poly[id_a]['pt']
        pt_b = poly[id_b]['pt']

        v = [pt_b[0]-pt_a[0], pt_b[1]-pt_a[1]]
        l = normalize(v)

        betw = []

        id_c = id_a

        while id_c != id_b:
            id_c = (id_c+1)%num
            betw.append(id_c)

        if len(betw) > 1:
            del betw[-1]

            #print id_a, id_b, betw

            for id_ in betw:
                pt_c = poly[id_]['pt']

                v_ = [pt_c[0]-pt_a[0], pt_c[1]-pt_a[1]]

                l_ = normalize(v_)

                print id_, l_/l

                poly[id_]['int'].update({ 't': l_/l, 'edge': (id_a, id_b) })

    for i, p in enumerate(poly):
        if 'int' in p:
            for j, q in enumerate(poly):
                if i == j:
                    continue

                if is_near(p['pt'], q['pt']):
                    del p['int']

def rm_trivials (pts, ind):
    global E

    _pts = [ { 'pt': pt, 'idx': i } for i, pt in enumerate(pts) ]

    rm_internals(_pts, ind)
    align_pts(_pts, ind)

    _poly = [ deepcopy(p) for p in _pts if 'int' not in p ]

    # cleanup
    for p in _pts:
        try:
            del p['int']
        except KeyError:
            pass

    _ind = next(( i for i, p in enumerate(_poly) if p['idx'] == ind ))

    print '_ind', _ind

    _num = len(_poly)

    pA, pB = _poly[(_ind+1)%_num], \
        _poly[(_ind+_num-1)%_num]

    pX = _poly[_ind]

    rA = [pA['pt'][0]-pX['pt'][0], pA['pt'][1]-pX['pt'][1]]
    rB = [pB['pt'][0]-pX['pt'][0], pB['pt'][1]-pX['pt'][1]]

    _test = defaultdict(list)

    poly = []

    for i in range(_num):
        poly.append(_poly[i])

        poly[-1].update({ 'src': None, 't': None })

        if i == _ind:
            continue

        j = (i+1)%_num

        new_pts = []
        # entlang rA

        if ld(pX['pt'], pA['pt'], _poly[i]['pt']) < 1e-2:
            t = get_t(pX['pt'], pA['pt'], _poly[i]['pt'])

            if t > E:
                poly[-1].update({ 'src': 'A', 't': t })

                _test[poly[-1]['idx']].append('A')
        else:
            try:
                d = intersect(pX['pt'], rA, _poly[i]['pt'], _poly[j]['pt'])

                if d['t1'] > E:
                    new_pts.append((d['t2'], { 'pt': d['s'], 'idx': None, 'src': 'A', 't': d['t1'] }))
            except:
                pass

        # und entlang rB

        if ld(pX['pt'], pB['pt'], _poly[i]['pt']) < 1e-2:
            t = get_t(pX['pt'], pB['pt'], _poly[i]['pt'])

            if t > E:
                poly[-1].update({ 'src': 'B', 't': t })

                _test[poly[-1]['idx']].append('B')
        else:
            try:
                d = intersect(pX['pt'], rB, _poly[i]['pt'], _poly[j]['pt'])

                if d['t1'] > E:
                    new_pts.append((d['t2'], { 'pt': d['s'], 'idx': None, 'src': 'B', 't': d['t1'] }))
            except:
                pass

        new_pts.sort(key=itemgetter(0))

        poly.extend([ p[1] for p in new_pts ])

        print new_pts

    for i, p in enumerate(poly):
        p['rm'] = False

        print i, p

    print '_test', _test

    _num = len(poly)
    _ind = next(( i for i, p in enumerate(poly) if p['idx'] == ind ))

    _nxt = poly[(_ind+1)%_num]

    print poly[_ind]['idx'], _nxt['idx']

    if _nxt['idx'] in _test \
        and len(_test[_nxt['idx']]) == 2:

        _rA = [-rA[0], -rA[1]]

        for i in range(_num):
            j = (i+1)%_num

            try:
                d = intersect(pX['pt'], _rA, poly[i]['pt'], poly[j]['pt'])
                if d['t1'] > E:
                    _bnd = i

                    break
            except:
                pass

        assert '_bnd' in locals()

        print _bnd, poly[_bnd]['idx']

        for i in range(_num):
            j = (_ind+i)%_num
            p = poly[j]

            if p['src']:
                p['src'] = 'A'

            if j == _bnd:
                break

        for i in range(_num):
            j = (_bnd+i)%_num
            p = poly[j]

            if p['src'] and p['src'] == 'A':
                p['rm'] = True

            if j == _ind:
                break

        poly = [ p for p in poly if not p['rm'] ]

    num_ = len(poly)

    def get_pocket (pair):
        id_ = pair['i']

        pocket = [id_]

        while id_ != pair['j']:
            id_ = (id_+1)%num_
            pocket.append(id_)

        print 'pocket', pocket, \
            [ poly[id_]['idx'] for id_ in pocket ]

        return pocket

    def assign_side (pair, src):
        pair['pocket'] = get_pocket(pair)

        if 'side' in pair:
            return

        if len(pair['pocket']) > 2:
            if src == 'B':
                pair['pocket'].reverse()

            p_poly = [ poly[id_]['pt'] for id_ in pair['pocket'] ]

            pair['side'] = Side.IN if not is_cw(p_poly) else Side.OUT

    def has_area (pocket):
        area = True
        l = len(pocket)

        if l%2 == 1:
            for i in range((l-1)/2):
                area = not(is_near(poly[pocket[i]]['pt'], poly[pocket[l-i-1]]['pt']))

        return area

    def remove_pockets (_good, rot, d, src):
        print rot, d

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

        print 'pairs', pairs

        if pairs:
            for i, p in enumerate(pairs):
                pocket = get_pocket(p)

                p_poly = [ poly[id_]['pt'] for id_ in pocket ]

                if len(set([ is_near(pX['pt'], pt) for pt in p_poly ])) == 1 \
                    and is_pip(p_poly, pX['pt']):

                    del pairs[i:]
                    print 'del..'
                    break

            for p in pairs:
                if p['dir'] == Dir.UNDEFINED:
                    pocket = get_pocket(p)

                    ss = set()

                    for id_ in pocket[1:-1]:
                        e = (poly[id_]['pt'][0]*rot[0]+poly[id_]['pt'][1]*rot[1])-d
                        if abs(e) > E:
                            print e
                            ss.add(int(e/abs(e)))

                    assert len(ss) == 1

                    if has_area(pocket):

                        p_poly = [ poly[id_]['pt'] for id_ in pocket ]

                        if src == 'B':
                            p_poly.reverse()

                        #print ss.pop() == 1, is_cw(p_poly),

                        is_backw = (ss.pop() == 1)^is_cw(p_poly)

                        if is_backw:
                            p['dir'] = Dir.BACKWARD
                        else:
                            p['dir'] = Dir.FORWARD
                    else:
                        # das kann gar nicht anders sein

                        if ss.pop() == 1:
                            p['dir'] = Dir.BACKWARD
                        else:
                            p['dir'] = Dir.FORWARD
                            p['side'] = Side.IN

        assert len([ p for p in pairs if p['dir'] == Dir.UNDEFINED ]) == 0

        def add_pair(pair):
            pairs.append(pair)
            return len(pairs)-1

        if pairs:
            new_pairs = []

            grps = []
            o = 0
            for k, g in groupby(pairs, key=itemgetter('dir')):
                n = len(list(g))
                grps.append((k, [ o+i for i in range(n) ]))

                o += n

            print ind, grps

            # TODO: kann sowas noch vorkommen?
            if grps[0][0] == Dir.BACKWARD:
                if len(grps) > 1:
                    pairs[grps[1][1][0]]['i'] = pairs[grps[0][1][0]]['i']
                    del grps[0]

            for p in pairs:
                if p['dir'] == Dir.FORWARD:
                    assign_side(p, src)

            if len(grps) > 0:
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
                                assign_side(pairs[grps[i-1][1][-1]], src)

                                _c = False

                                break

                        if _c and (i+1) in ids:
                            # es gibt einen danach

                            print '(2) new pair', \
                                pairs[grps[i-1][1][-1]]['j'], \
                                pairs[grps[i+1][1][0]]['j']

                            grps[i-1][1].append(add_pair({ 'i': pairs[grps[i-1][1][-1]]['j'], 'j': pairs[grps[i+1][1][0]]['j'] }))

                            assign_side(pairs[grps[i-1][1][-1]], src)

                            del grps[i+1][1][0]

                        del grps[i]

                assert len(grps) == 1

                print 'grps[0][1]', grps[0][1]

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

        #print 'M' + ' L'.join([ ','.join(map(str, p['pt'])) for p in poly ]) + 'Z'

        good = [ p for p in poly if p['src'] == 'A' ]

        print 'good', good

        first = next(( i for i, g in enumerate(good) if g['idx'] == pA['idx'] ))
        _good = deque(good)
        _good.rotate(-first)

        print [ (g['i'], g['t']) for g in _good ]

        rot = [-rA[1], rA[0]]
        normalize(rot)

        d = rot[0]*pX['pt'][0]+rot[1]*pX['pt'][1]

        remove_pockets(_good, rot, d, 'A')

    A()

    def B():
        print 'B()'
        poly.reverse()

        for i, p in enumerate(poly):
            p['i'] = i

        #print 'M' + ' L'.join([ ','.join(map(str, p['pt'])) for p in poly ]) + 'Z'

        good = [ p for p in poly if p['src'] == 'B' ]

        print 'good', good

        first = next(( i for i, g in enumerate(good) if g['idx'] == pB['idx'] ))
        _good = deque(good)
        _good.rotate(-first)

        print [ (g['i'], g['t']) for g in _good ]

        rot = [rB[1], -rB[0]]
        normalize(rot)

        d = rot[0]*pX['pt'][0]+rot[1]*pX['pt'][1]

        remove_pockets(_good, rot, d, 'B')

    B()

    res = deque([ p for p in poly if not p['rm'] ])

    ni = next((i for i, p in enumerate(res) if p['idx'] == ind))

    res.rotate(-ni)

    rm_internals(res, 0)

    _res = [ p for p in res if 'int' not in p ]

    assert _res[0]['idx'] == ind

    return [ { 'pt': p['pt'], 'idx': p['idx'] } for p in _res ], _pts


def add_internals (pts, poly):
    num = len(poly)

    _pts = deque(pts)
    _pts.reverse()

    first = next(( i for i, p in enumerate(_pts) if p['idx'] == poly[0]['idx'] ))

    _pts.rotate(-first)

    _pts = list(_pts)
    _num = len(_pts)

    for i, p in enumerate(poly):
        found = []

        for j, _p in enumerate(_pts):
            if is_near(_p['pt'], p['pt']):
                found.append(j)

        print i, found

    return poly
