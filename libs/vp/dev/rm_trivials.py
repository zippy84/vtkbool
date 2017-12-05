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

def mark_internals (poly, skip):
    num = len(poly)

    ids = set()

    for i in range(num):
        j = (i+1)%num
        k = (i+2)%num

        if j != skip \
            and not is_near(poly[i]['pt'], poly[k]['pt']) \
            and ld(poly[i]['pt'], poly[j]['pt'], poly[k]['pt']) < 1e-3:

            poly[j]['marked'] = { 'a': i, 'b': k }
            ids.add(j)

    invalids = set()

    for id_ in ids:
        for i in range(num):
            if i not in ids \
                and is_near(poly[id_]['pt'], poly[i]['pt']):

                invalids.add(id_)

    for id_ in invalids:
        ids.remove(id_)
        del poly[id_]['marked']

    col = defaultdict(list)

    for id_ in ids:
        m = poly[id_]['marked']

        while m['a'] in ids:
            m['a'] = poly[m['a']]['marked']['a']

        while m['b'] in ids:
            m['b'] = poly[m['b']]['marked']['b']

        pt_a, pt_b = poly[m['a']]['pt'], \
            poly[m['b']]['pt']

        v = [pt_b[0]-pt_a[0], pt_b[1]-pt_a[1]]
        w = [poly[id_]['pt'][0]-pt_a[0], poly[id_]['pt'][1]-pt_a[1]]

        m['t'] = normalize(w)/normalize(v)

        col[(m['a'], m['b'])].append(id_)

    if len(col) > 1:
        for edge, nodes in col.items():
            if len(nodes) > 1:
                nodes2 = sorted(nodes, key=lambda n: poly[n]['marked']['t'])

                nodes2.insert(0, edge[0])
                nodes2.append(edge[1])

                q = [0, len(nodes2)-1]

                e = 0

                while True:
                    e = 0

                    l = len(q)
                    for i in range(l-1):
                        j = i+1

                        a = q[i]
                        b = q[j]

                        if b-a > 1:
                            c = a+(b-a)//2

                            pt_a = poly[nodes2[a]]['pt']
                            pt_b = poly[nodes2[b]]['pt']
                            pt_c = poly[nodes2[c]]['pt']

                            if not(ld(pt_a, pt_b, pt_c) < 1e-3):
                                e += 1
                                q.insert(j, c)
                                break

                    if e == 0:
                        break

                if len(q) > 2:
                    for r in q[1:-1]:
                        del poly[nodes2[r]]['marked']

    for p in poly:
        if 'marked' in p:
            print p['marked']

def align_pts (poly, ind):
    verts = []

    pX = poly[ind]['pt']

    for i, p in enumerate(poly):
        if i != ind \
            and 'marked' not in p \
            and not is_near(p['pt'], pX):

            r = [p['pt'][0]-pX[0], p['pt'][1]-pX[1]]

            d = normalize(r)

            verts.append({ 'i': i, 'r': r, 'd': d })

    verts.sort(key=itemgetter('d'))

    for v in verts:
        _v = [-v['r'][1], v['r'][0]]
        _d = _v[0]*pX[0]+_v[1]*pX[1]

        for w in verts:
            if v is w:
                continue

            d = -(poly[w['i']]['pt'][0]*_v[0]+poly[w['i']]['pt'][1]*_v[1]-_d)

            _l = ld(pX, poly[v['i']]['pt'], poly[w['i']]['pt'])

            if _l < 1e-3:
                poly[w['i']]['pt'][0] += d*_v[0]
                poly[w['i']]['pt'][1] += d*_v[1]

    for p in poly:
        if 'marked' in p:
            m = p['marked']

            pt_a, pt_b = poly[m['a']]['pt'], \
                poly[m['b']]['pt']

            v = [pt_b[0]-pt_a[0], pt_b[1]-pt_a[1]]

            p['pt'][0] = pt_a[0]+m['t']*v[0]
            p['pt'][1] = pt_a[1]+m['t']*v[1]

def rm_internals (poly):
    poly[:] = [ p for p in poly if 'marked' not in p ]

class TrivialRm:
    def __init__ (self, poly, ind):
        self._poly = poly
        self.poly = deepcopy(poly)
        self.verts = []
        self.orig = ind

        self.pX = poly[ind]

    def get_simplified (self):
        global E

        mark_internals(self.poly, self.orig)
        align_pts(self.poly, self.orig)

        self._poly[:] = [ { 'pt': p['pt'], 'idx': p['idx'] } for p in self.poly ]

        rm_internals(self.poly)

        num = len(self.poly)
        ind = next(( i for i, p in enumerate(self.poly) if p['idx'] == self.orig ))

        pA, pB = self.poly[(ind+1)%num], \
            self.poly[(ind+num-1)%num]

        rA = [pA['pt'][0]-self.pX['pt'][0], pA['pt'][1]-self.pX['pt'][1]]
        rB = [pB['pt'][0]-self.pX['pt'][0], pB['pt'][1]-self.pX['pt'][1]]

        srcs = defaultdict(list)

        for i, p in enumerate(self.poly):
            j = (i+1)%num

            self.verts.append(deepcopy(p))

            self.verts[-1].update({ 'src': None, 't': None })

            if i == ind:
                continue

            new_pts = []

            # entlang rA

            if ld(self.pX['pt'], pA['pt'], p['pt']) < 1e-3:
                t = get_t(self.pX['pt'], pA['pt'], p['pt'])

                if t > E:
                    self.verts[-1].update({ 'src': 'A', 't': t })

                    srcs[self.verts[-1]['idx']].append('A')
            else:
                try:
                    d = intersect(self.pX['pt'], rA, p['pt'], self.poly[j]['pt'])

                    if d['t1'] > E:
                        new_pts.append((d['t2'], { 'pt': d['s'], 'idx': None, 'src': 'A', 't': d['t1'] }))
                except:
                    pass

            # und entlang rB

            if ld(self.pX['pt'], pB['pt'], p['pt']) < 1e-3:
                t = get_t(self.pX['pt'], pB['pt'], p['pt'])

                if t > E:
                    self.verts[-1].update({ 'src': 'B', 't': t })

                    srcs[self.verts[-1]['idx']].append('B')
            else:
                try:
                    d = intersect(self.pX['pt'], rB, p['pt'], self.poly[j]['pt'])

                    if d['t1'] > E:
                        new_pts.append((d['t2'], { 'pt': d['s'], 'idx': None, 'src': 'B', 't': d['t1'] }))
                except:
                    pass

            new_pts.sort(key=itemgetter(0))

            self.verts.extend([ p[1] for p in new_pts ])

            if new_pts:
                print new_pts

        # ...

        for i, p in enumerate(self.verts):
            p['rm'] = False

            print i, p

        print 'srcs', srcs

        num = len(self.verts)
        ind = next(( i for i, p in enumerate(self.verts) if p['idx'] == self.orig ))

        nxt = self.verts[(ind+1)%num]

        print self.verts[ind]['idx'], nxt['idx']

        if nxt['idx'] in srcs \
            and len(srcs[nxt['idx']]) == 2:

            _rA = [-rA[0], -rA[1]]

            for i in range(num):
                j = (i+1)%num

                try:
                    d = intersect(self.pX['pt'], _rA, self.verts[i]['pt'], self.verts[j]['pt'])
                    if d['t1'] > E:
                        bnd = i

                        break
                except:
                    pass

            assert 'bnd' in locals()

            print bnd, self.verts[bnd]['idx']

            for i in range(num):
                j = (ind+i)%num
                p = self.verts[j]

                if p['src']:
                    p['src'] = 'A'

                if j == bnd:
                    break

            for i in range(num):
                j = (bnd+i)%num
                p = self.verts[j]

                if p['src'] and p['src'] == 'A':
                    p['rm'] = True

                if j == ind:
                    break

            self.verts[:] = [ p for p in self.verts if not p['rm'] ]

        def A():
            print 'A()'
            for i, p in enumerate(self.verts):
                p['i'] = i

            #print 'M' + ' L'.join([ ','.join(map(str, p['pt'])) for p in self.verts ]) + 'Z'

            good = [ p for p in self.verts if p['src'] == 'A' ]

            print 'good', good

            first = next(( i for i, g in enumerate(good) if g['idx'] == pA['idx'] ))
            _good = deque(good)
            _good.rotate(-first)

            print [ (g['i'], g['t']) for g in _good ]

            rot = [-rA[1], rA[0]]
            normalize(rot)

            d = rot[0]*self.pX['pt'][0]+rot[1]*self.pX['pt'][1]

            self.rm_pockets(_good, rot, d, 'A')

        A()

        def B():
            print 'B()'
            self.verts.reverse()

            for i, p in enumerate(self.verts):
                p['i'] = i

            #print 'M' + ' L'.join([ ','.join(map(str, p['pt'])) for p in self.verts ]) + 'Z'

            good = [ p for p in self.verts if p['src'] == 'B' ]

            print 'good', good

            first = next(( i for i, g in enumerate(good) if g['idx'] == pB['idx'] ))
            _good = deque(good)
            _good.rotate(-first)

            print [ (g['i'], g['t']) for g in _good ]

            rot = [rB[1], -rB[0]]
            normalize(rot)

            d = rot[0]*self.pX['pt'][0]+rot[1]*self.pX['pt'][1]

            self.rm_pockets(_good, rot, d, 'B')

        B()

        # ergebnis

        res = deque([ p for p in self.verts if not p['rm'] ])

        ni = next((i for i, p in enumerate(res) if p['idx'] == self.orig))

        res.rotate(-ni)

        res = [ { 'pt': p['pt'], 'idx': p['idx'] } for p in res ]

        mark_internals(res, 0)
        rm_internals(res)

        assert res[0]['idx'] == self.orig

        return res

    def get_pocket (self, pair):
        num = len(self.verts)

        id_ = pair['i']

        pocket = [id_]

        while id_ != pair['j']:
            id_ = (id_+1)%num
            pocket.append(id_)

        print 'pocket', pocket, \
            [ self.verts[id_]['idx'] for id_ in pocket ]

        return pocket

    def assign_side (self, pair, src):
        pair['pocket'] = self.get_pocket(pair)

        if 'side' in pair:
            return

        if len(pair['pocket']) > 2:
            if src == 'B':
                pair['pocket'].reverse()

            p_poly = [ self.verts[id_]['pt'] for id_ in pair['pocket'] ]

            pair['side'] = Side.IN if not is_cw(p_poly) else Side.OUT

    def has_area (self, pocket):
        area = True
        l = len(pocket)

        if l%2 == 1:
            for i in range((l-1)/2):
                area = not(is_near(self.verts[pocket[i]]['pt'], self.verts[pocket[l-i-1]]['pt']))

        return area

    def rm_pockets (self, _good, rot, d, src):
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
                pocket = self.get_pocket(p)

                p_poly = [ self.verts[id_]['pt'] for id_ in pocket ]

                if is_pip(p_poly, self.pX['pt']):

                    if (src == 'A' and not is_cw(p_poly)) \
                        or (src == 'B' and is_cw(p_poly)):

                        break

                    print 'erasing from', i

                    del pairs[i:]

                    break

            for p in pairs:
                if p['dir'] == Dir.UNDEFINED:
                    pocket = self.get_pocket(p)

                    ss = set()

                    for id_ in pocket[1:-1]:
                        e = (self.verts[id_]['pt'][0]*rot[0]+self.verts[id_]['pt'][1]*rot[1])-d
                        if abs(e) > E:
                            ss.add(int(e/abs(e)))

                    assert len(ss) == 1

                    if self.has_area(pocket):

                        p_poly = [ self.verts[id_]['pt'] for id_ in pocket ]

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

            print grps

            if grps[0][0] == Dir.BACKWARD:
                if len(grps) > 1:
                    pairs[grps[1][1][0]]['i'] = pairs[grps[0][1][0]]['i']
                    del grps[0]

            for p in pairs:
                if p['dir'] == Dir.FORWARD:
                    self.assign_side(p, src)

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

                            a, b, c = self.verts[p_['i']]['t'], \
                                self.verts[p_['j']]['t'], \
                                self.verts[last_pair['j']]['t']

                            # liegt c zw (a, b)?

                            q = c-a > E and b-c > E and p_['side'] == Side.OUT

                            print j, (a, b, c), q

                            if q:
                                print '(1) new pair', (p_['i'], last_pair['j'])

                                del grps[i-1][1][-j-1:]
                                grps[i-1][1].append(add_pair({ 'i': p_['i'], 'j': last_pair['j'] }))
                                self.assign_side(pairs[grps[i-1][1][-1]], src)

                                _c = False

                                break

                        if _c:
                            if (i+1) in ids:
                                # es gibt einen danach

                                print '(2) new pair', \
                                    pairs[grps[i-1][1][-1]]['j'], \
                                    pairs[grps[i+1][1][0]]['j']

                                grps[i-1][1].append(add_pair({ 'i': pairs[grps[i-1][1][-1]]['j'], 'j': pairs[grps[i+1][1][0]]['j'] }))

                                self.assign_side(pairs[grps[i-1][1][-1]], src)

                                del grps[i+1][1][0]

                            else:
                                for j, p in enumerate(grps[i-1][1][::-1]):
                                    p_ = pairs[p]

                                    a = self.verts[p_['j']]['t']
                                    b = self.verts[last_pair['j']]['t']

                                    if abs(a-b) < E:
                                        print '(3) new pair', \
                                            p_['i'], \
                                            last_pair['j']

                                        del grps[i-1][1][-j-1:]

                                        grps[i-1][1].append(add_pair({ 'i': p_['i'], 'j': last_pair['j'] }))
                                        self.assign_side(pairs[grps[i-1][1][-1]], src)

                                        break

                        del grps[i]

                assert len(grps) == 1

                print 'grps[0][1]', grps[0][1]

                new_pairs.extend([ pairs[g] for g in grps[0][1] ])

            print 'new_pairs', new_pairs

            for p in new_pairs:
                assert 'pocket' in p

                if 'side' in p and p['side'] == Side.OUT:
                    for i in p['pocket'][1:-1]:
                        self.verts[i]['rm'] = True

def add_internals (orig, poly):
    global E

    num_a, num_b = len(poly), \
        len(orig)

    res = []

    for i in range(num_a):
        j = (i+1)%num_a

        p_a, p_b = poly[i], \
            poly[j]

        u = [p_b['pt'][0]-p_a['pt'][0], p_b['pt'][1]-p_a['pt'][1]]
        r = normalize(u)

        res.append(p_a)

        if is_near(p_a['pt'], p_b['pt']):
            assert p_a['idx'] is not None \
                and p_b['idx'] is not None

            k = p_b['idx']

            while k != p_a['idx']:
                k = (k+1)%num_b
                res.append(orig[k])

            del res[-1]

            continue

        verts = []

        for p in orig:
            if is_on_seg(p_a['pt'], p_b['pt'], p['pt']):
                v = [p['pt'][0]-p_a['pt'][0],
                    p['pt'][1]-p_a['pt'][1]]
                t = normalize(v)/r

                _v = { 'valid': True, 't': t }
                _v.update(deepcopy(p))

                verts.append(_v)

        verts.sort(key=itemgetter('t'))

        if len(verts) > 1:
            ids = { p_a['idx']: 0, p_b['idx']: 1 }

            for v in verts:
                ids[v['idx']] = v['t']

            for k in range(len(verts)-1):
                l = k+1
                if verts[l]['t']-verts[k]['t'] < E:
                    id_a, id_b = verts[k]['idx'], \
                        verts[l]['idx']

                    a = (id_a+1)%num_b
                    b = (id_a+num_b-1)%num_b

                    c_a = (a in ids and ids[a] < ids[id_a]) \
                        +(b in ids and ids[b] > ids[id_a])

                    a = (id_b+1)%num_b
                    b = (id_b+num_b-1)%num_b

                    c_b = (a in ids and ids[a] < ids[id_b]) \
                        +(b in ids and ids[b] > ids[id_b])

                    assert c_a > 0 or c_b > 0

                    if c_a > c_b:
                        verts[l]['valid'] = False
                    else:
                        verts[k]['valid'] = False

        for v in verts:
            if v['valid']:
                res.append({ 'pt': v['pt'], 'idx': v['idx'] })

    poly[:] = res
