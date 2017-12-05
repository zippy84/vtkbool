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

import sys
import json

from collections import deque, defaultdict
from copy import deepcopy

sys.path.append('../../vp/dev')

from vis_poly import vis_poly_wrapper
from tools import cross, ld, is_cw, to_abs_path, to_path, is_near, normalize

from rm_trivials import mark_internals, rm_internals, add_internals

from jinja2 import Environment, Template

class SubP:
    def __init__ (self):
        self.w = 99999
        self.S = []

        self.S_head = []
        self.S_tail = []
    def add_pair (self, pair, w):
        # bei gleichem w nur die narrowest pairs beibehalten

        if w > self.w:
            return

        if w < self.w:
            del self.S[:]
            del self.S_tail[:]

        if self.S:
            if pair[0] > self.S[0][0]:
                while self.S and self.S[0][1] >= pair[1]:
                    # das neue pair ist im intervall S[0] enthalten
                    del self.S[0]

        # durch einfügen an erster stelle ist S cw-ordered (gegenüber dem polygon)
        self.S.insert(0, pair)
        self.w = w

        del self.S_head[:]

    def restore_S (self):
        self.S[:0] = self.S_head
        self.S.extend(reversed(self.S_tail))

        del self.S_head[:]
        del self.S_tail[:]

def decompose (pts):
    poly = [ { 'pt': pt, 'idx': i } for i, pt in enumerate(pts) ]

    original = deepcopy(poly)

    scale(poly)
    simple_rm_internals(poly)

    poly = deque(poly)

    num = len(poly)

    def is_refl (a, b, c):
        return is_near(poly[b]['pt'], poly[c]['pt']) \
            or (not ld(poly[a]['pt'], poly[b]['pt'], poly[c]['pt']) < 1e-2 \
                and cross(poly[a]['pt'], poly[b]['pt'], poly[c]['pt']) < 0)

    for i in range(num):
        j = (i+1)%num
        k = (i+num-1)%num

        poly[i]['refl'] = is_refl(i, j, k)

    first = next(( i for i, p in enumerate(poly) if p['refl'] ))

    poly.rotate(-first)

    print poly[0]

    poly_ = [ p['pt'] for p in poly ]

    num_ = len(poly_)

    pairs = set()

    for i, p in enumerate(poly):
        if p['refl']:
            vp = filter(lambda v: v['idx'] is not None, vis_poly_wrapper(poly_, i))
            # idx bezieht sich auf poly_

            for v in vp[1:]:
                pairs.add(tuple(sorted([vp[0]['idx'], v['idx']])))

    print pairs

    # init um die reflexe

    subs = defaultdict(SubP)

    for i, p in enumerate(poly):
        if p['refl']:
            for j in [-2, -1, 1, 2]:
                k = i+j

                if k > 0 and k < num_:
                    a, b = sorted([i, k])

                    print (a, b)

                    _ = SubP()
                    _.w = 0

                    if b-a == 1:
                        print 'edge'

                    else:
                        c = a+1
                        print 'wedge', c, (a, b) in pairs

                        _.S.append((c, c))

                    subs[(a, b)] = _

    def forw (i, j, k):
        if (i, j) not in pairs:
            return

        a = j

        p = subs[(i, j)]
        w = p.w

        if k-j > 1:
            if (j, k) not in pairs:
                return

            w += subs[(j, k)].w+1

        if j-i > 1:
            assert len(p.S) > 0

            if not is_refl(j, k, p.S[-1][1]):
                try:
                    while not is_refl(j, k, p.S[-2][1]):
                        print 'pop'
                        p.S_tail.append(p.S[-1])
                        del p.S[-1]
                except IndexError:
                    pass

                if p.S and not is_refl(i, p.S[-1][0], k):
                    a = p.S[-1][0]
                else:
                    w += 1
            else:
                w += 1

        subs[(i, k)].add_pair((a, j), w)

    def backw (i, j, k):
        if (j, k) not in pairs:
            return

        a = j

        p = subs[(j, k)]
        w = p.w

        if j-i > 1:
            if (i, j) not in pairs:
                return

            w += subs[(i, j)].w+1

        if k-j > 1:
            assert len(p.S) > 0

            if not is_refl(j, p.S[0][0], i):
                try:
                    while not is_refl(j, p.S[1][0], i):
                        print 'pop'
                        p.S_head.append(p.S[0])
                        del p.S[0]
                except IndexError:
                    pass

                if p.S and not is_refl(k, i, p.S[0][1]):
                    a = p.S[0][1]
                else:
                    w += 1
            else:
                w += 1

        subs[(i, k)].add_pair((j, a), w)

    ids = range(num_)

    for l in range(3, num_):
        for i in ids[:-l]:
            if poly[i]['refl']:
                k = i+l

                if (i, k) in pairs:
                    if poly[k]['refl']:
                        for j in range(i+1, k):
                            forw(i, j, k)

                    else:
                        for j in range(i+1, k-1):
                            if poly[j]['refl']:
                                forw(i, j, k)

                        forw(i, k-1, k)

        for k in ids[l:]:
            if poly[k]['refl']:
                i = k-l

                if (i, k) in pairs:
                    if not poly[i]['refl']:
                        backw(i, i+1, k)

                        for j in range(i+2, k):
                            if poly[j]['refl']:
                                backw(i, j, k)


    def recover (i, k):
        if k-i < 2:
            return

        s_a = subs[(i, k)]

        if poly[i]['refl']:
            j = s_a.S[-1][1]

            recover(j, k)

            if j-i > 1:
                if s_a.S[-1][0] != s_a.S[-1][1]:
                    s_b = subs[(i, j)]
                    s_b.restore_S()

                    while s_b.S and s_a.S[-1][0] != s_b.S[-1][0]:
                        del s_b.S[-1]

                recover(i, j)

        else:
            j = s_a.S[0][0]

            recover(i, j)

            if k-j > 1:
                if s_a.S[0][0] != s_a.S[0][1]:
                    s_b = subs[(j, k)]
                    s_b.restore_S()

                    while s_b.S and s_a.S[0][1] != s_b.S[0][1]:
                        del s_b.S[0]

                recover(j, k)

    recover(0, num-1)

    diags = []

    def collect (i, k):
        if k-i < 2:
            return

        s = subs[(i, k)].S

        if poly[i]['refl']:
            j = s[-1][1]
            a = j == s[-1][0]
            b = True
        else:
            j = s[0][0]
            b = j == s[0][1]
            a = True

        if a and j-i > 1:
            diags.append((i, j))

        if b and k-j > 1:
            diags.append((j, k))

        collect(i, j)
        collect(j, k)

    collect(0, num-1)

    diags.sort(key=lambda v: (v[0], -v[1]))

    i = 0
    p = 0
    q = 0

    ps = []
    rs = []

    decs = [[]]

    while i < num:
        new_p = decs[p]

        if not new_p or new_p[-1] != i:
            new_p.append(i)

        if rs and i == diags[rs[-1]][1]:
            if new_p[0] != diags[rs[-1]][0]:
                new_p.append(diags[rs[-1]][0])

            del rs[-1]

            p = ps[-1]
            del ps[-1]

        elif q < len(diags) and i == diags[q][0]:
            if new_p[0] != diags[q][1]:
                new_p.append(diags[q][1])

            decs.append([])

            rs.append(q)
            ps.append(p)

            q += 1
            p = q

        else:
            i += 1

    print diags

    res = []

    for dec in decs:
        p = [ original[poly[d]['idx']] for d in dec ]

        add_internals(original, p)

        res.append(p)

    return res

def simple_rm_internals (poly):
    mark_internals(poly, None)
    rm_internals(poly)

def scale (poly):
    num = len(poly)

    ls = []

    for i in range(num):
        j = (i+1)%num

        pA = poly[i]
        pB = poly[j]

        v = [pB['pt'][0]-pA['pt'][0], pB['pt'][1]-pA['pt'][1]]
        l = normalize(v)

        ls.append(l)

    f = 1/min(ls)

    print 'f', f

    if f > 1:
        for p in poly:
            p['pt'][0] *= f
            p['pt'][1] *= f

if __name__ == '__main__':

    env = Environment()
    env.filters['to_path'] = to_path

    tmpl = env.from_string('''<?xml version="1.0" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
    version="1.1"
    height="{{ size[1] }}"
    width="{{ size[0] }}">
    <path d="{{ poly|to_path(-pos[0], -pos[1]) }}z" style="stroke:black; fill:none;"/>
    <g>
        {% for d in decs %}
        <path d="{{ d|map(attribute="pt")|list|to_path(-pos[0], -pos[1]) }}z" style="stroke:black; fill:#faa; stroke-width:1;"/>
        {% endfor %}
    </g>
</svg>
        ''')

    with open('../../vp/dev/complex.json', 'r') as f:
        polys = json.load(f)['polys']

        for i, poly in enumerate(polys):
            if i != 2:
                pass#continue

            num = len(poly)

            for j in range(1, num):
                poly[j][0] += poly[j-1][0]
                poly[j][1] += poly[j-1][1]

            assert is_cw(poly)

            decs = decompose(poly)

            xs = [ pt[0] for pt in poly ]
            ys = [ pt[1] for pt in poly ]

            min_x = min(xs)
            min_y = min(ys)

            max_x = max(xs)
            max_y = max(ys)

            with open('res/dec_%i.svg' % i, 'w') as out:
                out.write(tmpl.render({ 'poly': poly, 'decs': decs, 'pos': (min_x, min_y), 'size': (max_x-min_x, max_y-min_y) }))
