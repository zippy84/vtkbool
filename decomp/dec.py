#!/usr/bin/env python
# *-* coding: UTF-8 *-*

import sys
import json

from collections import deque, defaultdict

sys.path.append('../vp_new_try')

from vis_poly import vis_poly_wrapper
from rm_trivials import rm_internals
from tools import cross, to_abs_path, is_cw, from_path

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
    poly = deque([ { 'pt': pt, 'idx': i } for i, pt in enumerate(pts) ])

    rm_internals(poly)

    num = len(poly)

    def is_refl (a, b, c):
        return cross(poly[a]['pt'], poly[b]['pt'], poly[c]['pt']) < 0

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
                        subs[(a, b)] = _

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


    print subs[(0, num-1)].w

with open('../vp_new_try/complex.json', 'r') as f:
    polys = json.load(f)['polys']

    for i, poly in enumerate(polys):
        if i != 1:
            continue

        num = len(poly)

        for j in range(1, num):
            poly[j][0] += poly[j-1][0]
            poly[j][1] += poly[j-1][1]

        assert(is_cw(poly))

        decompose(poly)

        #print to_abs_path(poly)
