#!/usr/bin/env python

# Copyright 2012-2019 Ronald Römer

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

__all__ = ['parse_svg']

import re

pat = re.compile(r'([a-z]?)\s*([\d\-.,e]+)', flags=re.I)

def parse_svg(path):
    pts = []
    last_cmd = None

    for seg in pat.finditer(path):
        cmd = seg.group(1)
        data = seg.group(2)

        if cmd == '':
            cmd = last_cmd

        if len(pts) == 0:
            pts.append(list(map(lambda c: float(c), data.split(','))))
            last_cmd = 'l' if cmd == 'm' else 'L'
            continue
        elif cmd == 'h':
            x = float(data)
            pts.append([pts[-1][0]+x, pts[-1][1]])
        elif cmd == 'v':
            y = float(data)
            pts.append([pts[-1][0], pts[-1][1]+y])
        elif cmd == 'l':
            x, y = list(map(lambda c: float(c), data.split(',')))
            pts.append([pts[-1][0]+x, pts[-1][1]+y])
        elif cmd == 'H':
            x = float(data)
            pts.append([x, pts[-1][1]])
        elif cmd == 'V':
            y = float(data)
            pts.append([pts[-1][0], y])
        elif cmd == 'L':
            pts.append(list(map(lambda c: float(c), data.split(','))))

        last_cmd = cmd

    return pts

if __name__ == '__main__':
    path = 'm 5.96086550924,5.37468641552 -0.2590305926,1.08097294886 -1.38149649392,0.0186374646 0.25903059261,-1.584184494 2.45215627665,-0.44729915123 0.0863435313,-1.11824787813 -1.81321414823,-0.0186374646 -1.70960191125,1.34189745373 -0.31083671118,1.56554702937 1.03612237014,0.26092450488 1.38149649392,-0.0186374646 1.31242166954,-0.42866168674 0.0345374125,1.71464674644 -2.41761886442,0.50321154516 -1.7786767359,-0.76413605004 -0.60440471608,-2.21785829161 1.29515296308,-2.23649575624 2.6939181631,-0.98778562566 2.83206781257,1.47235970619 -1.15700331373,1.88238392817'

    #path = 'm 10.9425188625,2.62092799516 h 3.2695218222 l 2.8358097439,1.15660714062 -1.4012236382,2.44928570895 v 2.34723213787 l -2.5526047101,-0.94732330591 h -2.5518528288 l 1.8015732491,-2.11428383042 -2.1018354572,-0.7824107127 z'

    poly = parse_svg(path)

    #print(poly)

    print('M{}Z'.format(' '.join(map(lambda pt: ','.join(map(str, pt)), poly))))

    # wenn ccw
    #poly.reverse()

    # wenn cw
    poly.append(poly[0])
    del poly[0]

    poly_ = iter(poly)
    next(poly_)

    for a, b in zip(poly, poly_):
        a[0] -= b[0]
        a[1] -= b[1]

    poly.reverse()

    print('m{}z'.format(' '.join(map(lambda pt: ','.join(map(str, pt)), poly))))

    print(poly)
