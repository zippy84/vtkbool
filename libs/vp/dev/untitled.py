#!/usr/bin/env python3
import sys
import re
import json
from collections import defaultdict

from svg import *

data = []

for line in sys.stdin:
    ln = line.rstrip()
    # print('>', '"%s"' % ln)

    if ln.startswith('?'):
        if ln == '?':
            data.append(defaultdict(list))
        elif ln.startswith('?X'):
            data[-1]['X'] = int(ln[3:])
        elif ln.startswith('?_'):
            data[-1]['_'].append(ln[3:])
        else:
            m = re.match(r'([A-Z]\s*\d*)\s+(.+)', ln[1:])

            if m:
                poly = parse_svg(m.group(2))

                if len(poly) > 0:
                    data[-1][m.group(1).replace(' ', '')] = poly

print('var data = {};'.format(json.dumps(data)))
