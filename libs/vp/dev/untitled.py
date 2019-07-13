#!/usr/bin/env python3
# *-* coding: UTF-8 *-*

# Copyright 2012-2019 Ronald Römer
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
