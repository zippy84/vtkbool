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

import subprocess
import os
import signal

import urllib.request
import tarfile
import zipfile
import shlex
import shutil

from collections import deque

import sys

if not os.path.exists('VTK7.tar.gz'):
  print('Downloading VTK')

  f = urllib.request.urlopen('https://www.vtk.org/files/release/7.1/VTK-7.1.1.zip')

  with open('VTK-7.1.1.zip', 'wb') as zf:
    zf.write(f.read())

  with zipfile.ZipFile('VTK-7.1.1.zip') as zf:
    zf.extractall()

  os.mkdir('VTK-7.1.1/build')

  print('Configuring VTK')

  proc = subprocess.Popen(shlex.split('cmake -Wno-deprecated -Wno-dev -DVTK_PYTHON_VERSION=3 -DVTK_WRAP_PYTHON=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=C:/projects/vtkbool-1/VTK7 ..'), cwd='VTK-7.1.1/build')
  proc.wait()
else:
  print('VTK7.tar.gz already exists')

  with tarfile.open('VTK7.tar.gz', 'r:gz') as zf:
    zf.extractall('VTK-7.1.1')

  os.remove('VTK7.tar.gz')

with open('build.log', 'ab') as log:
  print('Building VTK')

  proc = subprocess.Popen(shlex.split('cmake --build . --config Release'), creationflags=subprocess.CREATE_NEW_PROCESS_GROUP, cwd='VTK-7.1.1/build')

  try:
    proc.wait(3600)
  except subprocess.TimeoutExpired:
    print('cmake will be killed')

    proc.send_signal(signal.CTRL_BREAK_EVENT)
    proc.kill()
    proc.wait()

  else:
    print('Ready :)')

    proc2 = subprocess.Popen(shlex.split('cmake --build . --config Release --target install'), cwd='VTK-7.1.1/build')
    proc2.wait()

print('Creating VTK7.tar.gz')

shutil.make_archive('VTK7', 'gztar', 'VTK-7.1.1')

try:
  with open('build.log') as log:
    for l in deque(log, 50):
      print(l.rstrip())
except FileNotFoundError:
  pass
