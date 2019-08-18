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

import sys

if __name__ == '__main__':
    assert(os.environ.get('APPVEYOR') == 'True')

    assert(len(sys.argv) == 3)

    arg_zip = sys.argv[1]
    arg_dir = sys.argv[2]

    filename = arg_zip.rsplit('/', 1)[-1]
    name = filename[:-4]

    if not os.path.exists('{}.tar'.format(arg_dir)):
        print('Downloading VTK')

        f = urllib.request.urlopen(arg_zip)

        with open(filename, 'wb') as zf:
            zf.write(f.read())

        with zipfile.ZipFile(filename) as zf:
            zf.extractall()

        os.mkdir('{}/build'.format(name))

        print('Configuring VTK')

        proc = subprocess.Popen(shlex.split('cmake -Wno-deprecated -Wno-dev -DVTK_PYTHON_VERSION=3 -DVTK_WRAP_PYTHON=ON -DCMAKE_BUILD_TYPE=Release ..'), cwd='{}/build'.format(name))
        proc.wait()
    else:
        print('{}.tar already exists'.format(arg_dir))

        with tarfile.open('{}.tar'.format(arg_dir), 'r:') as zf:
            zf.extractall(name)

        os.remove('{}.tar'.format(arg_dir))

        print('Building VTK')

        proc = subprocess.Popen(shlex.split('cmake --build . --config Release'), creationflags=subprocess.CREATE_NEW_PROCESS_GROUP, cwd='{}/build'.format(name))

        try:
            proc.wait(3600)
        except subprocess.TimeoutExpired:
            print('cmake will be killed')

            proc.send_signal(signal.CTRL_BREAK_EVENT)
            proc.kill()
            proc.wait()

        else:
            print('Ready :)')

            proc2 = subprocess.Popen(shlex.split('cmake --install . --prefix C:/projects/vtkbool/{}'.format(arg_dir)), cwd='{}/build'.format(name))
            proc2.wait()

    print('Creating {}.tar'.format(arg_dir))

    shutil.make_archive(arg_dir, 'tar', name)
