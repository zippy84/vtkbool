#!/usr/bin/env python
# *-* coding: UTF-8 *-*

# Copyright 2012-2022 Ronald Römer
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
sys.path.extend(['/home/zippy/vtkbool/build/lib/python3.10/site-packages/vtkbool'])

from vtkmodules.vtkIOLegacy import vtkPolyDataReader, vtkPolyDataWriter
from vtkBool import vtkPolyDataBooleanFilter

readerA = vtkPolyDataReader()
readerA.SetFileName('_test/AnteriorToInferiorModelSequence.vtk')

readerB = vtkPolyDataReader()
readerB.SetFileName('_test/InferiorAir.vtk')

bf = vtkPolyDataBooleanFilter()
bf.SetInputConnection(0, readerA.GetOutputPort())
bf.SetInputConnection(1, readerB.GetOutputPort())

bf.SetOperModeToNone()

writer_ = vtkPolyDataWriter()
writer_.SetFileName('_test/lines.vtk')
writer_.SetInputConnection(bf.GetOutputPort(2))
writer_.Update()

# das ursprüngliche problem war bei cell 151227

writer0 = vtkPolyDataWriter()
writer0.SetFileName('_test/result0.vtk')
writer0.SetInputConnection(bf.GetOutputPort())
writer0.Update()

writer1 = vtkPolyDataWriter()
writer1.SetFileName('_test/result1.vtk')
writer1.SetInputConnection(bf.GetOutputPort(2))
writer1.Update()

# fehler bei 544,545,546,547,548,690,691,692,693
# sind lines in holes
