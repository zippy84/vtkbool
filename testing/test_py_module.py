#!/usr/bin/env python

import sys

if sys.platform.startswith('linux'):
    sys.path.extend(sys.argv[1:])

elif sys.platform.startswith('win'):
    sys.path.extend([sys.argv[1],
        '/'.join(sys.argv[2:])])

import vtkboolPython
vtkboolPython.vtkPolyDataBooleanFilter()

#import vtk
