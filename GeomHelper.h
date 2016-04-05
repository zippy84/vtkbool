/*
   Copyright 2012-2016 Ronald Römer

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#ifndef __GeomHelper_h
#define __GeomHelper_h

#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkKdTreePointLocator.h>
#include <vector>

#define CPY(a, b) a[0] = b[0]; a[1] = b[1]; a[2] = b[2];

namespace GeomHelper {

    void ComputeNormal (vtkPoints *pts, vtkIdList *poly, double* n);
    double GetAngle (double *vA, double *vB, double *n);

    void RemoveCells (vtkPolyData *pd, std::vector<int> &cells);

    void FindPoints (vtkKdTreePointLocator *pl, const double *pt, vtkIdList *pts, double tol = 1e-7);

    void WriteVTK (const char *name, vtkPolyData *pd);

};

#endif

