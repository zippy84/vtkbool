/*
Copyright 2012-2020 Ronald Römer

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

#ifndef __Utilities_h
#define __Utilities_h

#include <iostream>

#include <vtkPolyData.h>
#include <vtkKdTreePointLocator.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkMath.h>

#include "Tools.h"

double GetAngle (double *vA, double *vB, double *n);
double GetD (double *a, double *b);

/* VTK */
void ComputeNormal (vtkPoints *pts, double *n, vtkIdList *poly = nullptr);
void FindPoints (vtkKdTreePointLocator *pl, const double *pt, vtkIdList *pts, double tol = 1e-6);
void WriteVTK (const char *name, vtkPolyData *pd);

inline void ComputeNormal2 (vtkPolyData *pd, double *n, vtkIdType num, const vtkIdType *poly) {
    n[0] = 0; n[1] = 0; n[2] = 0;

    if (num == 3) {
        double p0[3], p1[3], p2[3], a[3], b[3];

        pd->GetPoint(poly[0], p0);
        pd->GetPoint(poly[1], p1);
        pd->GetPoint(poly[2], p2);

        vtkMath::Subtract(p1, p0, a);
        vtkMath::Subtract(p2, p0, b);

        vtkMath::Cross(a, b, n);
    } else {
        double p0[3], p1[3];

        for (int i = 0; i < num; i++) {
            vtkIdType a = poly[i],
                b = poly[(i+1)%num];

            pd->GetPoint(a, p0);
            pd->GetPoint(b, p1);

            n[0] += (p0[1]-p1[1])*(p0[2]+p1[2]);
            n[1] += (p0[2]-p1[2])*(p0[0]+p1[0]);
            n[2] += (p0[0]-p1[0])*(p0[1]+p1[1]);
        }
    }

    vtkMath::Normalize(n);
}

/* Misc */
double Mod (int a, int b);

class Base {
public:
    Base (vtkPoints *pts, vtkIdList *poly);
    Base () {}
    double n[3], ei[3], ej[3], d;
};

void Transform (const double *in, double *out, Base &base);
// void BackTransform (const double *in, double *out, Base &base);

#endif
