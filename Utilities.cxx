/*
Copyright 2012-2022 Ronald Römer

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

#include "Utilities.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkDataWriter.h>

#include <vtkPolyDataWriter.h>

void ComputeNormal (vtkPoints *pts, double *n, vtkIdType num, const vtkIdType *poly) {
    n[0] = 0; n[1] = 0; n[2] = 0;

    if (num == 3) {
        double pA[3], pB[3], pC[3], a[3], b[3];

        pts->GetPoint(poly[0], pA);
        pts->GetPoint(poly[1], pB);
        pts->GetPoint(poly[2], pC);

        vtkMath::Subtract(pB, pA, a);
        vtkMath::Subtract(pC, pA, b);

        vtkMath::Cross(a, b, n);
    } else {
        double pA[3], pB[3];

        vtkIdType i, a, b;

        for (i = 0; i < num; i++) {
            a = poly[i];
            b = poly[i+1 == num ? 0 : i+1];

            pts->GetPoint(a, pA);
            pts->GetPoint(b, pB);

            n[0] += (pA[1]-pB[1])*(pA[2]+pB[2]);
            n[1] += (pA[2]-pB[2])*(pA[0]+pB[0]);
            n[2] += (pA[0]-pB[0])*(pA[1]+pB[1]);
        }
    }

    vtkMath::Normalize(n);
}

void FindPoints (vtkKdTreePointLocator *pl, const double *pt, vtkIdList *pts, double tol) {
    pts->Reset();

    vtkPolyData *pd = vtkPolyData::SafeDownCast(pl->GetDataSet());

    vtkIdList *closest = vtkIdList::New();

    // vtkKdTree.cxx#L2505
    // arbeitet mit single-precision
    pl->FindPointsWithinRadius(std::max(1e-3, tol), pt, closest);

    vtkIdType i, numPts = closest->GetNumberOfIds();

    double c[3], v[3];

    for (i = 0; i < numPts; i++) {
        pd->GetPoint(closest->GetId(i), c);
        vtkMath::Subtract(pt, c, v);

        if (vtkMath::Norm(v) < tol) {
            pts->InsertNextId(closest->GetId(i));
        }
    }

    closest->Delete();
}

void WriteVTK (const char *name, vtkPolyData *pd) {
    vtkPolyDataWriter *w = vtkPolyDataWriter::New();
    w->SetInputData(pd);
    w->SetFileName(name);
    w->Update();
    w->Delete();
}

double GetAngle (double *vA, double *vB, double *n) {
    // http://math.stackexchange.com/questions/878785/how-to-find-an-angle-in-range0-360-between-2-vectors

    double _vA[3];

    vtkMath::Cross(n, vA, _vA);
    double ang = std::atan2(vtkMath::Dot(_vA, vB), vtkMath::Dot(vA, vB));

    if (ang < 0) {
        ang += 2*M_PI;
    }

    return ang;
}

Base::Base (vtkPoints *pts, vtkIdType num, const vtkIdType *poly) {
    ComputeNormal(pts, n, num, poly);

    double ptA[3],
        ptB[3];

    pts->GetPoint(poly[0], ptA);
    pts->GetPoint(poly[1], ptB);

    ei[0] = ptB[0]-ptA[0];
    ei[1] = ptB[1]-ptA[1];
    ei[2] = ptB[2]-ptA[2];

    vtkMath::Normalize(ei);

    ej[0] = n[1]*ei[2]-n[2]*ei[1];
    ej[1] = -n[0]*ei[2]+n[2]*ei[0];
    ej[2] = n[0]*ei[1]-n[1]*ei[0];

    vtkMath::Normalize(ej);

    d = n[0]*ptA[0]+n[1]*ptA[1]+n[2]*ptA[2];
}

void Transform (const double *in, double *out, const Base &base) {
    double x = base.ei[0]*in[0]+base.ei[1]*in[1]+base.ei[2]*in[2],
        y = base.ej[0]*in[0]+base.ej[1]*in[1]+base.ej[2]*in[2];

    out[0] = x;
    out[1] = y;
}

void BackTransform (const double *in, double *out, const Base &base) {
    double x = in[0]*base.ei[0]+in[1]*base.ej[0]+base.d*base.n[0],
        y = in[0]*base.ei[1]+in[1]*base.ej[1]+base.d*base.n[1],
        z = in[0]*base.ei[2]+in[1]*base.ej[2]+base.d*base.n[2];

    out[0] = x;
    out[1] = y;
    out[2] = z;
}
