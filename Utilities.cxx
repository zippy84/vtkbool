/*
Copyright 2012-2019 Ronald Römer

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

#include <cmath>

#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkDataWriter.h>

void ComputeNormal (vtkPoints *pts, double *n, vtkIdList *poly) {
    n[0] = 0; n[1] = 0; n[2] = 0;
    double p0[3], p1[3];

    int numPts = pts->GetNumberOfPoints();

    vtkIdList *_poly = poly;

    if (poly == nullptr) {
        _poly = vtkIdList::New();
        _poly->SetNumberOfIds(numPts);
        for (int i = 0; i < numPts; i++) {
            _poly->SetId(i, i);
        }
    } else {
        numPts = poly->GetNumberOfIds();
    }

    pts->GetPoint(_poly->GetId(0), p0);

    for (int i = 0; i < numPts; i++) {
        pts->GetPoint(_poly->GetId((i+1)%numPts), p1);

        n[0] += (p0[1]-p1[1])*(p0[2]+p1[2]);
        n[1] += (p0[2]-p1[2])*(p0[0]+p1[0]);
        n[2] += (p0[0]-p1[0])*(p0[1]+p1[1]);

        Cpy(p0, p1, 3);
    }

    vtkMath::Normalize(n);

    if (poly == nullptr) {
        _poly->Delete();
    }
}

void FindPoints (vtkKdTreePointLocator *pl, const double *pt, vtkIdList *pts, double tol) {
    pts->Reset();

    vtkPolyData *pd = vtkPolyData::SafeDownCast(pl->GetDataSet());

    vtkIdList *closest = vtkIdList::New();

    // vtkKdTree.cxx#L2505
    // arbeitet mit single-precision
    pl->FindPointsWithinRadius(std::max(1e-3, tol), pt, closest);

    int numPts = closest->GetNumberOfIds();

    double c[3], v[3];

    for (int i = 0; i < numPts; i++) {
        pd->GetPoint(closest->GetId(i), c);
        vtkMath::Subtract(pt, c, v);

        if (vtkMath::Norm(v) < tol) {
            pts->InsertNextId(closest->GetId(i));
        }
    }

    closest->Delete();
}

void WriteVTK (const char *name, vtkPolyData *pd) {
    vtkDataWriter *w = vtkDataWriter::New();

    vtkPoints *pts = pd->GetPoints();

    int numPts = pts->GetNumberOfPoints();

    std::ofstream f(name);
    f << "# vtk DataFile Version 3.0\n"
      << "vtk output\n"
      << "ASCII\n"
      << "DATASET POLYDATA\n"
      << "POINTS " << numPts << " double\n";

    f << std::setprecision(8);

    double pt[3];
    for (int i = 0; i < numPts; i++) {
        pts->GetPoint(i, pt);
        f << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    }

    w->WriteCells(&f, pd->GetLines(), "LINES");
    w->WriteCells(&f, pd->GetPolys(), "POLYGONS");
    w->WriteCells(&f, pd->GetStrips(), "TRIANGLE_STRIPS");
    w->WriteCellData(&f, pd);
    w->WritePointData(&f, pd);

    f.close();

    w->Delete();
}

double GetAngle (double *vA, double *vB, double *n) {
    // http://math.stackexchange.com/questions/878785/how-to-find-an-angle-in-range0-360-between-2-vectors

    double _vA[3];

    vtkMath::Cross(n, vA, _vA);
    double ang = std::atan2(vtkMath::Dot(_vA, vB), vtkMath::Dot(vA, vB));

    if (ang < 0) {
        ang += 2*PI;
    }

    return ang;
}

double GetD (double *ptA, double *ptB) {
    double v[] = {ptA[0]-ptB[0], ptA[1]-ptB[1], ptA[2]-ptB[2]};
    return std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

// ermöglicht mod mit neg. int
double Mod (int a, int b) {
    return (double) (((a%b)+b)%b);
}

Base::Base (vtkPoints *pts, vtkIdList *poly) {
    ComputeNormal(pts, n, poly);

    double ptA[3],
        ptB[3];

    pts->GetPoint(poly->GetId(0), ptA);
    pts->GetPoint(poly->GetId(1), ptB);

    ei[0] = ptB[0]-ptA[0];
    ei[1] = ptB[1]-ptA[1];
    ei[2] = ptB[2]-ptA[2];

    Normalize(ei, 3);

    ej[0] = n[1]*ei[2]-n[2]*ei[1];
    ej[1] = -n[0]*ei[2]+n[2]*ei[0];
    ej[2] = n[0]*ei[1]-n[1]*ei[0];

    Normalize(ej, 3);

    d = n[0]*ptA[0]+n[1]*ptA[1]+n[2]*ptA[2];
}

void Transform (const double *in, double *out, Base &base) {
    double x = base.ei[0]*in[0]+base.ei[1]*in[1]+base.ei[2]*in[2],
        y = base.ej[0]*in[0]+base.ej[1]*in[1]+base.ej[2]*in[2];

    out[0] = x;
    out[1] = y;
}

/*void BackTransform (const double *in, double *out, Base &base) {
    double x = in[0]*base.ei[0]+in[1]*base.ej[0]+base.d*base.n[0],
        y = in[0]*base.ei[1]+in[1]*base.ej[1]+base.d*base.n[1],
        z = in[0]*base.ei[2]+in[1]*base.ej[2]+base.d*base.n[2];

    out[0] = x;
    out[1] = y;
    out[2] = z;
}
*/
