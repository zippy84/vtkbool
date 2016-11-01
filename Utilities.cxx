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

#include <cmath>

#include <vector>
#include <algorithm>

#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkDataWriter.h>

#include "Utilities.h"

namespace Utilities {

    void ComputeNormal (vtkPoints *pts, vtkIdList *poly, double* n) {
        n[0] = 0; n[1] = 0; n[2] = 0;
        double p0[3], p1[3];

        pts->GetPoint(poly->GetId(0), p0);

        int numPts = poly->GetNumberOfIds();

        for (int i = 0; i < numPts; i++) {
            pts->GetPoint(poly->GetId((i+1)%numPts), p1);

            n[0] += (p0[1]-p1[1])*(p0[2]+p1[2]);
            n[1] += (p0[2]-p1[2])*(p0[0]+p1[0]);
            n[2] += (p0[0]-p1[0])*(p0[1]+p1[1]);

            CPY(p0, p1)
        }

        vtkMath::Normalize(n);
    }

    double GetAngle (double *vA, double *vB, double *n) {
        double _vA[3];

        vtkMath::Cross(n, vA, _vA);
        double ang = std::atan2(vtkMath::Dot(_vA, vB), vtkMath::Dot(vA, vB));

        if (ang < 0) {
            ang += 2*pi;
        }
        return ang;
    }

    void FindPoints (vtkKdTreePointLocator *pl, const double *pt, vtkIdList *pts, double tol) {
        pts->Reset();

        vtkPolyData *pd = vtkPolyData::SafeDownCast(pl->GetDataSet());

        vtkIdList *closest = vtkIdList::New();

        // vtkKdTree.cxx#L2529
        // arbeitet mit single-precision
        pl->FindPointsWithinRadius(std::max(1e-5, tol), pt, closest);

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

        f << std::setprecision(12);

        double pt[3];
        for (int i = 0; i < numPts; i++) {
            pts->GetPoint(i, pt);
            f << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
        }

        w->WriteCells(&f, pd->GetLines(), "LINES");
        w->WriteCells(&f, pd->GetPolys(), "POLYGONS");
        w->WriteCells(&f, pd->GetStrips(), "TRIANGLE_STRIPS");
        w->WriteCellData(&f, pd);

        f.close();

        w->Delete();
    }

    // ermöglicht mod mit neg. int
    double Mod (int a, int b) {
        return (double) (((a%b)+b)%b);
    }

    double GetD (double *ptA, double *ptB) {
        double v[] = {ptA[0]-ptB[0], ptA[1]-ptB[1]};
        return sqrt(v[0]*v[0]+v[1]*v[1]);
    }

    double GetD3 (double *ptA, double *ptB) {
        double v[] = {ptA[0]-ptB[0], ptA[1]-ptB[1], ptA[2]-ptB[2]};
        return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    }

    double Normalize (double *p) {
        double n = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
        p[0] /= n;
        p[1] /= n;
        p[2] /= n;

        return n;
    }

    void GetNormal (double pts[][3], double *n, const int num) {
        n[0] = 0; n[1] = 0, n[2] = 0;
        double *p0 = pts[0],
            *p1;

        for (int i = 0; i < num; i++) {
            p1 = pts[(i+1)%num];

            n[0] += (p0[1]-p1[1])*(p0[2]+p1[2]);
            n[1] += (p0[2]-p1[2])*(p0[0]+p1[0]);
            n[2] += (p0[0]-p1[0])*(p0[1]+p1[1]);

            p0 = p1;
        }

        Normalize(n);
    }

}
