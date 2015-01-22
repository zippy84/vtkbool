/*
   Copyright 2012-2014 Ronald Römer

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

#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>
#include <algorithm>

#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkTriangleStrip.h>
#include <vtkCellArray.h>
#include <vtkIntArray.h>
#include <vtkDataWriter.h>

#include "GeomHelper.h"

void GeomHelper::ComputeNormal (vtkPoints *pts, vtkIdList *poly, double* n) {

    n[0] = 0; n[1] = 0; n[2] = 0;
    double pt0[3], pt1[3];

    pts->GetPoint(poly->GetId(0), pt0);

    unsigned int nbr = poly->GetNumberOfIds();

    for (unsigned int i = 0; i < nbr; i++) {
        pts->GetPoint(poly->GetId((i+1)%nbr), pt1);

        n[0] += (pt0[1]-pt1[1])*(pt0[2]+pt1[2]);
        n[1] += (pt0[2]-pt1[2])*(pt0[0]+pt1[0]);
        n[2] += (pt0[0]-pt1[0])*(pt0[1]+pt1[1]);

        pt0[0] = pt1[0];
        pt0[1] = pt1[1];
        pt0[2] = pt1[2];
    }

    vtkMath::Normalize(n);
}

double GeomHelper::GetAngle (double *vA, double *vB, double *n) {
    // vA drehen
    double _vA[3];

    vtkMath::Cross(n, vA, _vA);

    double ang = std::atan2(vtkMath::Dot(_vA, vB), vtkMath::Dot(vA, vB));

    if (ang < 0) {
        ang += 2*M_PI;
    }

    // der winkel ist pos. und bewegt sich zw. 0 und 2pi
    return ang;
}

void GeomHelper::RemoveCells (vtkPolyData *pd, std::vector<int> &cells) {

    vtkIntArray *types = vtkIntArray::New();

    for (unsigned int i = 0; i < pd->GetNumberOfCells(); i++) {
        types->InsertNextValue(pd->GetCellType(i));
    }

    vtkCellArray *verts = vtkCellArray::New();
    verts->DeepCopy(pd->GetVerts());
    pd->GetVerts()->Initialize();

    vtkCellArray *lines = vtkCellArray::New();
    lines->DeepCopy(pd->GetLines());
    pd->GetLines()->Initialize();

    vtkCellArray *polys = vtkCellArray::New();
    polys->DeepCopy(pd->GetPolys());
    pd->GetPolys()->Initialize();

    vtkCellArray *strips = vtkCellArray::New();
    strips->DeepCopy(pd->GetStrips());
    pd->GetStrips()->Initialize();

    pd->DeleteCells();

    vtkCellData *cellData = vtkCellData::New();
    cellData->CopyAllocate(pd->GetCellData());

    int type;
    vtkIdType *pts, n;

    verts->InitTraversal();
    lines->InitTraversal();
    polys->InitTraversal();
    strips->InitTraversal();

    int cellId;

    for (unsigned int i = 0; i < types->GetNumberOfTuples(); i++) {
        type = types->GetValue(i);

        bool keep = std::count(cells.begin(), cells.end(), i) == 0;

        if (type == VTK_VERTEX || type == VTK_POLY_VERTEX) {
            verts->GetNextCell(n, pts);

            if (keep) {
                cellId = pd->InsertNextCell(type, n, pts);
                cellData->CopyData(pd->GetCellData(), i, cellId);
            }

        } else if (type == VTK_LINE || type == VTK_POLY_LINE) {
            lines->GetNextCell(n, pts);

            if (keep) {
                cellId = pd->InsertNextCell(type, n, pts);
                cellData->CopyData(pd->GetCellData(), i, cellId);
            }

        } else if (type == VTK_POLYGON || type == VTK_TRIANGLE || type == VTK_QUAD) {
            polys->GetNextCell(n, pts);

            if (keep) {
                cellId = pd->InsertNextCell(type, n, pts);
                cellData->CopyData(pd->GetCellData(), i, cellId);
            }

        } else if (type == VTK_TRIANGLE_STRIP) {
            strips->GetNextCell(n, pts);

            if (keep) {
                cellId = pd->InsertNextCell(type, n, pts);
                cellData->CopyData(pd->GetCellData(), i, cellId);
            }

        }

    }

    cellData->Squeeze();

    pd->GetCellData()->ShallowCopy(cellData);

    cellData->Delete();

    strips->Delete();
    polys->Delete();
    lines->Delete();
    verts->Delete();
    types->Delete();

}

void GeomHelper::FindPoints (vtkKdTreePointLocator *pl, const double *pt, vtkIdList *pts, double tol) {
    pts->Reset();

    vtkPolyData *pd = vtkPolyData::SafeDownCast(pl->GetDataSet());

    vtkIdList *closest = vtkIdList::New();

    // vtkKdTree.cxx#L2529
    // arbeitet mit single-precision
    pl->FindPointsWithinRadius(std::max(1e-5, tol), pt, closest);

    double c[3], v[3];

    for (unsigned int i = 0; i < closest->GetNumberOfIds(); i++) {
        pd->GetPoint(closest->GetId(i), c);
        vtkMath::Subtract(pt, c, v);

        if (vtkMath::Norm(v) < tol) {
            pts->InsertNextId(closest->GetId(i));
        }
    }

    closest->Delete();
}

void GeomHelper::WriteVTK (const char *name, vtkPolyData *pd) {
    vtkDataWriter *w = vtkDataWriter::New();

    vtkPoints *pts = pd->GetPoints();

    std::ofstream f(name);
    f << "# vtk DataFile Version 3.0\n"
      << "vtk output\n"
      << "ASCII\n"
      << "DATASET POLYDATA\n"
      << "POINTS " << pts->GetNumberOfPoints() << " double\n";

    f << std::setprecision(12);

    double pt[3];
    for (unsigned int i = 0; i < pts->GetNumberOfPoints(); i++) {
        pts->GetPoint(i, pt);
        f << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    }

    w->WriteCells(&f, pd->GetLines(), "LINES");
    w->WriteCells(&f, pd->GetPolys(), "POLYGONS");
    w->WriteCellData(&f, pd);

    f.close();

    w->Delete();
}
