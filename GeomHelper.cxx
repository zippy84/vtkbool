/*
   Copyright 2012, 2013 Ronald Römer

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

// C++0x/C++11
#include <random>

#define _USE_MATH_DEFINES
#include <cmath>

#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkTriangleStrip.h>
#include <vtkCellArray.h>

#include "GeomHelper.h"

#define CPY(a, b) a[0] = b[0]; a[1] = b[1]; a[2] = b[2];

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

void GeomHelper::PreparePolyData (vtkPolyData *pd) {
    // nicht erwünscht sind Lines, Verts, Triangle-Strips, ... sowie CellData und PointData

    if (pd->GetCellData() != NULL) { pd->GetCellData()->Initialize(); }
    if (pd->GetPointData() != NULL) { pd->GetPointData()->Initialize(); }

    vtkCellArray *strips = pd->GetStrips();

    vtkIdType n;
    vtkIdType *pts;

    for (strips->InitTraversal(); strips->GetNextCell(n, pts);) {
        vtkTriangleStrip::DecomposeStrip(n, pts, pd->GetPolys());
    }

    int type;

    for (unsigned int i = 0; i < pd->GetNumberOfCells(); i++) {
        type = pd->GetCellType(i);

        if (type != VTK_POLYGON && type != VTK_QUAD && type != VTK_TRIANGLE) {
            pd->DeleteCell(i);
        }
    }

    pd->RemoveDeletedCells();

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

