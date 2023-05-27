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

#include <vtkCellData.h>
#include <vtkCellArrayIterator.h>
#include <vtkIdTypeArray.h>

#include "vtkPolyDataBooleanFilter.h"

Poly Draw (double r, vtkIdType step, double x, double y, double rotate = 0) {
    Poly poly;

    double phi = 2*M_PI/static_cast<double>(step);

    vtkIdType i;

    double x1, y1, x2, y2;

    double alpha = rotate*M_PI/180;

    for (i = 0; i < step; i++) {
        double _i = static_cast<double>(i);

        x1 = r*std::cos(_i*phi);
        y1 = r*std::sin(_i*phi);

        x2 = std::cos(alpha)*x1-std::sin(alpha)*y1;
        y2 = std::sin(alpha)*x1+std::cos(alpha)*y1;

        poly.emplace_back(x2+x, y2+y, 0);
    }

    return poly;
}

bool Test (PolysType &polys, vtkIdType numCells) {
    auto pts = vtkSmartPointer<vtkPoints>::New();

    auto pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(pts);
    pd->Allocate(1);

    auto cell = vtkSmartPointer<vtkIdList>::New();

    for (const auto &p : polys[0]) {
        cell->InsertNextId(pts->InsertNextPoint(p.x, p.y, p.z));
    }

    pd->InsertNextCell(VTK_POLYGON, cell);

    auto cellIds = vtkSmartPointer<vtkIdTypeArray>::New();
    cellIds->SetName("OrigCellIds");
    cellIds->InsertNextValue(0);
    pd->GetCellData()->SetScalars(cellIds);

    PStrips pStrips {pd, 0};

    Base &base = pStrips.base;

    base.ei[0] = 1; base.ei[1] = 0; base.ei[2] = 0;
    base.ej[0] = 0; base.ej[1] = 1; base.ej[2] = 0;
    base.n[0] = 0; base.n[1] = 0; base.n[2] = 1;
    base.d = 0;

    StripsType holes;

    vtkIdType ind {0};
    std::size_t stripId {0};

    PolysType::iterator itr;
    for (itr = polys.begin()+1; itr != polys.end(); ++itr) {
        StripType strip;

        Poly &poly = *itr;

        poly.push_back(poly.front());

        for (const auto &p : poly) {
            StripPt sp;
            sp.ind = ind;
            sp.polyId = 0;

            sp.pt[0] = p.x;
            sp.pt[1] = p.y;
            sp.pt[2] = p.z;

            pStrips.pts.emplace(ind, std::move(sp));

            strip.emplace_back(ind, stripId);

            ind++;
        }

        holes.push_back(strip);

        stripId++;
    }

    IdsType descIds {0};

    try {
        Merger(pd, pStrips, holes, descIds, 0).Run();
    } catch (const std::runtime_error &e) {
        return false;
    }

    // dafür ist der Merger nicht zuständig
    pd->RemoveDeletedCells();

    {
        // schneiden sich die einzelnen polygone selbst?

        vtkIdType i, num;
        const vtkIdType *cell;

        double pt[3];

        auto cellItr = vtk::TakeSmartPointer(pd->GetPolys()->NewIterator());

        for (cellItr->GoToFirstCell(); !cellItr->IsDoneWithTraversal(); cellItr->GoToNextCell()) {
            cellItr->GetCurrentCell(num, cell);

            std::set<Point3d> poly;

            for (i = 0; i < num; i++) {
                pd->GetPoint(cell[i], pt);
                if (!poly.emplace(pt[0], pt[1], pt[2]).second) {
                    // schlägt fehl, wenn bereits vorhanden
                    return false;
                }
            }

        }
    }

    // WriteVTK("merged.vtk", pd);

    if (pd->GetNumberOfCells() != numCells) {
        return false;
    }

    return true;
}

int main() {
    PolysType polysA {
        Draw(8, 18, 0, 0),
        Draw(.25, 6, 3.5, 0),
        Draw(.25, 6, -3.5, 0),

        Draw(.5, 6, 1, 0),
        Draw(.5, 6, -1, 0),
        Draw(.5, 6, 0, 1, 30),
        Draw(.5, 6, 0, -1, 30)
    };

    WritePolys("polysA.vtk", polysA);

    Poly a = Draw(4, 18, 0, 0, 10),
        b = Draw(3, 18, 0, 0, 10);

    std::copy(b.rbegin(), b.rend(), std::back_inserter(a));

    PolysType polysB {
        Poly{{5, 5, 0}, {-5, 5, 0}, {-5, -5, 0}, {5, -5, 0}},
        Draw(2, 18, 0, 0, 5),
        a
    };

    WritePolys("polysB.vtk", polysB);

    if (!Test(polysA, 10)) {
        return EXIT_FAILURE;
    }

    if (!Test(polysB, 5)) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
