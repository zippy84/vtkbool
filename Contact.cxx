/*
Copyright 2012-2025 Ronald Römer

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

#include "Contact.h"

#include <vtkCleanPolyData.h>
#include <vtkCellIterator.h>
#include <vtkCellData.h>
#include <vtkTriangleStrip.h>

// #define _debA 0
// #define _debB 0

#if (defined(_debA) && defined(_debB))
vtkIdType _idA, _idB;
#endif

vtkSmartPointer<vtkPolyData> Clean (vtkPolyData *pd) {

    auto clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetOutputPointsPrecision(vtkAlgorithm::DOUBLE_PRECISION);
    clean->SetTolerance(1e-6);

    clean->ConvertLinesToPointsOff();
    clean->ConvertPolysToLinesOff();
    clean->ConvertStripsToPolysOff();

    clean->SetInputData(pd);

    clean->Update();

    auto cleaned = clean->GetOutput();

    vtkIdType numCells = cleaned->GetNumberOfCells();

    auto newPd = vtkSmartPointer<vtkPolyData>::New();
    newPd->SetPoints(cleaned->GetPoints());
    newPd->Allocate(numCells);

    auto cellIds = vtkSmartPointer<vtkIdTypeArray>::New();
    cellIds->SetName("OrigCellIds");
    cellIds->Allocate(numCells);

    vtkCellIterator *cellItr = cleaned->NewCellIterator();

    vtkIdType cellId;
    vtkIdList *ptIds;

    vtkIdType num;
    const vtkIdType *pts;

    for (cellItr->InitTraversal(); !cellItr->IsDoneWithTraversal(); cellItr->GoToNextCell()) {
        cellId = cellItr->GetCellId();
        ptIds = cellItr->GetPointIds();

        if (cellItr->GetCellType() == VTK_TRIANGLE || cellItr->GetCellType() == VTK_POLYGON) {
            newPd->InsertNextCell(cellItr->GetCellType(), ptIds);
            cellIds->InsertNextValue(cellId);

        } else if (cellItr->GetCellType() == VTK_TRIANGLE_STRIP) {
            auto cells = vtkSmartPointer<vtkCellArray>::New();

            vtkTriangleStrip::DecomposeStrip(cellItr->GetNumberOfPoints(), ptIds->GetPointer(0), cells);

            for (cells->InitTraversal(); cells->GetNextCell(num, pts);) {
                if (pts[0] != pts[1] && pts[1] != pts[2] && pts[2] != pts[0]) {
                    newPd->InsertNextCell(VTK_TRIANGLE, num, pts);
                    cellIds->InsertNextValue(cellId);
                }
            }

        } else if (cellItr->GetCellType() == VTK_QUAD) {

            double pA[3], pB[3], pC[3], pD[3];

            cleaned->GetPoint(ptIds->GetId(0), pA);
            cleaned->GetPoint(ptIds->GetId(1), pB);
            cleaned->GetPoint(ptIds->GetId(2), pC);
            cleaned->GetPoint(ptIds->GetId(3), pD);

            double det = vtkMath::Determinant3x3(pB[0]-pA[0], pC[0]-pA[0], pD[0]-pA[0],
                pB[1]-pA[1], pC[1]-pA[1], pD[1]-pA[1],
                pB[2]-pA[2], pC[2]-pA[2], pD[2]-pA[2]);

            if (std::abs(det) < 1e-10) {
                newPd->InsertNextCell(VTK_POLYGON, ptIds);
                cellIds->InsertNextValue(cellId);
            } else {
                const vtkIdType cellA[] = {ptIds->GetId(0), ptIds->GetId(1), ptIds->GetId(2)};
                const vtkIdType cellB[] = {ptIds->GetId(2), ptIds->GetId(3), ptIds->GetId(0)};

                newPd->InsertNextCell(VTK_TRIANGLE, 3, cellA);
                cellIds->InsertNextValue(cellId);

                newPd->InsertNextCell(VTK_TRIANGLE, 3, cellB);
                cellIds->InsertNextValue(cellId);
            }
        }
    }

    cellItr->Delete();

    newPd->GetCellData()->SetScalars(cellIds);
    newPd->Squeeze();

    return newPd;
}

Contact::Contact (vtkPolyData *newPdA, vtkPolyData *newPdB) : newPdA(newPdA), newPdB(newPdB) {
    pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetDataTypeToDouble();

    lines = vtkSmartPointer<vtkPolyData>::New();
    lines->SetPoints(pts);
    lines->Allocate(1000);

    contA = vtkSmartPointer<vtkIdTypeArray>::New();
    contB = vtkSmartPointer<vtkIdTypeArray>::New();

    contA->Allocate(1000);
    contB->Allocate(1000);

    contA->SetName("cA");
    contB->SetName("cB");

    sourcesA = vtkSmartPointer<vtkIdTypeArray>::New();
    sourcesB = vtkSmartPointer<vtkIdTypeArray>::New();

    sourcesA->Allocate(1000);
    sourcesB->Allocate(1000);

    sourcesA->SetNumberOfComponents(2);
    sourcesB->SetNumberOfComponents(2);

    sourcesA->SetName("sourcesA");
    sourcesB->SetName("sourcesB");

    touchesEdgesA = false;
    touchesEdgesB = false;

    GetNonManifoldEdges(newPdA, edgesA);
    GetNonManifoldEdges(newPdB, edgesB);
}

vtkSmartPointer<vtkPolyData> Contact::GetLines () {
    auto treeA = vtkSmartPointer<vtkOBBTree>::New();
    treeA->SetDataSet(newPdA);
    treeA->BuildLocator();

    auto treeB = vtkSmartPointer<vtkOBBTree>::New();
    treeB->SetDataSet(newPdB);
    treeB->BuildLocator();

    auto matrix = vtkSmartPointer<vtkMatrix4x4>::New();

    treeA->IntersectWithOBBTree(treeB, matrix, InterNodes, this);

    if (touchesEdgesA || touchesEdgesB) {
        throw std::runtime_error("Intersection goes through non-manifold edges.");
    }

    lines->GetCellData()->AddArray(contA);
    lines->GetCellData()->AddArray(contB);

    lines->GetCellData()->AddArray(sourcesA);
    lines->GetCellData()->AddArray(sourcesB);

    lines->RemoveDeletedCells();
    lines->Squeeze();

    auto clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInputData(lines);
    clean->ToleranceIsAbsoluteOn();
    clean->SetAbsoluteTolerance(1e-5);
    clean->Update();

    auto cleaned = clean->GetOutput();

    vtkCellIterator *cellItr = cleaned->NewCellIterator();

    vtkIdType cellId;

    for (cellItr->InitTraversal(); !cellItr->IsDoneWithTraversal(); cellItr->GoToNextCell()) {
        cellId = cellItr->GetCellId();

        if (cellItr->GetCellType() != VTK_LINE) {
            cleaned->DeleteCell(cellId);
        }
    }

    cleaned->RemoveDeletedCells();

    return cleaned;
}

void Contact::GetNonManifoldEdges (vtkPolyData *pd, NonManifoldEdgesType &edges) {

    vtkCellIterator *cellItr = pd->NewCellIterator();

    vtkIdType cellId;
    vtkIdList *ptIds;

    auto neigs = vtkSmartPointer<vtkIdList>::New();

    vtkIdType i, j;

    vtkIdType idA, idB;

    for (cellItr->InitTraversal(); !cellItr->IsDoneWithTraversal(); cellItr->GoToNextCell()) {
        cellId = cellItr->GetCellId();
        ptIds = cellItr->GetPointIds();

        for (i = 0; i < ptIds->GetNumberOfIds(); i++) {
            j = i+1;

            if (j == ptIds->GetNumberOfIds()) {
                j = 0;
            }

            idA = ptIds->GetId(i);
            idB = ptIds->GetId(j);

            pd->GetCellEdgeNeighbors(cellId, idA, idB, neigs);

            if (neigs->GetNumberOfIds() > 1) {
                edges.emplace(idA, idB);
                edges.emplace(idB, idA);
            }
        }
    }

    cellItr->Delete();
}

void Contact::InterEdgeLine (InterPtsType &interPts, const Point3d &pA, const Point3d &pB, Src src) {
    // schnitt mit x-achse

    double v[3];
    Point3d::GetVec(pA, pB, v);

    double w[] = {-v[1], v[0]};

    double c = w[0]*pA.x+w[1]*pA.y;

    if (std::abs(w[0]) < 1e-10) {
        if (std::abs(pA.y) < 1e-5 && std::abs(pB.y) < 1e-5) {
            // interPts.emplace_back(pA.x, 0, 0, pA.x, pA.id, pB.id, End::A, src, PointSrc::Copied);
            // interPts.emplace_back(pB.x, 0, 0, pB.x, pA.id, pB.id, End::B, src, PointSrc::Copied);
        }
    } else {
        double yA = pA.y-1e-6*v[1];
        double yB = pB.y+1e-6*v[1];

        if (std::signbit(yA) != std::signbit(yB)) {
            double x = c/w[0];

            double sA[] = {pA.x-x, pA.y};
            double sB[] = {pB.x-x, pB.y};

            double dA = sA[0]*sA[0]+sA[1]*sA[1];
            double dB = sB[0]*sB[0]+sB[1]*sB[1];

            End end = dA < 1e-12 ? End::A : (dB < 1e-12 ? End::B : End::None);

            interPts.emplace_back(x, 0, 0, x, pA.id, pB.id, end, src, PointSrc::Calculated);
        }
    }

}

bool Contact::InterPolyLine (InterPtsType &interPts, vtkPolyData *pd, const Base2 &base, const Poly &poly, Src src) {

#if (defined(_debA) && defined(_debB))
    if (_idA == _debA && _idB == _debB) {
        std::cout << "InterPolyLine()" << std::endl;
    }
#endif

    Poly::const_iterator itrA, itrB;

    double q[3];

    for (itrA = poly.begin(); itrA != poly.end(); ++itrA) {
        itrB = itrA+1;

        if (itrB == poly.end()) {
            itrB = poly.begin();
        }

        InterEdgeLine(interPts, *itrA, *itrB, src); // schnitt mit x-achse
    }

    for (auto& p : interPts) {
        base.BackTransform(p.pt, q);

        std::copy_n(q, 3, p.pt);
    }

    std::sort(interPts.begin(), interPts.end(), [](auto &a, auto &b) { return a.t < b.t; });

#if (defined(_debA) && defined(_debB))
    for (auto &p : interPts) {
        std::cout << p << '\n';
    }
#endif

    return true;
}

void Contact::InterPolys (vtkIdType idA, vtkIdType idB) {

#if (defined(_debA) && defined(_debB))
    _idA = idA; _idB = idB;

    if (_idA == _debA && _idB == _debB) {
        std::cout << "InterPolys(" << idA << ", " << idB << ")" << std::endl;
    }
#endif

    vtkIdType numA, numB;
    const vtkIdType *ptsA, *ptsB;

    newPdA->GetCellPoints(idA, numA, ptsA);
    newPdB->GetCellPoints(idB, numB, ptsB);

    Poly polyA, polyB;
    GetPoly(newPdA->GetPoints(), numA, ptsA, polyA);
    GetPoly(newPdB->GetPoints(), numB, ptsB, polyB);

    double nA[3], nB[3], r[3], s[3];

    ComputeNormal(polyA, nA);
    ComputeNormal(polyB, nB);

#if (defined(_debA) && defined(_debB))
    if (_idA == _debA && _idB == _debB) {
        std::cout << "nA [" << nA[0] << ", " << nA[1] << ", " << nA[2] << "]" << std::endl;
        std::cout << "nB [" << nB[0] << ", " << nB[1] << ", " << nB[2] << "]" << std::endl;
    }
#endif

    if (vtkMath::Dot(nA, nB) > .9999999999) {
        return;
    }

    double ptA[3], ptB[3];

    newPdA->GetPoint(ptsA[0], ptA);
    newPdB->GetPoint(ptsB[0], ptB);

    double dA = vtkMath::Dot(nA, ptA);
    double dB = vtkMath::Dot(nB, ptB);

    vtkMath::Cross(nA, nB, r);
    vtkMath::Normalize(r);

    std::array<std::tuple<int, int, int, double>, 3> dets {
        std::make_tuple(0, 1, 2, nA[0]*nB[1]-nB[0]*nA[1]),
        std::make_tuple(0, 2, 1, nA[0]*nB[2]-nB[0]*nA[2]),
        std::make_tuple(1, 2, 0, nA[1]*nB[2]-nB[1]*nA[2])
    };

    const auto& [i, j, k, det] = *std::max_element(dets.begin(), dets.end(), [](auto &a, auto &b) { return std::abs(std::get<3>(a)) < std::abs(std::get<3>(b)); });

    s[i] = (dA*nB[j]-dB*nA[j])/det;
    s[j] = (dB*nA[i]-dA*nB[i])/det;
    s[k] = 0;

#if (defined(_debA) && defined(_debB))
    if (_idA == _debA && _idB == _debB) {
        std::cout << "det " << det << std::endl;
        std::cout << "r [" << r[0] << ", " << r[1] << ", " << r[2] << "]" << std::endl;
        std::cout << "s [" << s[0] << ", " << s[1] << ", " << s[2] << "]" << std::endl;
    }
#endif

    Base2 baseA(s, r, nA);
    Base2 baseB(s, r, nB);

    Poly _polyA, _polyB;

    FlattenPoly2(polyA, _polyA, baseA);
    FlattenPoly2(polyB, _polyB, baseB);

#if (defined(_debA) && defined(_debB))
    if (_idA == _debA && _idB == _debB) {
        std::cout << "_polyA {";
        for (auto& p : _polyA) {
            std::cout << p << ", ";
        }
        std::cout << std::endl;

        std::cout << "_polyB {";
        for (auto& p : _polyB) {
            std::cout << p << ", ";
        }
        std::cout << std::endl;
    }
#endif

    InterPtsType intersPtsA, intersPtsB;

    if (!InterPolyLine(intersPtsA, newPdA, baseA, _polyA, Src::A)) {
        throw std::runtime_error("Found invalid intersection points.");
    }

    if (!InterPolyLine(intersPtsB, newPdB, baseB, _polyB, Src::B)) {
        throw std::runtime_error("Found invalid intersection points.");
    }

    if (!CheckInters(intersPtsA, newPdA)) {
        std::stringstream ss;
        ss << "polys (" << idA << ", " << idB << "): Intersection points do not lie on the edges.";

        if (std::find_if(_polyA.begin(), _polyA.end(), [](auto &p) { return std::abs(p.z) > 1e-6; }) != _polyA.end()) {
            ss << " Normal is inaccurate.";
        }

        throw std::runtime_error(ss.str());
    }

    if (!CheckInters(intersPtsB, newPdB)) {
        std::stringstream ss;
        ss << "polys (" << idA << ", " << idB << "): Intersection points do not lie on the edges.";

        if (std::find_if(_polyB.begin(), _polyB.end(), [](auto &p) { return std::abs(p.z) > 1e-6; }) != _polyB.end()) {
            ss << " Normal is inaccurate.";
        }

        throw std::runtime_error(ss.str());
    }

    if ((intersPtsA.size() & 1) == 0 && (intersPtsB.size() & 1) == 0) {
        AddContactLines(intersPtsA, intersPtsB, idA, idB);
    }

}

bool Contact::CheckInters (const InterPtsType &interPts, vtkPolyData *pd) {
#if (defined(_debA) && defined(_debB))
    if (_idA == _debA && _idB == _debB) {
        std::cout << "CheckInters()" << std::endl;
    }
#endif


    double ptA[3],
        ptB[3],
        v[3],
        w[3],
        k,
        l,
        alpha,
        d;

    for (auto &p : interPts) {

#if (defined(_debA) && defined(_debB))
        if (_idA == _debA && _idB == _debB) {
            std::cout << p << std::endl;
        }
#endif

        pd->GetPoint(p.edge.f, ptA);
        pd->GetPoint(p.edge.g, ptB);

        vtkMath::Subtract(ptA, ptB, v);
        vtkMath::Normalize(v);
        vtkMath::Subtract(ptA, p.pt, w);

        k = vtkMath::Norm(w);
        l = vtkMath::Dot(v, w);
        alpha = std::acos(l/k);

#if (defined(_debA) && defined(_debB))
        if (_idA == _debA && _idB == _debB) {
            std::cout << "alpha " << alpha << std::endl;
        }
#endif

        if (std::isnan(alpha)) {
            continue;
        }

        d = std::sin(alpha)*k;

#if (defined(_debA) && defined(_debB))
        if (_idA == _debA && _idB == _debB) {
            std::cout << "d " << d << std::endl;
        }
#endif

        if (d < 1e-5) {
            continue;
        }

        return false;

    }

    return true;

}

void Contact::OverlapLines (OverlapsType &overlaps, InterPtsType &intersA, InterPtsType &intersB) {

    auto Add = [](InterPt &a, InterPt &b, InterPt &c, InterPt &d) {
        a.Merge(c);
        b.Merge(d);

        return std::make_tuple(a, b);
    };

    InterPtsType::iterator itr, itr2;

    for (itr = intersA.begin(); itr != intersA.end(); itr += 2) {
        for (itr2 = intersB.begin(); itr2 != intersB.end(); itr2 += 2) {
            if (itr->t <= itr2->t && (itr+1)->t > itr2->t) {
                if ((itr2+1)->t < (itr+1)->t) {
                    overlaps.push_back(Add(*itr2, *(itr2+1), *itr, *(itr+1)));
                } else {
                    overlaps.push_back(Add(*itr2, *(itr+1), *itr, *(itr2+1)));
                }
            } else if (itr2->t <= itr->t && (itr2+1)->t > itr->t) {
                if ((itr+1)->t < (itr2+1)->t) {
                    overlaps.push_back(Add(*itr, *(itr+1), *itr2, *(itr2+1)));
                } else {
                    overlaps.push_back(Add(*itr, *(itr2+1), *itr2, *(itr+1)));
                }
            }
        }
    }

}

void Contact::AddContactLines (InterPtsType &intersA, InterPtsType &intersB, vtkIdType idA, vtkIdType idB) {

    if (intersA.size() == 0 || intersB.size() == 0) {
        return;
    }

    OverlapsType overlaps;
    OverlapLines(overlaps, intersA, intersB);

    OverlapsType::const_iterator itr;

    for (itr = overlaps.begin(); itr != overlaps.end(); ++itr) {
        auto &f = std::get<0>(*itr);
        auto &s = std::get<1>(*itr);

        if (f.src == Src::A) {
            if (edgesA.count(f.edge) == 1) {
                touchesEdgesA = true;
            }
        }

        if (s.src == Src::A) {
            if (edgesA.count(s.edge) == 1) {
                touchesEdgesA = true;
            }
        }

        if (f.src == Src::B) {
            if (edgesB.count(f.edge) == 1) {
                touchesEdgesB = true;
            }
        }

        if (s.src == Src::B) {
            if (edgesB.count(s.edge) == 1) {
                touchesEdgesB = true;
            }
        }

        vtkIdList *linePts = vtkIdList::New();

        linePts->InsertNextId(pts->InsertNextPoint(f.pt));
        linePts->InsertNextId(pts->InsertNextPoint(s.pt));

        lines->InsertNextCell(VTK_LINE, linePts);

        linePts->Delete();

        const vtkIdType tupleA[] = {f.srcA, s.srcA};
        const vtkIdType tupleB[] = {f.srcB, s.srcB};

        sourcesA->InsertNextTypedTuple(tupleA);
        sourcesB->InsertNextTypedTuple(tupleB);

        contA->InsertNextValue(idA);
        contB->InsertNextValue(idB);
    }

}

int Contact::InterNodes (vtkOBBNode *nodeA, vtkOBBNode *nodeB, vtkMatrix4x4 *vtkNotUsed(matrix), void *ptr) {
    auto _this = reinterpret_cast<Contact*>(ptr);

    vtkIdList *cellsA = nodeA->Cells;
    vtkIdList *cellsB = nodeB->Cells;

    vtkIdType numCellsA = cellsA->GetNumberOfIds();
    vtkIdType numCellsB = cellsB->GetNumberOfIds();

    vtkIdType i, j, cellA, cellB;

    for (i = 0; i < numCellsA; i++) {
        cellA = cellsA->GetId(i);

        for (j = 0; j < numCellsB; j++) {
            cellB = cellsB->GetId(j);

            _this->InterPolys(cellA, cellB);
        }
    }

    return 0;
}
