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
#include "Optimize.h"

#include <vtkCleanPolyData.h>
#include <vtkCellIterator.h>
#include <vtkCellData.h>
#include <vtkTriangleStrip.h>
#include <vtkArrayIteratorTemplate.h>

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
    if (newPdA->GetCellData()->GetScalars("OrigCellIds") == nullptr) {
        throw std::invalid_argument("OrigCellIds is missing.");
    }

    if (newPdB->GetCellData()->GetScalars("OrigCellIds") == nullptr) {
        throw std::invalid_argument("OrigCellIds is missing.");
    }

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

    lines->GetCellData()->AddArray(contA);
    lines->GetCellData()->AddArray(contB);

    lines->GetCellData()->AddArray(sourcesA);
    lines->GetCellData()->AddArray(sourcesB);

    IntersectReplacements();

    if (touchesEdgesA || touchesEdgesB) {
        throw std::runtime_error("Intersection goes through non-manifold edges.");
    }

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

    vtkIdType i, j, k;

    vtkIdType idA, idB;

    vtkIdType n;

    vtkIdType num;
    const vtkIdType *pts;

    vtkIdType a, b;

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

                n = 0;

                for (k = 0; k < neigs->GetNumberOfIds(); k++) {
                    pd->GetCellPoints(neigs->GetId(k), num, pts);

                    for (a = 0; a < num; a++) {
                        b = a+1;

                        if (b == num) {
                            b = 0;
                        }

                        if (pts[a] == idB && pts[b] == idA) {
                            n++;
                        }
                    }
                }

                if (n > 1) {
                    edges.emplace(idA, idB);
                    edges.emplace(idB, idA);
                }
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
            interPts.emplace_back(pA.x, 0, 0, pA.x, pA.id, pB.id, End::A, src, PointSrc::Copied);
            interPts.emplace_back(pB.x, 0, 0, pB.x, pA.id, pB.id, End::B, src, PointSrc::Copied);
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

bool Contact::InterPolyLine (InterPtsType &interPts, const Base2 &base, const Poly &poly, Src src) {

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

    struct Cmp {
        const double tol = 1e-5;

        bool operator() (const InterPt &lhs, const InterPt &rhs) const {
            if (lhs.pointSrc == PointSrc::Copied && rhs.pointSrc == PointSrc::Copied) {
                if (lhs.GetEnd() == rhs.GetEnd()) {
                    return false;
                }
            } else if (lhs.pointSrc == PointSrc::Copied) {
                const vtkIdType end = lhs.GetEnd();

                if (end == rhs.edge.f || end == rhs.edge.g) {
                    return false;
                }
            } else if (rhs.pointSrc == PointSrc::Copied) {
                const vtkIdType end = rhs.GetEnd();

                if (end == lhs.edge.f || end == lhs.edge.g) {
                    return false;
                }
            }

            const double d = lhs.t-rhs.t;

            if (std::abs(d) < tol) {
                return false;
            } else {
                return d < 0;
            }
        }
    };

    std::map<InterPt, InterPtsType, Cmp> grouped;

    std::vector<InterPtsType> sortedPts;

    for (auto &p : interPts) {
        grouped[p].push_back(p);
    }

    for (auto& [k, v] : grouped) {
        std::sort(v.begin(), v.end(), [](const InterPt &lhs, const InterPt &rhs) { return lhs.pointSrc > rhs.pointSrc; });

        sortedPts.push_back(v);
    }

    std::map<vtkIdType, double> allEnds;

    decltype(sortedPts)::iterator itr;

    for (itr = sortedPts.begin(); itr != sortedPts.end(); ++itr) {
        if (itr == sortedPts.begin()) {
            if (itr->size() == 2) {
                itr->pop_back();
            }
        } else if (itr == sortedPts.end()-1) {
            if (itr->size() == 2) {
                itr->pop_back();
            }
        } else if (itr->size() == 1 && itr->back().end != End::None) {
            itr->push_back(itr->back());
        }

        if (itr->back().end != End::None) {
            auto p = itr->back();

            allEnds.emplace(p.end == End::A ? p.edge.f : p.edge.g, p.t);
        }

        // hier schneidet sich das polygon selbst

        if (itr->size() > 2) {
            return false;
        }

        if (itr->size() == 2) {
            // schnitt durch kongruente kanten

            if (itr->back().end == End::None) {
                return false;
            }

            auto &edgeA = itr->front().edge;
            auto &edgeB = itr->back().edge;

            if (edgeA.f == edgeB.g && edgeB.f == edgeA.g) {
                return false;
            }

        }

    }

    std::map<vtkIdType, std::reference_wrapper<const Point3d>> pts;

    for (auto &p : poly) {
        pts.emplace(p.id, p);
    }

    vtkIdType ind;

    vtkIdType prev, next;

    for (itr = sortedPts.begin(); itr != sortedPts.end(); ++itr) {
        auto p = itr->back();

        if (p.end == End::None) {
            continue;
        }

        if (p.end == End::A) {
            ind = p.edge.f;

            next = p.edge.g;

            auto _itr = std::find_if(poly.begin(), poly.end(), [&ind](auto &p) { return p.id == ind; });

            if (_itr == poly.begin()) {
                prev = poly.back().id;
            } else {
                prev = std::prev(_itr)->id;
            }

        } else {
            ind = p.edge.g;

            prev = p.edge.f;

            auto _itr = std::find_if(poly.begin(), poly.end(), [&ind](auto &p) { return p.id == ind; });

            if (_itr == poly.end()-1) {
                next = poly.front().id;
            } else {
                next = std::next(_itr)->id;
            }
        }

        if (itr->size() == 2) {
            if (allEnds.count(prev) == 0 && allEnds.count(next) == 1) {
                const Point3d &q = pts.at(prev);

                if ((allEnds.at(next) < p.t && q.y < 0) || (allEnds.at(next) > p.t && q.y > 0)) {
                    itr->pop_back();
                }

            } else if (allEnds.count(prev) == 1 && allEnds.count(next) == 0) {
                const Point3d &q = pts.at(next);

                if ((allEnds.at(prev) < p.t && q.y > 0) || (allEnds.at(prev) > p.t && q.y < 0)) {
                    itr->pop_back();
                }
            }
        }

        if (allEnds.count(prev) == 0 && allEnds.count(next) == 0) {
            const Point3d &a = pts.at(prev);
            const Point3d &b = pts.at(next);

            if (std::signbit(a.y) != std::signbit(b.y)) {
                if (itr->size() == 2) {
                    itr->pop_back();
                }
            } else {
                if ((a.x > b.x) == std::signbit(a.y)) {
                    itr->clear();
                }
            }

        }

    }

    InterPtsType _interPts;

    for (const auto &pts : sortedPts) {
        std::copy(pts.begin(), pts.end(), std::back_inserter(_interPts));
    }

    interPts.swap(_interPts);

#if (defined(_debA) && defined(_debB))
    for (auto &p : interPts) {
        std::cout << p << std::endl;
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

    Poly transA, transB;

    FlattenPoly2(polyA, transA, baseA);
    FlattenPoly2(polyB, transB, baseB);

#if (defined(_debA) && defined(_debB))
    if (_idA == _debA && _idB == _debB) {
        std::cout << "transA {";
        for (auto& p : transA) {
            std::cout << p << ", ";
        }
        std::cout << std::endl;

        std::cout << "transB {";
        for (auto& p : transB) {
            std::cout << p << ", ";
        }
        std::cout << std::endl;
    }
#endif

    bool isPlanarA = std::find_if(transA.begin(), transA.end(), [](auto &p) { return std::abs(p.z) > 1e-6; }) == transA.end();
    bool isPlanarB = std::find_if(transB.begin(), transB.end(), [](auto &p) { return std::abs(p.z) > 1e-6; }) == transB.end();

    bool hasReplA = replsA.count(idA) == 1;
    bool hasReplB = replsB.count(idB) == 1;

    if (!isPlanarA && !hasReplA && newPdA->GetCellType(idA) != VTK_TRIANGLE) {
        auto newIds = PreventEqualCaptPoints::TriangluteCell(newPdA, idA, {});
        replsA.emplace(idA, newIds);
        hasReplA = true;
    }

    if (!isPlanarB && !hasReplB && newPdB->GetCellType(idB) != VTK_TRIANGLE) {
        auto newIds = PreventEqualCaptPoints::TriangluteCell(newPdB, idB, {});
        replsB.emplace(idB, newIds);
        hasReplB = true;
    }

    if (hasReplA || hasReplB) {
        pairs.emplace_back(idA, idB);
        return;
    }

    InterPtsType intersPtsA, intersPtsB;

    if (!InterPolyLine(intersPtsA, baseA, transA, Src::A)) {
        throw std::runtime_error("Found invalid intersection points.");
    }

    if (!InterPolyLine(intersPtsB, baseB, transB, Src::B)) {
        throw std::runtime_error("Found invalid intersection points.");
    }

    if (!CheckInters(intersPtsA, newPdA)) {
        std::stringstream ss;
        ss << "Intersection points do not lie on the edges (cells " << idA << ", " << idB << ").";

        throw std::runtime_error(ss.str());
    }

    if (!CheckInters(intersPtsB, newPdB)) {
        std::stringstream ss;
        ss << "Intersection points do not lie on the edges (cells " << idA << ", " << idB << ").";

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

void Contact::IntersectReplacements () {
    {
        std::vector<std::tuple<vtkIdType, vtkIdType, vtkIdType>> invalid;

        vtkIdType i;

        auto iterA = vtkArrayIteratorTemplate<vtkIdType>::New();
        iterA->Initialize(contA);

        auto iterB = vtkArrayIteratorTemplate<vtkIdType>::New();
        iterB->Initialize(contB);

        for (i = 0; i < iterA->GetNumberOfValues(); i++) {
            if (replsA.count(iterA->GetValue(i)) == 1 || replsB.count(iterB->GetValue(i)) == 1) {
                invalid.emplace_back(i, iterA->GetValue(i), iterB->GetValue(i));
            }
        }

        for (auto& [lineId, idA, idB] : invalid) {
            lines->DeleteCell(lineId);

            pairs.emplace_back(idA, idB);
        }

        lines->RemoveDeletedCells();

        contA = vtkIdTypeArray::SafeDownCast(lines->GetCellData()->GetScalars("cA"));
        contB = vtkIdTypeArray::SafeDownCast(lines->GetCellData()->GetScalars("cB"));

        sourcesA = vtkIdTypeArray::SafeDownCast(lines->GetCellData()->GetScalars("sourcesA"));
        sourcesB = vtkIdTypeArray::SafeDownCast(lines->GetCellData()->GetScalars("sourcesB"));
    }

    for (auto &p : pairs) {
        auto itrA = replsA.find(p.f);
        auto itrB = replsB.find(p.g);

        IdsType cellsA, cellsB;

        if (itrA == replsA.end()) {
            cellsA.push_back(p.f);
        } else {
            auto ids = itrA->second;
            std::copy(ids.begin(), ids.end(), std::back_inserter(cellsA));

            newPdA->DeleteCell(p.f);
        }

        if (itrB == replsB.end()) {
            cellsB.push_back(p.g);
        } else {
            auto ids = itrB->second;
            std::copy(ids.begin(), ids.end(), std::back_inserter(cellsB));

            newPdB->DeleteCell(p.g);
        }

        for (auto &idA : cellsA) {
            for (auto &idB : cellsB) {
                InterPolys(idA, idB);
            }
        }

    }

    // contA und contB aktualisieren

    auto oldCellIdsA = vtkSmartPointer<vtkIdTypeArray>::New();
    auto oldCellIdsB = vtkSmartPointer<vtkIdTypeArray>::New();

    oldCellIdsA->SetName("OldCellIds");
    oldCellIdsB->SetName("OldCellIds");

    vtkIdType numCellsA = newPdA->GetNumberOfCells();
    vtkIdType numCellsB = newPdB->GetNumberOfCells();

    oldCellIdsA->SetNumberOfValues(numCellsA);
    oldCellIdsB->SetNumberOfValues(numCellsB);

    vtkIdType i;

    for (i = 0; i < numCellsA; i++) {
        oldCellIdsA->SetValue(i, i);
    }

    for (i = 0; i < numCellsB; i++) {
        oldCellIdsB->SetValue(i, i);
    }

    newPdA->GetCellData()->AddArray(oldCellIdsA);
    newPdB->GetCellData()->AddArray(oldCellIdsB);

    newPdA->RemoveDeletedCells();
    newPdB->RemoveDeletedCells();

    numCellsA = newPdA->GetNumberOfCells();
    numCellsB = newPdB->GetNumberOfCells();

    oldCellIdsA = vtkIdTypeArray::SafeDownCast(newPdA->GetCellData()->GetScalars("OldCellIds"));
    oldCellIdsB = vtkIdTypeArray::SafeDownCast(newPdB->GetCellData()->GetScalars("OldCellIds"));

    std::map<vtkIdType, vtkIdType> newCellIdsA, newCellIdsB;

    auto iterA = vtkArrayIteratorTemplate<vtkIdType>::New();
    iterA->Initialize(oldCellIdsA);

    for (i = 0; i < numCellsA; i++) {
        newCellIdsA.emplace(iterA->GetValue(i), i);
    }

    auto iterB = vtkArrayIteratorTemplate<vtkIdType>::New();
    iterB->Initialize(oldCellIdsB);

    for (i = 0; i < numCellsB; i++) {
        newCellIdsB.emplace(iterB->GetValue(i), i);
    }

    vtkIdType numLines = lines->GetNumberOfCells();

    try {

        auto _iterA = vtkArrayIteratorTemplate<vtkIdType>::New();
        _iterA->Initialize(contA);

        for (i = 0; i < numLines; i++) {
            _iterA->SetValue(i, newCellIdsA.at(_iterA->GetValue(i)));
        }

        auto _iterB = vtkArrayIteratorTemplate<vtkIdType>::New();
        _iterB->Initialize(contB);

        for (i = 0; i < numLines; i++) {
            _iterB->SetValue(i, newCellIdsB.at(_iterB->GetValue(i)));
        }

    } catch (const std::out_of_range &e) {
        throw std::runtime_error("");
    }

    newPdA->GetCellData()->RemoveArray("OldCellIds");
    newPdB->GetCellData()->RemoveArray("OldCellIds");

}
