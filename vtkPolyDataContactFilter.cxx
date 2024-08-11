/*
Copyright 2012-2024 Ronald Römer

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
#include <map>
#include <set>
#include <algorithm>

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkObjectFactory.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkPolyData.h>
#include <vtkOBBTree.h>
#include <vtkMatrix4x4.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkMath.h>
#include <vtkIdTypeArray.h>
#include <vtkCellData.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleStrip.h>
#include <vtkSmartPointer.h>
#include <vtkFeatureEdges.h>
#include <vtkCellIterator.h>
#include <vtkCellArrayIterator.h>
#include <vtkCellArray.h>

#include "vtkPolyDataContactFilter.h"
#include "Utilities.h"

// #define _debA 18340
// #define _debB 1659

#if (defined(_debA) && defined(_debB))
vtkIdType _idA, _idB;
#endif

vtkStandardNewMacro(vtkPolyDataContactFilter);

vtkPolyDataContactFilter::vtkPolyDataContactFilter () {
    SetNumberOfInputPorts(2);
    SetNumberOfOutputPorts(3);

    invalidA = false;
    invalidB = false;

    aborted = false;
}

vtkPolyDataContactFilter::~vtkPolyDataContactFilter () {
    // nix mehr
}

int vtkPolyDataContactFilter::RequestData (vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

    if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA())) {

        vtkInformation *inInfoA = inputVector[0]->GetInformationObject(0);
        vtkInformation *inInfoB = inputVector[1]->GetInformationObject(0);

        vtkPolyData *_pdA = vtkPolyData::SafeDownCast(inInfoA->Get(vtkDataObject::DATA_OBJECT()));
        vtkPolyData *_pdB = vtkPolyData::SafeDownCast(inInfoB->Get(vtkDataObject::DATA_OBJECT()));

        vtkInformation *outInfoA = outputVector->GetInformationObject(0);
        vtkInformation *outInfoB = outputVector->GetInformationObject(1);
        vtkInformation *outInfoC = outputVector->GetInformationObject(2);

        vtkPolyData *resultA = vtkPolyData::SafeDownCast(outInfoA->Get(vtkDataObject::DATA_OBJECT()));
        vtkPolyData *resultB = vtkPolyData::SafeDownCast(outInfoB->Get(vtkDataObject::DATA_OBJECT()));
        vtkPolyData *resultC = vtkPolyData::SafeDownCast(outInfoC->Get(vtkDataObject::DATA_OBJECT()));

        // contLines anlegen

        contLines = vtkSmartPointer<vtkPolyData>::New();
        contLines->Allocate(1000);

        auto contPts = vtkSmartPointer<vtkPoints>::New();
        contPts->SetDataTypeToDouble();

        contLines->SetPoints(contPts);

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

        // durchführung der aufgabe

        newPdA = vtkSmartPointer<vtkPolyData>::New();
        newPdB = vtkSmartPointer<vtkPolyData>::New();

        CopyPolyData(_pdA, newPdA);
        CopyPolyData(_pdB, newPdB);

        GetInvalidEdges(newPdA, edgesA);
        GetInvalidEdges(newPdB, edgesB);

        auto obbA = vtkSmartPointer<vtkOBBTree>::New();
        obbA->SetDataSet(newPdA);
        obbA->SetNumberOfCellsPerNode(10);
        obbA->BuildLocator();

        auto obbB = vtkSmartPointer<vtkOBBTree>::New();
        obbB->SetDataSet(newPdB);
        obbB->SetNumberOfCellsPerNode(10);
        obbB->BuildLocator();

        auto mat = vtkSmartPointer<vtkMatrix4x4>::New();

        obbA->IntersectWithOBBTree(obbB, mat, InterOBBNodes, this);

        if (aborted) {
            vtkErrorMacro("Bad shaped cells detected.");
            return 1;
        }

        if (invalidA) {
            vtkErrorMacro("First input has non-manifold edges.");
            return 1;
        }

        if (invalidB) {
            vtkErrorMacro("Second input has non-manifold edges.");
            return 1;
        }

        if (contLines->GetNumberOfCells() == 0) {
            vtkErrorMacro("There is no contact.");
            return 1;
        }

        contLines->GetCellData()->AddArray(contA);
        contLines->GetCellData()->AddArray(contB);

        contLines->GetCellData()->AddArray(sourcesA);
        contLines->GetCellData()->AddArray(sourcesB);

        contLines->RemoveDeletedCells();

        contLines->Squeeze();

        auto clean = vtkSmartPointer<vtkCleanPolyData>::New();
        clean->SetInputData(contLines);
        clean->ToleranceIsAbsoluteOn();
        clean->SetAbsoluteTolerance(1e-5);
        clean->Update();

        resultA->DeepCopy(clean->GetOutput());

        vtkIdType i, numCellsA = resultA->GetNumberOfCells();

        for (i = 0; i < numCellsA; i++) {
            if (resultA->GetCellType(i) != VTK_LINE) {
                resultA->DeleteCell(i);
            }
        }

        resultA->RemoveDeletedCells();

        resultB->DeepCopy(newPdA);
        resultC->DeepCopy(newPdB);
    }

    return 1;

}

void vtkPolyDataContactFilter::CopyPolyData (vtkPolyData *pd, vtkPolyData *newPd) {
    auto newPoints = vtkSmartPointer<vtkPoints>::New();

    newPoints->SetDataType(pd->GetPoints()->GetDataType());
    newPoints->DeepCopy(pd->GetPoints());

    newPd->SetPoints(newPoints);

    auto newPolys = vtkSmartPointer<vtkCellArray>::New();

    newPolys->Allocate(pd->GetNumberOfPolys());

    newPd->SetPolys(newPolys);

    auto cellIds = vtkSmartPointer<vtkIdTypeArray>::New();
    cellIds->Allocate(pd->GetNumberOfPolys());

    vtkCellIterator *cellItr = pd->NewCellIterator();

    vtkIdType cellId;
    vtkIdList *ptIds;

    double p0[3], p1[3], p2[3], p3[3];

    double n[3], d;

    double dA, dB;

    for (cellItr->InitTraversal(); !cellItr->IsDoneWithTraversal(); cellItr->GoToNextCell()) {
        cellId = cellItr->GetCellId();
        ptIds = cellItr->GetPointIds();

        if (cellItr->GetCellType() == VTK_TRIANGLE || cellItr->GetCellType() == VTK_POLYGON) {
            newPd->InsertNextCell(cellItr->GetCellType(), ptIds);

            cellIds->InsertNextValue(cellId);

        } else if (cellItr->GetCellType() == VTK_TRIANGLE_STRIP) {
            auto cells = vtkSmartPointer<vtkCellArray>::New();

            vtkTriangleStrip::DecomposeStrip(cellItr->GetNumberOfPoints(), ptIds->GetPointer(0), cells);

            vtkIdType n;
            const vtkIdType *pts;

            for (cells->InitTraversal(); cells->GetNextCell(n, pts);) {
                if (pts[0] != pts[1] && pts[1] != pts[2] && pts[2] != pts[0]) {
                    newPd->InsertNextCell(VTK_TRIANGLE, n, pts);
                    cellIds->InsertNextValue(cellId);
                }
            }

        } else if (cellItr->GetCellType() == VTK_QUAD) {

            newPoints->GetPoint(ptIds->GetId(0), p0);
            newPoints->GetPoint(ptIds->GetId(1), p1);
            newPoints->GetPoint(ptIds->GetId(2), p2);
            newPoints->GetPoint(ptIds->GetId(3), p3);

            ComputeNormal(newPoints, n, 4, ptIds->GetPointer(0));

            d = vtkMath::Dot(n, p0);

            if (CheckNormal(newPoints, 4, ptIds->GetPointer(0), n, d)) {
                newPd->InsertNextCell(VTK_POLYGON, ptIds);
                cellIds->InsertNextValue(cellId);

            } else {
                dA = vtkMath::Distance2BetweenPoints(p0, p2);
                dB = vtkMath::Distance2BetweenPoints(p1, p3);

                if (dA < dB) {
                    // 0, 2
                    const vtkIdType cellA[] = {ptIds->GetId(0), ptIds->GetId(1), ptIds->GetId(2)};
                    const vtkIdType cellB[] = {ptIds->GetId(2), ptIds->GetId(3), ptIds->GetId(0)};

                    newPd->InsertNextCell(VTK_TRIANGLE, 3, cellA);
                    cellIds->InsertNextValue(cellId);

                    newPd->InsertNextCell(VTK_TRIANGLE, 3, cellB);
                    cellIds->InsertNextValue(cellId);

                } else {
                    // 1, 3
                    const vtkIdType cellA[] = {ptIds->GetId(1), ptIds->GetId(2), ptIds->GetId(3)};
                    const vtkIdType cellB[] = {ptIds->GetId(3), ptIds->GetId(0), ptIds->GetId(1)};

                    newPd->InsertNextCell(VTK_TRIANGLE, 3, cellA);
                    cellIds->InsertNextValue(cellId);

                    newPd->InsertNextCell(VTK_TRIANGLE, 3, cellB);
                    cellIds->InsertNextValue(cellId);
                }
            }
        }
    }

    cellItr->Delete();

    cellIds->SetName("OrigCellIds");
    newPd->GetCellData()->SetScalars(cellIds);

    newPd->Squeeze();

}

void vtkPolyDataContactFilter::GetInvalidEdges (vtkPolyData *pd, InvalidEdgesType &edges) {
    auto features = vtkSmartPointer<vtkFeatureEdges>::New();
    features->SetInputData(pd);

    features->BoundaryEdgesOff();
    features->FeatureEdgesOff();
    features->ManifoldEdgesOff();
    features->NonManifoldEdgesOn();

    features->Update();

    vtkPolyData *lines = features->GetOutput();

    vtkIdType num;
    const vtkIdType *line;

    double ptA[3], ptB[3];

    vtkIdType idA, idB;

    auto cellsA = vtkSmartPointer<vtkIdList>::New();
    auto cellsB = vtkSmartPointer<vtkIdList>::New();

    auto cell = vtkSmartPointer<vtkIdList>::New();

    vtkIdType indexA, indexB, nextA, nextB, numNeigs;

    auto lineItr = vtk::TakeSmartPointer(lines->GetLines()->NewIterator());

    for (lineItr->GoToFirstCell(); !lineItr->IsDoneWithTraversal(); lineItr->GoToNextCell()) {
        lineItr->GetCurrentCell(num, line);

        lines->GetPoint(line[0], ptA);
        lines->GetPoint(line[1], ptB);

        idA = pd->FindPoint(ptA);
        idB = pd->FindPoint(ptB);

        pd->GetPointCells(idA, cellsA);
        pd->GetPointCells(idB, cellsB);

        cellsA->IntersectWith(cellsB);

        std::set<vtkIdType> cellIds(cellsA->begin(), cellsA->end());

        numNeigs = 0;

        for (auto &id : cellIds) {
            pd->GetCellPoints(id, cell);

            indexA = cell->IsId(idA);
            indexB = cell->IsId(idB);

            nextA = indexA+1;

            if (nextA == cell->GetNumberOfIds()) {
                nextA = 0;
            }

            nextB = indexB+1;

            if (nextB == cell->GetNumberOfIds()) {
                nextB = 0;
            }

            if (cell->GetId(nextA) == idB || cell->GetId(nextB) == idA) {
                numNeigs++;
            }
        }

        if (numNeigs > 2) {
            edges.emplace(idA, idB);
            edges.emplace(idB, idA);
        }

    }
}

void vtkPolyDataContactFilter::InterEdgeLine (InterPtsType &interPts, vtkPolyData *pd, vtkIdType idA, vtkIdType idB, const double *r, const double *ptA, Src src) {
    double eA[3], eB[3];

    pd->GetPoint(idA, eA);
    pd->GetPoint(idB, eB);

    double ptB[3];
    vtkMath::Add(ptA, r, ptB);

    // richtungsvektor der kante bestimmen

    double e[3];
    vtkMath::Subtract(eB, eA, e);
    double l = vtkMath::Normalize(e);

    double p[3];
    vtkMath::Subtract(eA, ptA, p);

    double w = std::abs(vtkMath::Determinant3x3(r, e, p));

    if (w < 1e-4) { // ~89.995deg
        // schnittpunkt ermitteln

        double v[3];
        vtkMath::Cross(r, e, v);

        double n = vtkMath::Norm(v);

        if (n > 1e-4) { // ~0.0057deg

            double s = vtkMath::Determinant3x3(p, r, v)/(n*n);

            if (s > -1e-6 && s < l+1e-6) {

                double t = vtkMath::Determinant3x3(p, e, v)/(n*n);

                InterPt q(ptA[0]+t*r[0], ptA[1]+t*r[1], ptA[2]+t*r[2], t, idA, idB, End::NONE, src, PointSrc::CALCULATED);

                if (s > -1e-6 && s < 1e-6) {
                    q.end = End::BEGIN;
                } else if (s > l-1e-6 && s < l+1e-6) {
                    q.end = End::END;
                }

#if (defined(_debA) && defined(_debB))
                if (_idA == _debA && _idB == _debB) {
                    std::cout << "s " << s << ", " << q << std::endl;
                }

#endif

                interPts.push_back(std::move(q));

            }

        } else {
            // parallel

            double vA[3], vB[3];

            vtkMath::Subtract(eA, ptA, vA);
            vtkMath::Subtract(eB, ptA, vB);

            double dotA = vtkMath::Dot(vA, r),
                dotB = vtkMath::Dot(vB, r);

            double a = vtkMath::Norm(vA),
                b = vtkMath::Norm(vB);

            double angA = std::acos(dotA/a),
                angB = std::acos(dotB/b);

            double dA = std::sin(angA)*a,
                dB = std::sin(angB)*b;

            if ((std::isnan(dA) || dA < 1e-5) && (std::isnan(dB) || dB < 1e-5)) {
                InterPt qA(ptA[0]+dotA*r[0], ptA[1]+dotA*r[1], ptA[2]+dotA*r[2], dotA, idA, idB, End::BEGIN, src, PointSrc::COPIED);
                InterPt qB(ptA[0]+dotB*r[0], ptA[1]+dotB*r[1], ptA[2]+dotB*r[2], dotB, idA, idB, End::END, src, PointSrc::COPIED);

#if (defined(_debA) && defined(_debB))
                if (_idA == _debA && _idB == _debB) {
                    std::cout << "dA " << dA << ", " << qA << std::endl;
                    std::cout << "dB " << dB << ", " << qB << std::endl;
                }
#endif

                interPts.push_back(std::move(qA));
                interPts.push_back(std::move(qB));

            }

        }

    } else {
        // windschief
    }

}

bool vtkPolyDataContactFilter::InterPolyLine (InterPtsType &interPts, vtkPolyData *pd, vtkIdType num, const vtkIdType *poly, const double *r, const double *pt, Src src, const double *n) {

#if (defined(_debA) && defined(_debB))
    if (_idA == _debA && _idB == _debB) {
        std::cout << "InterPolyLine()" << std::endl;
    }
#endif

    std::vector<Pair> edges;
    edges.reserve(static_cast<std::size_t>(num));

    vtkIdType i, j;

    for (i = 0; i < num; i++) {
        j = i+1;

        if (j == num) {
            j = 0;
        }

        const vtkIdType &a = poly[i],
            &b = poly[j];

        vtkPolyDataContactFilter::InterEdgeLine(interPts, pd, a, b, r, pt, src);

        edges.emplace_back(a, b);

    }

    if (interPts.empty()) {
        return true;
    }

    struct Cmp {
        const double tol = 1e-5;

        bool operator() (const InterPt &lhs, const InterPt &rhs) const {
            if (lhs.pointSrc == PointSrc::COPIED && rhs.pointSrc == PointSrc::COPIED) {
                const vtkIdType indA = lhs.getEnd(),
                    indB = rhs.getEnd();

                if (indA == indB) {
                    return false;
                }
            } else if (lhs.pointSrc == PointSrc::COPIED) {
                const vtkIdType ind = lhs.getEnd();

                if (ind == rhs.edge.f || ind == rhs.edge.g) {
                    return false;
                }
            } else if (rhs.pointSrc == PointSrc::COPIED) {
                const vtkIdType ind = rhs.getEnd();

                if (ind == lhs.edge.f || ind == lhs.edge.g) {
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

    std::map<InterPt, InterPtsType, Cmp> paired;

    for (auto &p : interPts) {
        paired[p].push_back(p);
    }

    std::vector<InterPtsType> sortedPts;

    std::transform(paired.begin(), paired.end(), std::back_inserter(sortedPts), [](auto &kv) {
        std::sort(kv.second.begin(), kv.second.end(), [](const InterPt &lhs, const InterPt &rhs) { return lhs.pointSrc > rhs.pointSrc; });

        return kv.second;
    });

    std::map<vtkIdType, double> ends;

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
        } else if (itr->size() == 1 && itr->back().end != End::NONE) {
            itr->push_back(itr->back());
        }

        // validierung

        if (itr->size() > 2) {
            return false;
        }

        if (itr->size() == 2) {
            if (itr->back().end == End::NONE) {
                return false;
            }

            const Pair &edgeA = itr->front().edge,
                &edgeB = itr->back().edge;

            if (edgeA.f == edgeB.g && edgeB.f == edgeA.g) {
                return false;
            }
        }

        const auto &p = itr->back();

        if (p.end != End::NONE) {
            ends.emplace(p.end == End::BEGIN ? p.edge.f : p.edge.g, p.t);
        }

    }

    double s[3], d;

    vtkMath::Cross(n, r, s);
    d = vtkMath::Dot(s, pt);

    double ptA[3], ptB[3], dA, dB, vA[3], vB[3], tA, tB;

    vtkIdType id, prev, next;

    for (itr = sortedPts.begin(); itr != sortedPts.end(); ++itr) {
        const auto &p = itr->back();

        if (p.end != End::NONE) {
            if (p.end == End::BEGIN) {
                id = p.edge.f;
                next = p.edge.g;
                prev = std::find_if(edges.begin(), edges.end(), [&id](const Pair &edge) { return edge.g == id; })->f;
            } else {
                id = p.edge.g;
                prev = p.edge.f;
                next = std::find_if(edges.begin(), edges.end(), [&id](const Pair &edge) { return edge.f == id; })->g;
            }

            if (itr->size() == 2) {
                if (ends.count(next) == 0 && ends.count(prev) == 1) {
                    pd->GetPoint(next, ptA);
                    dA = vtkMath::Dot(s, ptA)-d;

                    if ((p.t > ends.at(prev) && dA > 0) || (p.t < ends.at(prev) && dA < 0)) {
                        // tasche
                        itr->pop_back();
                    }

                    continue;

                } else if (ends.count(next) == 1 && ends.count(prev) == 0) {
                    pd->GetPoint(prev, ptB);
                    dB = vtkMath::Dot(s, ptB)-d;

                    if ((p.t > ends.at(next) && dB < 0) || (p.t < ends.at(next) && dB > 0)) {
                        // tasche
                        itr->pop_back();
                    }

                    continue;

                }
            }

            if (ends.count(prev) == 0 && ends.count(next) == 0) {
                pd->GetPoint(next, ptA);
                pd->GetPoint(prev, ptB);

                dA = vtkMath::Dot(s, ptA)-d;
                dB = vtkMath::Dot(s, ptB)-d;

                if (std::signbit(dA) != std::signbit(dB)) {
                    if (itr->size() == 2) {
                        itr->pop_back();
                    }

                } else {
                    vtkMath::Subtract(ptA, pt, vA);
                    vtkMath::Subtract(ptB, pt, vB);

                    tA = vtkMath::Dot(vA, r);
                    tB = vtkMath::Dot(vB, r);

                    if ((tB > tA) == std::signbit(dA)) {
                        itr->clear();
                    }

                }
            }
        }
    }

    InterPtsType _interPts;

    for (const auto &pts : sortedPts) {
        _interPts.insert(_interPts.end(), pts.begin(), pts.end());
    }

    // probe, ob die schnittpunkte auf den kanten liegen

    if (!vtkPolyDataContactFilter::CheckInters(_interPts, pd)) {
        return false;
    }

    interPts.swap(_interPts);

    return true;

}

void vtkPolyDataContactFilter::InterPolys (vtkIdType idA, vtkIdType idB) {

#if (defined(_debA) && defined(_debB))
    _idA = idA; _idB = idB;

    if (_idA == _debA && _idB == _debB) {
        std::cout << "InterPolys() -> idA " << idA << ", idB " << idB << std::endl;
    }
#endif

    vtkIdType numA, numB;
    const vtkIdType *polyA, *polyB;

    newPdA->GetCellPoints(idA, numA, polyA);
    newPdB->GetCellPoints(idB, numB, polyB);

#if (defined(_debA) && defined(_debB))
    if (_idA == _debA && _idB == _debB) {
        double p[3];

        vtkIdType i;

        std::cout << _idA << " [";

        for (i = 0; i < numA; i++) {
            newPdA->GetPoint(polyA[i], p);

            std::cout << polyA[i] << ": [" << p[0] << ", " << p[1] << ", " << p[2] << "], ";
        }

        std::cout << "]" << std::endl;

        std::cout << _idB << " [";

        for (i = 0; i < numB; i++) {
            newPdB->GetPoint(polyB[i], p);

            std::cout << polyB[i] << ": [" << p[0] << ", " << p[1] << ", " << p[2] << "], ";
        }

        std::cout << "]" << std::endl;
    }

#endif

    // ebenen aufstellen

    double nA[3], nB[3], ptA[3], ptB[3], dA, dB;

    ComputeNormal(newPdA->GetPoints(), nA, numA, polyA);
    ComputeNormal(newPdB->GetPoints(), nB, numB, polyB);

    if (vtkMath::Dot(nA, nB) > .999999) {
        return;
    }

    newPdA->GetPoint(polyA[0], ptA);
    newPdB->GetPoint(polyB[0], ptB);

    dA = vtkMath::Dot(nA, ptA);
    dB = vtkMath::Dot(nB, ptB);

    if (!CheckNormal(newPdA->GetPoints(), numA, polyA, nA, dA)) {
        aborted = true;
        return;
    }

    if (!CheckNormal(newPdB->GetPoints(), numB, polyB, nB, dB)) {
        aborted = true;
        return;
    }

    // richtungsvektor

    double r[3];
    vtkMath::Cross(nA, nB, r);
    vtkMath::Normalize(r);

    // lsg. des lin. gls. mittels cramerscher regel

    int i = 0;

    for (int j = 1; j < 3; j++) {
        if (std::abs(r[j]) > std::abs(r[i])) {
            i = j;
        }
    }

    int inds[] = {1, 2};

    if (i == 1) {
        inds[0] = 0; inds[1] = 2;
    } else if (i == 2) {
        inds[0] = 0; inds[1] = 1;
    }

    double det = nA[inds[0]]*nB[inds[1]]-nB[inds[0]]*nA[inds[1]];

    if (std::abs(det) < 1e-10) {
        return;
    }

    // ein punkt auf der schnittgeraden der beiden ebenen

    double s[3];
    s[inds[0]] = (dA*nB[inds[1]]-dB*nA[inds[1]])/det;
    s[inds[1]] = (nA[inds[0]]*dB-nB[inds[0]]*dA)/det;
    s[i] = 0;

#if (defined(_debA) && defined(_debB))
    if (_idA == _debA && _idB == _debB) {
        std::cout << "r [" << r[0] << ", " << r[1] << ", " << r[2] << "]"
            << ", s [" << s[0] << ", " << s[1] << ", " << s[2] << "]"
            << std::endl;
    }
#endif

    InterPtsType intersA, intersB;

    if (!vtkPolyDataContactFilter::InterPolyLine(intersA, newPdA, numA, polyA, r, s, Src::A, nA)) {
        aborted = true;
        return;
    }

    if (!vtkPolyDataContactFilter::InterPolyLine(intersB, newPdB, numB, polyB, r, s, Src::B, nB)) {
        aborted = true;
        return;
    }

#if (defined(_debA) && defined(_debB))
    if (_idA == _debA && _idB == _debB) {

        std::cout << "intersA" << std::endl;

        for (auto inter : intersA) {
            std::cout << inter << std::endl;
        }

        std::cout << "intersB" << std::endl;

        for (auto inter : intersB) {
            std::cout << inter << std::endl;
        }

    }

#endif

    if (intersA.size() != 0 && intersB.size() != 0
        && intersA.size()%2 == 0 && intersB.size()%2 == 0) {

        AddContactLines(intersA, intersB, idA, idB);
    }

}

void vtkPolyDataContactFilter::OverlapLines (OverlapsType &ols, InterPtsType &intersA, InterPtsType &intersB, vtkIdType vtkNotUsed(idA), vtkIdType vtkNotUsed(idB)) {

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
                    ols.push_back(Add(*itr2, *(itr2+1), *itr, *(itr+1)));
                } else {
                    ols.push_back(Add(*itr2, *(itr+1), *itr, *(itr2+1)));
                }
            } else if (itr2->t <= itr->t && (itr2+1)->t > itr->t) {
                if ((itr+1)->t < (itr2+1)->t) {
                    ols.push_back(Add(*itr, *(itr+1), *itr2, *(itr2+1)));
                } else {
                    ols.push_back(Add(*itr, *(itr2+1), *itr2, *(itr+1)));
                }
            }
        }
    }

}

void vtkPolyDataContactFilter::AddContactLines (InterPtsType &intersA, InterPtsType &intersB, vtkIdType idA, vtkIdType idB) {

    OverlapsType overlaps;
    OverlapLines(overlaps, intersA, intersB, idA, idB);

    OverlapsType::const_iterator itr;

    for (itr = overlaps.begin(); itr != overlaps.end(); ++itr) {
        auto &f = std::get<0>(*itr),
            &s = std::get<1>(*itr);

        if (f.src == Src::A) {
            if (edgesA.count(f.edge) == 1) {
                invalidA = true;
            }
        }

        if (s.src == Src::A) {
            if (edgesA.count(s.edge) == 1) {
                invalidA = true;
            }
        }

        if (f.src == Src::B) {
            if (edgesB.count(f.edge) == 1) {
                invalidB = true;
            }
        }

        if (s.src == Src::B) {
            if (edgesB.count(s.edge) == 1) {
                invalidB = true;
            }
        }

        vtkIdList *linePts = vtkIdList::New();

        linePts->InsertNextId(contLines->GetPoints()->InsertNextPoint(f.pt));
        linePts->InsertNextId(contLines->GetPoints()->InsertNextPoint(s.pt));

        contLines->InsertNextCell(VTK_LINE, linePts);

        const vtkIdType tupleA[] = {f.srcA, s.srcA};
        const vtkIdType tupleB[] = {f.srcB, s.srcB};

        sourcesA->InsertNextTypedTuple(tupleA);
        sourcesB->InsertNextTypedTuple(tupleB);

        linePts->Delete();

        contA->InsertNextValue(idA);
        contB->InsertNextValue(idB);
    }

}

int vtkPolyDataContactFilter::InterOBBNodes (vtkOBBNode *nodeA, vtkOBBNode *nodeB, vtkMatrix4x4 *vtkNotUsed(mat), void *caller) {
    vtkPolyDataContactFilter *self = reinterpret_cast<vtkPolyDataContactFilter*>(caller);

    vtkIdList *cellsA = nodeA->Cells;
    vtkIdList *cellsB = nodeB->Cells;

    vtkIdType numCellsA = cellsA->GetNumberOfIds(),
        numCellsB = cellsB->GetNumberOfIds();

    vtkIdType i, j, ci, cj;

    for (i = 0; i < numCellsA && !self->aborted; i++) {
        ci = cellsA->GetId(i);

        for (j = 0; j < numCellsB && !self->aborted; j++) {
            cj = cellsB->GetId(j);

            self->InterPolys(ci, cj);
        }
    }

    return 0;
}

bool vtkPolyDataContactFilter::CheckInters (InterPtsType &interPts, vtkPolyData *pd) {
    double ptA[3],
        ptB[3],
        v[3],
        w[3],
        k,
        l,
        alpha,
        d;

    for (auto &p : interPts) {
        pd->GetPoint(p.edge.f, ptA);
        pd->GetPoint(p.edge.g, ptB);

        vtkMath::Subtract(ptA, ptB, v);
        vtkMath::Normalize(v);
        vtkMath::Subtract(ptA, p.pt, w);

        k = vtkMath::Norm(w);
        l = vtkMath::Dot(v, w);
        alpha = std::acos(l/k);

        if (std::isnan(alpha)) {
            continue;
        }

        d = std::sin(alpha)*k;

        if (d < 1e-5) {
            continue;
        }

        return false;

    }

    return true;

}
