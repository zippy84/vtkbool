/*
Copyright 2012-2023 Ronald Römer

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
#include <vtkPointData.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleStrip.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkFeatureEdges.h>
#include <vtkCellIterator.h>
#include <vtkCellArrayIterator.h>

#include <vtkCellArray.h>

#include "vtkPolyDataContactFilter.h"
#include "Utilities.h"

#undef DEBUG

vtkStandardNewMacro(vtkPolyDataContactFilter);

vtkPolyDataContactFilter::vtkPolyDataContactFilter () {

    contLines = vtkPolyData::New();
    contLines->Allocate(1000);

    contPts = vtkPoints::New();
    contPts->SetDataTypeToDouble();
    contLines->SetPoints(contPts);

    contA = vtkIdTypeArray::New();
    contB = vtkIdTypeArray::New();

    contA->SetName("cA");
    contB->SetName("cB");

    sourcesA = vtkIdTypeArray::New();
    sourcesA->SetNumberOfComponents(2);

    sourcesB = vtkIdTypeArray::New();
    sourcesB->SetNumberOfComponents(2);

    sourcesA->SetName("sourcesA");
    sourcesB->SetName("sourcesB");

    SetNumberOfInputPorts(2);
    SetNumberOfOutputPorts(3);

    invalidA = false;
    invalidB = false;
}

vtkPolyDataContactFilter::~vtkPolyDataContactFilter () {
    sourcesB->Delete();
    sourcesA->Delete();

    contB->Delete();
    contA->Delete();

    contPts->Delete();
    contLines->Delete();
}

int vtkPolyDataContactFilter::ProcessRequest (vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

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

        // durchführung der aufgabe

        pdA = vtkPolyData::New();
        pdA->DeepCopy(_pdA);

        pdB = vtkPolyData::New();
        pdB->DeepCopy(_pdB);

        PreparePolyData(pdA);
        PreparePolyData(pdB);

        if (pdA->GetNumberOfCells() == 0) {
            vtkErrorMacro("First input does not contain any supported cells.");
            return 1;
        }

        if (pdB->GetNumberOfCells() == 0) {
            vtkErrorMacro("Second input does not contain any supported cells.");
            return 1;
        }

        GetInvalidEdges(pdA, edgesA);
        GetInvalidEdges(pdB, edgesB);

        // anlegen der obb-trees

        vtkOBBTree *obbA = vtkOBBTree::New();
        obbA->SetDataSet(pdA);
        obbA->SetNumberOfCellsPerNode(1);
        obbA->BuildLocator();

        vtkOBBTree *obbB = vtkOBBTree::New();
        obbB->SetDataSet(pdB);
        obbB->SetNumberOfCellsPerNode(1);
        obbB->BuildLocator();

        vtkMatrix4x4 *mat = vtkMatrix4x4::New();

        obbA->IntersectWithOBBTree(obbB, mat, InterOBBNodes, this);

        if (invalidA) {
            vtkErrorMacro("First input has non-manifold edges.");
            return 1;
        }

        if (invalidB) {
            vtkErrorMacro("Second input has non-manifold edges.");
            return 1;
        }

        contLines->GetCellData()->AddArray(contA);
        contLines->GetCellData()->AddArray(contB);

        contLines->GetCellData()->AddArray(sourcesA);
        contLines->GetCellData()->AddArray(sourcesB);

        contLines->RemoveDeletedCells();

        vtkCleanPolyData *clean = vtkCleanPolyData::New();
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

        clean->Delete();
        mat->Delete();
        obbB->Delete();
        obbA->Delete();

        resultB->DeepCopy(pdA);
        resultC->DeepCopy(pdB);

        pdB->Delete();
        pdA->Delete();

    }

    return 1;

}

void vtkPolyDataContactFilter::PreparePolyData (vtkPolyData *pd) {

    pd->GetCellData()->Initialize();
    pd->GetPointData()->Initialize();

    vtkIdTypeArray *cellIds = vtkIdTypeArray::New();

    vtkCellIterator *cellItr = pd->NewCellIterator();
    for (cellItr->InitTraversal(); !cellItr->IsDoneWithTraversal(); cellItr->GoToNextCell()) {
        cellIds->InsertNextValue(cellItr->GetCellId());
    }

    vtkIdType cellId;

    for (cellItr->InitTraversal(); !cellItr->IsDoneWithTraversal(); cellItr->GoToNextCell()) {
        cellId = cellItr->GetCellId();

        if (cellItr->GetCellType() == VTK_QUAD) {
            vtkIdList *ptIds = cellItr->GetPointIds();

            vtkPoints *pts = cellItr->GetPoints();

            double n[3];

            ComputeNormal(pd->GetPoints(), n, 4, ptIds->GetPointer(0));

            double dA = vtkMath::Dot(n, pts->GetPoint(0)),
                dB = vtkMath::Dot(n, pts->GetPoint(1))-dA;

            if (std::abs(dB) > 1e-6) {
                // nur wenn nicht auf einer ebene

                dA = vtkMath::Distance2BetweenPoints(pts->GetPoint(0), pts->GetPoint(2));
                dB = vtkMath::Distance2BetweenPoints(pts->GetPoint(1), pts->GetPoint(3));

                auto newCellA = vtkSmartPointer<vtkIdList>::New();
                auto newCellB = vtkSmartPointer<vtkIdList>::New();

                newCellA->SetNumberOfIds(3);
                newCellB->SetNumberOfIds(3);

                if (dA < dB) {
                    newCellA->SetId(0, ptIds->GetId(1));
                    newCellA->SetId(1, ptIds->GetId(2));
                    newCellA->SetId(2, ptIds->GetId(3));

                    newCellB->SetId(0, ptIds->GetId(3));
                    newCellB->SetId(1, ptIds->GetId(0));
                    newCellB->SetId(2, ptIds->GetId(1));
                } else {
                    newCellA->SetId(0, ptIds->GetId(0));
                    newCellA->SetId(1, ptIds->GetId(1));
                    newCellA->SetId(2, ptIds->GetId(2));

                    newCellB->SetId(0, ptIds->GetId(2));
                    newCellB->SetId(1, ptIds->GetId(3));
                    newCellB->SetId(2, ptIds->GetId(0));
                }

                pd->InsertNextCell(VTK_TRIANGLE, newCellA);
                cellIds->InsertNextValue(cellId);

                pd->InsertNextCell(VTK_TRIANGLE, newCellB);
                cellIds->InsertNextValue(cellId);

                pd->DeleteCell(cellId);
            }

        } else if (cellItr->GetCellType() == VTK_TRIANGLE_STRIP) {
            vtkIdList *ptIds = cellItr->GetPointIds();

            vtkCellArray *cells = vtkCellArray::New();

            vtkTriangleStrip::DecomposeStrip(cellItr->GetNumberOfPoints(), ptIds->GetPointer(0), cells);

            vtkIdType n;
            const vtkIdType *pts;

            for (cells->InitTraversal(); cells->GetNextCell(n, pts);) {
                if (pts[0] != pts[1] && pts[1] != pts[2] && pts[2] != pts[0]) {
                    pd->InsertNextCell(VTK_TRIANGLE, n, pts);
                    cellIds->InsertNextValue(cellId);
                }

            }

            cells->Delete();

            pd->DeleteCell(cellId);

        } else if (cellItr->GetCellType() != VTK_TRIANGLE && cellItr->GetCellType() != VTK_POLYGON) {
            pd->DeleteCell(cellId);
        }
    }

    cellItr->Delete();

    cellIds->SetName("OrigCellIds");
    pd->GetCellData()->SetScalars(cellIds);

    cellIds->Delete();

    pd->RemoveDeletedCells();
    pd->BuildLinks();

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

    vtkIdType i, j, indexA, indexB, indexA_, indexB_;

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

        j = 0;

        for (i = 0; i < cellsA->GetNumberOfIds(); i++) {
            pd->GetCellPoints(cellsA->GetId(i), cell);

            indexA = cell->IsId(idA);
            indexB = cell->IsId(idB);

            indexA_ = indexA+1;

            if (indexA_ == cell->GetNumberOfIds()) {
                indexA_ = 0;
            }

            indexB_ = indexB+1;

            if (indexB_ == cell->GetNumberOfIds()) {
                indexB_ = 0;
            }

            if (cell->GetId(indexA_) == idB || cell->GetId(indexB_) == idA) {
                j++;
            }
        }

        if (j > 2) {
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

                End end = End::NONE;

                if (s > -1e-6 && s < 1e-6) {
                    end = End::BEGIN;
                } else if (s > l-1e-6 && s < l+1e-6) {
                    end = End::END;
                }

                interPts.emplace_back(ptA[0]+t*r[0], ptA[1]+t*r[1], ptA[2]+t*r[2], t, idA, idB, end, src);

            }

        } else {
            // parallel

            double vA[3], vB[3], cA[3], cB[3], dA, dB;

            vtkMath::Subtract(eA, ptA, vA);
            vtkMath::Subtract(eA, ptB, vB);
            vtkMath::Cross(vA, vB, cA);

            double dotA = vtkMath::Dot(vA, r);

            vtkMath::Subtract(eB, ptA, vA);
            vtkMath::Subtract(eB, ptB, vB);
            vtkMath::Cross(vA, vB, cB);

            double dotB = vtkMath::Dot(vA, r);

            dA = vtkMath::Norm(cA);
            dB = vtkMath::Norm(cB);

            if (dA < 1e-4 || dB < 1e-4) {
#ifdef DEBUG
                std::cout << "congruent lines with d (" << dA << ", " << dB << ") and t (" << dotA << ", " << dotB << ") and l " << l << std::endl;
#endif

                interPts.emplace_back(ptA[0]+dotA*r[0], ptA[1]+dotA*r[1], ptA[2]+dotA*r[2], dotA, idA, idB, End::BEGIN, src);
                interPts.emplace_back(ptA[0]+dotB*r[0], ptA[1]+dotB*r[1], ptA[2]+dotB*r[2], dotB, idA, idB, End::END, src);

            }

        }

    } else {
        // windschief
    }

}

void vtkPolyDataContactFilter::InterPolyLine (InterPtsType &interPts, vtkPolyData *pd, vtkIdType num, const vtkIdType *poly, const double *r, const double *pt, Src src, const double *n) {

#ifdef DEBUG
    std::cout << "InterPolyLine()" << std::endl;
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
        return;
    }

    struct Cmp {
        bool operator() (const double &l, const double &r) const {
            long a = std::lround(l*1e5),
                b = std::lround(r*1e5);

            return a < b;
        }
    };

    std::map<double, InterPtsType, Cmp> paired;

    for (auto &p : interPts) {
        paired[p.t].push_back(p);
    }

    std::vector<InterPtsType> _paired;

    for (auto &p : paired) {
        InterPtsType &pts = p.second;

        if (pts.size() == 1 && pts.front().end != End::NONE) {
            // hier fehlt der zweite punkt
            pts.push_back(pts.back());
        }

        _paired.push_back(pts);
    }

    // trivial

    if (_paired.front().size() == 2) {
        _paired.front().pop_back();
    }

    if (_paired.back().size() == 2) {
        _paired.back().pop_back();
    }

    // ...

    std::map<vtkIdType, double> ends;

    for (const auto &pts : _paired) {
        auto &last = pts.back();

        if (last.end != End::NONE) {
            ends.emplace(last.end == End::BEGIN ? last.edge.f : last.edge.g, last.t);
        }
    }

    double s[3], d;

    vtkMath::Cross(n, r, s);
    d = vtkMath::Dot(s, pt);

    double ptA[3], ptB[3], dA, dB, vA[3], vB[3], tA, tB;

    vtkIdType id, prev, next;

    for (auto &pts : _paired) {
        auto &last = pts.back();

        if (last.end != End::NONE) {
            if (last.end == End::BEGIN) {
                id = last.edge.f;
                next = last.edge.g;
                prev = std::find_if(edges.begin(), edges.end(), [&id](const Pair &edge) { return edge.g == id; })->f;
            } else {
                id = last.edge.g;
                prev = last.edge.f;
                next = std::find_if(edges.begin(), edges.end(), [&id](const Pair &edge) { return edge.f == id; })->g;
            }

            if (pts.size() == 2) {
                if (ends.count(next) == 0 && ends.count(prev) == 1) {
                    pd->GetPoint(next, ptA);
                    dA = vtkMath::Dot(s, ptA)-d;

                    if ((last.t > ends.at(prev) && dA > 0) || (last.t < ends.at(prev) && dA < 0)) {
                        // tasche
                        pts.pop_back();
                    }

                    continue;

                } else if (ends.count(next) == 1 && ends.count(prev) == 0) {
                    pd->GetPoint(prev, ptB);
                    dB = vtkMath::Dot(s, ptB)-d;

                    if ((last.t > ends.at(next) && dB < 0) || (last.t < ends.at(next) && dB > 0)) {
                        // tasche
                        pts.pop_back();
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
                    if (pts.size() == 2) {
                        pts.pop_back();
                    }

                } else {
                    vtkMath::Subtract(ptA, pt, vA);
                    vtkMath::Subtract(ptB, pt, vB);

                    tA = vtkMath::Dot(vA, r);
                    tB = vtkMath::Dot(vB, r);

                    if ((tB > tA) == std::signbit(dA)) {
                        pts.clear();
                    }

                }
            }
        }
    }

    // ...

    InterPtsType _interPts;

    for (const auto &pts : _paired) {
        _interPts.insert(_interPts.end(), pts.begin(), pts.end());
    }

    interPts.swap(_interPts);

}

void vtkPolyDataContactFilter::InterPolys (vtkIdType idA, vtkIdType idB) {

#ifdef DEBUG
    std::cout << "InterPolys() -> idA " << idA << ", idB " << idB << std::endl;
#endif

    vtkIdType numA, numB;
    const vtkIdType *polyA, *polyB;

    pdA->GetCellPoints(idA, numA, polyA);
    pdB->GetCellPoints(idB, numB, polyB);

    // ebenen aufstellen

    double nA[3], nB[3], ptA[3], ptB[3], dA, dB;

    ComputeNormal(pdA->GetPoints(), nA, numA, polyA);
    ComputeNormal(pdB->GetPoints(), nB, numB, polyB);

    pdA->GetPoint(polyA[0], ptA);
    pdB->GetPoint(polyB[0], ptB);

    dA = vtkMath::Dot(nA, ptA);
    dB = vtkMath::Dot(nB, ptB);

    // richtungsvektor

    double r[3];
    vtkMath::Cross(nA, nB, r);
    vtkMath::Normalize(r);

    // std::cout << r[0] << ", "
    //     << r[1] << ", "
    //     << r[2] << std::endl;

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

#ifdef DEBUG
    std::cout << "r [" << r[0] << ", " << r[1] << ", " << r[2] << "]"
        << ", s [" << s[0] << ", " << s[1] << ", " << s[2] << "]"
        << std::endl;
#endif

    InterPtsType intersA, intersB;

    vtkPolyDataContactFilter::InterPolyLine(intersA, pdA, numA, polyA, r, s, Src::A, nA);
    vtkPolyDataContactFilter::InterPolyLine(intersB, pdB, numB, polyB, r, s, Src::B, nB);

    // probe, ob die schnittpunkte auf den kanten liegen
    // bei ungenauen normalen ist das manchmal nicht der fall

    vtkPolyDataContactFilter::CheckInters(intersA, pdA, idA, idB);
    vtkPolyDataContactFilter::CheckInters(intersB, pdB, idA, idB);

#ifdef DEBUG
    std::cout << "intersA" << std::endl;

    for (auto inter : intersA) {
        std::cout << inter << std::endl;
    }

    std::cout << "intersB" << std::endl;

    for (auto inter : intersB) {
        std::cout << "B " << inter << std::endl;
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

#ifdef DEBUG
        std::cout << "f " << f << std::endl;
        std::cout << "s " << s << std::endl;
#endif

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

        linePts->InsertNextId(contPts->InsertNextPoint(f.pt));
        linePts->InsertNextId(contPts->InsertNextPoint(s.pt));

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

    for (i = 0; i < numCellsA; i++) {
        ci = cellsA->GetId(i);

        for (j = 0; j < numCellsB; j++) {
            cj = cellsB->GetId(j);

            self->InterPolys(ci, cj);
        }
    }

    return 0;
}

void vtkPolyDataContactFilter::CheckInters (InterPtsType &interPts, vtkPolyData *pd, vtkIdType vtkNotUsed(idA), vtkIdType vtkNotUsed(idB)) {
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

        // if (p.src == Src::A) {
        //     std::cout << "? A ";
        // } else {
        //     std::cout << "? B ";
        // }

        // std::cout << idA << ", " << idB << ": "
        //     << std::fixed
        //     << d
        //     << ", [" << p.pt[0] << ", " << p.pt[1] << ", " << p.pt[2] << "], "
        //     << p.edge
        //     << std::endl;

    }

}
