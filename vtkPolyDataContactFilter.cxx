/*
Copyright 2012-2020 Ronald Römer

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
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleStrip.h>
#include <vtkDoubleArray.h>

#include <vtkCellArray.h>

#include "vtkPolyDataContactFilter.h"
#include "Utilities.h"

// #undef DEBUG

vtkStandardNewMacro(vtkPolyDataContactFilter);

vtkPolyDataContactFilter::vtkPolyDataContactFilter () {

    contLines = vtkPolyData::New();
    contLines->Allocate(1000);

    contPts = vtkPoints::New();
    contPts->SetDataTypeToDouble();
    contLines->SetPoints(contPts);

    contA = vtkIntArray::New();
    contB = vtkIntArray::New();

    contA->SetName("cA");
    contB->SetName("cB");

    sourcesA = vtkIntArray::New();
    sourcesA->SetNumberOfComponents(2);

    sourcesB = vtkIntArray::New();
    sourcesB->SetNumberOfComponents(2);

    sourcesA->SetName("sourcesA");
    sourcesB->SetName("sourcesB");

    SetNumberOfInputPorts(2);
    SetNumberOfOutputPorts(3);

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

        if (pdA->GetNumberOfCells() == 0 || pdB->GetNumberOfCells() == 0) {
            vtkErrorMacro("One of the inputs does not contain any supported cells.");

            return 1;
        }

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

        contLines->GetCellData()->AddArray(contA);
        contLines->GetCellData()->AddArray(contB);

        contLines->GetCellData()->AddArray(sourcesA);
        contLines->GetCellData()->AddArray(sourcesB);

        vtkCleanPolyData *clean = vtkCleanPolyData::New();
        clean->SetInputData(contLines);
        clean->ToleranceIsAbsoluteOn();
        clean->SetAbsoluteTolerance(1e-5);
        clean->Update();

        resultA->DeepCopy(clean->GetOutput());

        vtkIdType numCellsA = resultA->GetNumberOfCells();

        for (vtkIdType i = 0; i < numCellsA; i++) {
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

    vtkIdType numCells = pd->GetNumberOfCells();

    vtkIntArray *cellIds = vtkIntArray::New();
    vtkIntArray *stripIds = vtkIntArray::New();

    for (vtkIdType i = 0; i < numCells; i++) {
        cellIds->InsertNextValue(i);

        if (pd->GetCellType(i) == VTK_TRIANGLE_STRIP) {
            stripIds->InsertNextValue(i);
        }
    }

    vtkCellArray *cells = vtkCellArray::New();

    vtkCellArray *strips = pd->GetStrips();

    vtkIdType n;
    vtkIdType *pts;

    vtkIdType i = 0;

    for (strips->InitTraversal(); strips->GetNextCell(n, pts);) {
        cells->Reset();

        vtkTriangleStrip::DecomposeStrip(n, pts, cells);

        for (cells->InitTraversal(); cells->GetNextCell(n, pts);) {
            if (pts[0] != pts[1] && pts[1] != pts[2] && pts[2] != pts[0]) {
                pd->InsertNextCell(VTK_TRIANGLE, n, pts);
                cellIds->InsertNextValue(stripIds->GetValue(i));

                numCells++;
            }

        }

        i++;

    }

    int type;
    for (vtkIdType i = 0; i < numCells; i++) {
        type = pd->GetCellType(i);

        if (type != VTK_POLYGON && type != VTK_QUAD && type != VTK_TRIANGLE) {
            pd->DeleteCell(i);
        }

    }

    cellIds->SetName("OrigCellIds");

    pd->GetCellData()->SetScalars(cellIds);

    cells->Delete();
    stripIds->Delete();
    cellIds->Delete();

    pd->RemoveDeletedCells();

}

void vtkPolyDataContactFilter::InterEdgeLine (InterPtsType &interPts, const double *eA, const double *eB, const double *r, const double *ptA) {

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

                vtkIdType end {NO_USE};

                if (s > -1e-6 && s < 1e-6) {
                    end = 0;
                } else if (s > l-1e-6 && s < l+1e-6) {
                    end = 1;
                }

                interPts.emplace_back(t, end, ptA[0]+t*r[0], ptA[1]+t*r[1], ptA[2]+t*r[2]);

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

                interPts.emplace_back(dotA, 0, ptA[0]+dotA*r[0], ptA[1]+dotA*r[1], ptA[2]+dotA*r[2]);
                interPts.emplace_back(dotB, 1, ptA[0]+dotB*r[0], ptA[1]+dotB*r[1], ptA[2]+dotB*r[2]);

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

    InterPtsType interPtsA;

    // durchläuft die kanten und ermittelt die schnittpunkte

    double ptA[3],
        ptB[3],
        dA,
        dB;

    vtkIdType i, j;

    for (i = 0; i < num; i++) {
        j = i == num-1 ? 0 : i+1;

#ifdef DEBUG
        std::cout << "(" << poly[i] << ", " << poly[j] << ")" << std::endl;
#endif

        // kante ermitteln
        pd->GetPoint(poly[i], ptA);
        pd->GetPoint(poly[j], ptB);

        // schnittpunkt

        InterPtsType interPtsB;
        vtkPolyDataContactFilter::InterEdgeLine(interPtsB, ptA, ptB, r, pt);

        for (InterPt &p : interPtsB) {
            p.src = src;

            if (p.onEdge) {
                p.edge[0] = i;
                p.edge[1] = j;

                if (p.end != NO_USE) {
                    p.end = p.end == 0 ? i : j;
                }

#ifdef DEBUG
                std::cout << "-> " << p << std::endl;
#endif

                interPtsA.push_back(p);
            }
        }
    }

    if (!interPtsA.empty()) {
        struct Cmp {
            bool operator() (const double &l, const double &r) const {
                long a = std::lround(l*1e5),
                    b = std::lround(r*1e5);

                return a < b;
            }
        };

        std::map<double, InterPtsType, Cmp> paired;

        for (auto &p : interPtsA) {
            paired[p.t].push_back(p);

#ifdef DEBUG
            std::cout << p << std::endl;
#endif
        }

        std::vector<InterPtsType> _paired;

        for (auto &p : paired) {
            InterPtsType &pts = p.second;

            if (pts.size() == 1 && pts.front().end != NO_USE) {
                // hier fehlt der zweite punkt
                pts.push_back(pts.back());
            }

            _paired.push_back(pts);
        }

        // trivial

        auto &pairA = _paired.front(),
            &pairB = _paired.back();

        if (pairA.size() == 2) {
            pairA.pop_back();
        }

        if (pairB.size() == 2) {
            pairB.pop_back();
        }

        double m[3], q[3], d, e, t;

        vtkMath::Cross(n, r, m);
        d = vtkMath::Dot(m, pt);

        std::map<vtkIdType, double> ends;

        for (const auto &p : _paired) {
            if (p.back().end != NO_USE) {
                ends.emplace(p.back().end, p.back().t);
            }
        }

        double vA[3],
            vB[3],
            tA,
            tB;

        vtkIdType before, after;

        for (auto &p : _paired) {
            const InterPt &dupl = p.back();

            if (dupl.end != NO_USE) {
                before = dupl.end == 0 ? num-1 : dupl.end-1;
                after = dupl.end == num-1 ? 0 : dupl.end+1;

                if (p.size() == 2) {
                    if (ends.count(after) == 0 && ends.count(before) == 1) {
                        pd->GetPoint(poly[after], q);
                        e = vtkMath::Dot(m, q)-d;

                        t = ends[before];

                        if ((dupl.t > t && e > 0) || (dupl.t < t && e < 0)) {
                            // tasche
                            p.pop_back();
                        }

                        continue;

                    } else if (ends.count(before) == 0 && ends.count(after) == 1) {
                        pd->GetPoint(poly[before], q);
                        e = vtkMath::Dot(m, q)-d;

                        t = ends[after];

                        if ((dupl.t > t && e < 0) || (dupl.t < t && e > 0)) {
                            // tasche
                            p.pop_back();
                        }

                        continue;

                    }
                }

                if (ends.count(before) == 0 && ends.count(after) == 0) {
                    pd->GetPoint(poly[after], ptA);
                    pd->GetPoint(poly[before], ptB);

                    dA = vtkMath::Dot(m, ptA)-d;
                    dB = vtkMath::Dot(m, ptB)-d;

                    if (std::signbit(dA) != std::signbit(dB)) {
                        if (p.size() == 2) {
                            p.pop_back();
                        }

                    } else {
                        vtkMath::Subtract(ptA, pt, vA);
                        vtkMath::Subtract(ptB, pt, vB);

                        tA = vtkMath::Dot(vA, r);
                        tB = vtkMath::Dot(vB, r);

                        if ((tB > tA) == std::signbit(dA)) {
                            p.clear();
                        }

                    }
                }
            }
        }

        for (const auto &pair : _paired) {
            for (const InterPt &p : pair) {
#ifdef DEBUG
                std::cout << "* " << p << std::endl;
#endif

                interPts.push_back(p);
            }
        }

    }

}

void vtkPolyDataContactFilter::InterPolys (vtkIdType idA, vtkIdType idB) {

#ifdef DEBUG
    std::cout << "InterPolys() -> idA " << idA << ", idB " << idB << std::endl;
#endif

    vtkIdType numA, numB;
    vtkIdType *polyA, *polyB;

    pdA->GetCellPoints(idA, numA, polyA);
    pdB->GetCellPoints(idB, numB, polyB);

    // ebenen aufstellen

    double nA[3], nB[3], ptA[3], ptB[3], dA, dB;

    ComputeNormal2(pdA, nA, numA, polyA);
    ComputeNormal2(pdB, nB, numB, polyB);

    pdA->GetPoint(polyA[0], ptA);
    pdB->GetPoint(polyB[0], ptB);

    dA = vtkMath::Dot(nA, ptA);
    dB = vtkMath::Dot(nB, ptB);

    // sind die ebenen parallel?

    double p = std::abs(vtkMath::Dot(nA, nB));

    if (p < 0.999999) {

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

#ifdef DEBUG
        std::cout << "intersA " << intersA.size()
            << ", intersB " << intersB.size()
            << std::endl;
#endif


        if (intersA.size() != 0 && intersB.size() != 0
            && intersA.size()%2 == 0 && intersB.size()%2 == 0) {

            OverlapsType overlaps;
            vtkPolyDataContactFilter::OverlapLines(overlaps, intersA, intersB);

            OverlapsType::const_iterator itr;

            for (itr = overlaps.begin(); itr != overlaps.end(); ++itr) {

                auto &f = itr->first,
                    &s = itr->second;

#ifdef DEBUG
                std::cout << "f " << f << std::endl;
                std::cout << "s " << s << std::endl;
#endif

                vtkIdList *linePts = vtkIdList::New();

                linePts->InsertNextId(contPts->InsertNextPoint(f.pt));
                linePts->InsertNextId(contPts->InsertNextPoint(s.pt));

                contLines->InsertNextCell(VTK_LINE, linePts);

                sourcesA->InsertNextTuple2(f.srcA, s.srcA);
                sourcesB->InsertNextTuple2(f.srcB, s.srcB);

                linePts->Delete();

                contA->InsertNextValue(idA);
                contB->InsertNextValue(idB);

            }

        }

    }

}

void vtkPolyDataContactFilter::OverlapLines (OverlapsType &ols, InterPtsType &intersA, InterPtsType &intersB) {

    auto Add = [](InterPt &a, InterPt &b, InterPt &c, InterPt &d) {
        auto p = std::make_pair(a, b);

        p.first.Merge(c);
        p.second.Merge(d);

        return p;
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
