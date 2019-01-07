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

#undef DEBUG

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

        /*
        for (int i = 0; i < pdA->GetNumberOfCells(); i++) {
            for (int j = 0; j < pdB->GetNumberOfCells(); j++) {
                IntersectPolys(i, j);
            }
        }
        */

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

        int numCellsA = resultA->GetNumberOfCells();

        for (int i = 0; i < numCellsA; i++) {
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

    for (int i = 0; i < numCells; i++) {
        cellIds->InsertNextValue(i);

        if (pd->GetCellType(i) == VTK_TRIANGLE_STRIP) {
            stripIds->InsertNextValue(i);
        }
    }

    vtkCellArray *cells = vtkCellArray::New();

    vtkCellArray *strips = pd->GetStrips();

    vtkIdType n;
    vtkIdType *pts;

    int i = 0;

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
    for (int i = 0; i < numCells; i++) {
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

InterPtType vtkPolyDataContactFilter::InterEdgeLine (double *eA, double *eB, double *r, double *pt, int _pid) {

    InterPtType inter;

    // richtungsvektor der kante bestimmen

    double e[3];
    vtkMath::Subtract(eB, eA, e);
    double l = vtkMath::Normalize(e);

    double p[3];
    vtkMath::Subtract(eA, pt, p);

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

                inter.pt[0] = pt[0]+t*r[0];
                inter.pt[1] = pt[1]+t*r[1];
                inter.pt[2] = pt[2]+t*r[2];

                inter.onEdge = true;

                inter.t = t;

                if (s > -1e-6 && s < 1e-6) {
                    inter.end = 0;
                } else if (s > l-1e-6 && s < l+1e-6) {
                    inter.end = 1;
                }

            }

        } else {
            // parallel
        }

    } else {
        // windschief
    }

    return inter;

}

InterPtsType vtkPolyDataContactFilter::InterPolyLine (vtkPoints *pts, vtkIdList *poly, double *r, double *pt, _Src src, int _pid) {

#ifdef DEBUG
    std::cout << "InterPolyLine()" << std::endl;
#endif

    InterPtsType interPts;

    // durchläuft die kanten und ermittelt die schnittpunkte

    int numPts = poly->GetNumberOfIds();

    for (int i = 0; i < numPts; i++) {
        int j = (i+1)%numPts;

        // kante ermitteln

        double ptA[3], ptB[3];

        pts->GetPoint(poly->GetId(i), ptA);
        pts->GetPoint(poly->GetId(j), ptB);

        // schnittpunkt

        InterPtType inter(InterEdgeLine(ptA, ptB, r, pt, _pid));

        inter.src = src;

        if (inter.onEdge) {
            inter.edge[0] = i;
            inter.edge[1] = j;

            if (inter.end != NO_USE) {
                inter.end = inter.end == 0 ? i : j;
            }

            // if (_pid == 34026) {
            //     std::cout << inter << std::endl;
            // }

            interPts.push_back(inter);
        }

    }

    std::sort(interPts.begin(), interPts.end());

    // jetzt müssen noch spezielle gelöscht werden

    InterPtsType interPts2, interPts3;

    InterPtsType::const_iterator itr;

#ifdef DEBUG
    for (itr = interPts.begin(); itr != interPts.end(); ++itr) {
        std::cout << *itr << std::endl;
    }
#endif

    if (interPts.size() > 1) {

        double n[3];
        ComputeNormal(pts, n, poly);

        double p[3], d;

        vtkMath::Cross(r, n, p);
        vtkMath::Normalize(p);

        d = vtkMath::Dot(p, pt);

        std::vector<const InterPtType*> omits;

        for (itr = interPts.begin(); itr != interPts.end()-1; ++itr) {

            int shared = -1;

            if ((itr->edge[0] == (itr+1)->edge[0])
                || (itr->edge[0] == (itr+1)->edge[1])) {
                shared = itr->edge[0];
            } else if ((itr->edge[1] == (itr+1)->edge[0])
                || (itr->edge[1] == (itr+1)->edge[1])) {
                shared = itr->edge[1];
            }

            if (shared > -1 && (itr->end == shared || (itr+1)->end == shared)) {

                int prev = (shared+numPts-1)%numPts;
                int next = (shared+1)%numPts;

                double prevPt[3], nextPt[3];

                pts->GetPoint(poly->GetId(prev), prevPt);
                pts->GetPoint(poly->GetId(next), nextPt);

                double prevD = vtkMath::Dot(p, prevPt)-d;
                double nextD = vtkMath::Dot(p, nextPt)-d;

                if ((prevD < -1e-6 && nextD < -1e-6)
                    || (prevD > 1e-6 && nextD > 1e-6)) {
                    omits.push_back(&*itr);
                    omits.push_back(&*(itr+1));
                } else {
                    if (itr->end == shared) {
                        omits.push_back(&*(itr+1));
                    } else {
                        omits.push_back(&*itr);
                    }
                }

            }

        }

#ifdef DEBUG
        std::vector<const InterPtType*>::const_iterator itr2;

        std::cout << "omits [";
        for (itr2 = omits.begin(); itr2 != omits.end(); ++itr2) {
            std::cout << (*itr2)->ind << ", ";
        }
        std::cout << "]" << std::endl;
#endif

        std::set<int> ends;

        for (itr = interPts.begin(); itr != interPts.end(); ++itr) {
            if (std::count(omits.begin(), omits.end(), &*itr) == 0) {
                interPts2.push_back(*itr);

                ends.insert(itr->end);
            }
        }

        // löst probleme mit kongruenten kanten

        // punkte die durch jeweils zwei kongruente kanten begrenzt sind

        InterPtsType::iterator itr3;

        for (int i = 0; i < numPts; i++) {
            if (ends.count(i) == 0) {
                double endPt[3], endD;

                pts->GetPoint(poly->GetId(i), endPt);

                endD = std::abs(vtkMath::Dot(p, endPt)-d);

                if (endD < 1e-6) {
                    InterPtType inter;

                    double v[3];
                    vtkMath::Subtract(endPt, pt, v);
                    inter.t = vtkMath::Dot(v, r);

                    double w[3];
                    Cpy(w, r, 3);

                    vtkMath::MultiplyScalar(w, inter.t);

                    vtkMath::Add(pt, w, inter.pt);

                    inter.count++;

                    inter.end = i;

                    inter.src = src;

                    itr3 = std::lower_bound(interPts2.begin(), interPts2.end(), inter);
                    interPts2.insert(itr3, inter);

                }
            }
        }

        std::map<int, int> locs;

        for (itr = interPts2.begin(); itr != interPts2.end(); ++itr) {
            if (itr->end != NO_USE) {
                locs[itr->end] = itr-interPts2.begin();
            }
        }

        // punkte die durch jeweils eine kongruente kante begrenzt sind

        for (int i = 0; i < numPts; i++) {
            int j = (i+1)%numPts;

            if (locs.count(i) == 1 && locs.count(j) == 1) {
                double ptA[3], ptB[3];

                pts->GetPoint(poly->GetId(i), ptA);
                pts->GetPoint(poly->GetId(j), ptB);

                double e[3];

                vtkMath::Subtract(ptB, ptA, e);
                vtkMath::Normalize(e);

                vtkMath::Cross(e, n, p);
                vtkMath::Normalize(p);

                // p zeigt nach außen

                d = vtkMath::Dot(p, ptA);

                // erstes ende

                if (locs[i] > 0 && locs[i] < interPts2.size()-1
                    && interPts2[locs[i]].count == 1) {

                    int prev = poly->GetId((i+numPts-1)%numPts);

                    double prevPt[3];
                    pts->GetPoint(prev, prevPt);

                    double prevD = vtkMath::Dot(p, prevPt)-d;

                    if (prevD > 0) {
                        interPts2[locs[i]].count++;
                    }

                }

                // zweites ende

                if (locs[j] > 0 && locs[j] < interPts2.size()-1
                    && interPts2[locs[j]].count == 1) {

                    int next = poly->GetId((j+1)%numPts);

                    double nextPt[3];
                    pts->GetPoint(next, nextPt);

                    double nextD = vtkMath::Dot(p, nextPt)-d;

                    if (nextD > 0) {
                        interPts2[locs[j]].count++;
                    }

                }

            }
        }

        InterPtsType::iterator itr4;

        for (itr4 = interPts2.begin(); itr4 != interPts2.end(); ++itr4) {
            if (itr4->end != NO_USE) {
                pts->GetPoint(poly->GetId(itr4->end), itr4->pt);
            }

            for (int i = 0; i < itr4->count; i++) {
                interPts3.push_back(*itr4);
            }

#ifdef DEBUG
            if (itr4->count > 1) {
                std::cout << "ind " << itr4->ind << ", count " << itr4->count << std::endl;
            }
#endif
        }

    }

    // wenn alles richtig ist, müsste die anzahl gerade sein

    return interPts3;

}

void vtkPolyDataContactFilter::InterPolys (vtkIdType idA, vtkIdType idB) {

#ifdef DEBUG
    std::cout << "InterPolys() -> idA " << idA << ", idB " << idB << std::endl;
#endif

    vtkIdList *polyA = vtkIdList::New();
    vtkIdList *polyB = vtkIdList::New();

    pdA->GetCellPoints(idA, polyA);
    pdB->GetCellPoints(idB, polyB);

    // ebenen aufstellen

    double nA[3], nB[3], ptA[3], ptB[3], dA, dB;

    ComputeNormal(pdA->GetPoints(), nA, polyA);
    ComputeNormal(pdB->GetPoints(), nB, polyB);

    pdA->GetPoint(polyA->GetId(0), ptA);
    pdB->GetPoint(polyB->GetId(0), ptB);

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

        InterPtsType intersA(InterPolyLine(pdA->GetPoints(), polyA, r, s, _Src::A, idA));
        InterPtsType intersB(InterPolyLine(pdB->GetPoints(), polyB, r, s, _Src::B, idB));

#ifdef DEBUG
        std::cout << "intersA " << intersA.size()
            << ", intersB " << intersB.size()
            << std::endl;
#endif

        if (intersA.size() != 0 && intersB.size() != 0
            && intersA.size()%2 == 0 && intersB.size()%2 == 0) {

            OverlapsType overlaps(OverlapLines(intersA, intersB));

            OverlapsType::const_iterator itr;

            for (itr = overlaps.begin(); itr != overlaps.end(); ++itr) {

                auto &f = itr->first,
                    &s = itr->second;

#ifdef DEBUG
                std::cout << "first " << f.ind
                    << ", second " << s.ind
                    << std::endl;
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

    polyB->Delete();
    polyA->Delete();

}

OverlapsType vtkPolyDataContactFilter::OverlapLines (InterPtsType &intersA, InterPtsType &intersB) {

    auto Add = [](InterPtType &a, InterPtType &b, InterPtType &c, InterPtType &d) {
        auto p = std::make_pair(a, b);

        p.first.Merge(c);
        p.second.Merge(d);

        return p;
    };

    OverlapsType ols;

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

    return ols;

}

int vtkPolyDataContactFilter::InterOBBNodes (vtkOBBNode *nodeA, vtkOBBNode *nodeB, vtkMatrix4x4 *mat, void *caller) {
    vtkPolyDataContactFilter *self = reinterpret_cast<vtkPolyDataContactFilter*>(caller);

    vtkIdList *cellsA = nodeA->Cells;
    vtkIdList *cellsB = nodeB->Cells;

    int numCellsA = cellsA->GetNumberOfIds(),
        numCellsB = cellsB->GetNumberOfIds();

    for (int i = 0; i < numCellsA; i++) {
        vtkIdType ci = cellsA->GetId(i);

        for (int j = 0; j < numCellsB; j++) {
            vtkIdType cj = cellsB->GetId(j);

            self->InterPolys(ci, cj);
        }
    }

    return 0;
}
