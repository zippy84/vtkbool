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

#include <cmath>
#include <vector>
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
#include "GeomHelper.h"

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

    SetNumberOfInputPorts(2);
    SetNumberOfOutputPorts(3);
}


vtkPolyDataContactFilter::~vtkPolyDataContactFilter () {
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

        vtkCleanPolyData *clean = vtkCleanPolyData::New();

#if (VTK_MAJOR_VERSION == 5)
        clean->SetInput(contLines);
#else
        clean->SetInputData(contLines);
#endif
        clean->ToleranceIsAbsoluteOn();
        clean->SetAbsoluteTolerance(1e-5);
        clean->Update();

        resultA->DeepCopy(clean->GetOutput());

        std::vector<int> toRemove;

        for (unsigned int i = 0; i < resultA->GetNumberOfCells(); i++) {
            if (resultA->GetCellType(i) != VTK_LINE) {

                //resultA->DeleteCell(i);

                toRemove.push_back(i);

            }

        }

        //resultA->RemoveDeletedCells();

        GeomHelper::RemoveCells(resultA, toRemove);

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

    vtkIdType cellNbr = pd->GetNumberOfCells();

    vtkIntArray *cellIds = vtkIntArray::New();

    vtkIntArray *stripIds = vtkIntArray::New();

    for (unsigned int i = 0; i < cellNbr; i++) {
        cellIds->InsertNextValue(i);

        if (pd->GetCellType(i) == VTK_TRIANGLE_STRIP) {
            stripIds->InsertNextValue(i);
        }
    }

    vtkCellArray *cells = vtkCellArray::New();

    vtkCellArray *strips = pd->GetStrips();

    vtkIdType n;
    vtkIdType *pts;

    unsigned int i = 0;

    for (strips->InitTraversal(); strips->GetNextCell(n, pts);) {
        cells->Reset();

        vtkTriangleStrip::DecomposeStrip(n, pts, cells);

        for (cells->InitTraversal(); cells->GetNextCell(n, pts);) {
            pd->InsertNextCell(VTK_TRIANGLE, n, pts);
            cellIds->InsertNextValue(stripIds->GetValue(i));

            cellNbr++;
        }

        i++;

    }

    std::vector<int> toRemove;

    int type;

    for (unsigned int i = 0; i < cellNbr; i++) {
        type = pd->GetCellType(i);

        if (type != VTK_POLYGON && type != VTK_QUAD && type != VTK_TRIANGLE) {

            //pd->DeleteCell(i);

            toRemove.push_back(i);
        }

    }

    cellIds->SetName("OrigCellIds");

    pd->GetCellData()->SetScalars(cellIds);

    cells->Delete();
    stripIds->Delete();
    cellIds->Delete();

    vtkDoubleArray *doublePts = vtkDoubleArray::New();
    doublePts->DeepCopy(pd->GetPoints()->GetData());

    pd->GetPoints()->SetData(doublePts);
    doublePts->Delete();

    //pd->RemoveDeletedCells();

    GeomHelper::RemoveCells(pd, toRemove);

}

#ifdef DEBUG
int ind = 0;
#endif

InterPtType vtkPolyDataContactFilter::InterEdgeLine (double *eA, double *eB, double *r, double *pt) {

    InterPtType inter;

    // richtungsvektor der kante bestimmen

    double e[3];
    vtkMath::Subtract(eB, eA, e);
    double l = vtkMath::Normalize(e);

    // schnittpunkt ermitteln

    double v[3];
    vtkMath::Cross(r, e, v);

    double n = vtkMath::Norm(v);

    if (n*n > 1e-9) {
        double p[3];

        vtkMath::Subtract(eA, pt, p);

        double s = vtkMath::Determinant3x3(p, r, v)/(n*n);

        if (s > -1e-6 && s < l+1e-6) {
            double t = vtkMath::Determinant3x3(p, e, v)/(n*n);

            inter.pt[0] = pt[0]+t*r[0];
            inter.pt[1] = pt[1]+t*r[1];
            inter.pt[2] = pt[2]+t*r[2];

            inter.onEdge = true;

            inter.t = t;

            inter.outer = s < 0 || s > l;

#ifdef DEBUG
            inter.ind = ind++;
#endif

        }

    } else {
        // parallel
    }

    return inter;

}

InterPtsType vtkPolyDataContactFilter::InterPolyLine (vtkPoints *pts, vtkIdList *poly, double *r, double *pt) {

#ifdef DEBUG
    std::cout << "InterPolyLine()" << std::endl;
#endif

    InterPtsType interPts;

    // durchläuft die kanten und ermittelt die schnittpunkte

    for (unsigned int i = 0; i < poly->GetNumberOfIds(); i++) {
        // kante ermitteln

        double ptA[3], ptB[3];

        pts->GetPoint(poly->GetId(i), ptA);
        pts->GetPoint(poly->GetId((i+1)%poly->GetNumberOfIds()), ptB);

        // schnittpunkt

        InterPtType inter(InterEdgeLine(ptA, ptB, r, pt));

        if (inter.onEdge) {
            interPts.push_back(inter);
        }

    }

    std::sort(interPts.begin(), interPts.end());

    // jetzt müssen noch spezielle gelöscht werden

    InterPtsType interPts2;

    InterPtsType::const_iterator itr;

#ifdef DEBUG
    for (itr = interPts.begin(); itr != interPts.end(); itr++) {
        std::cout << "ind " << itr->ind
            << ", pt [" << itr->pt[0] << ", " << itr->pt[1] << ", " << itr->pt[2] << "]"
            << ", onEdge " << itr->onEdge
            << ", t " << itr->t
            << ", outer " << itr->outer
            << std::endl;
    }
#endif

    if (interPts.size() > 1) {

        std::vector<const InterPtType*> omits;

        for (itr = interPts.begin(); itr != interPts.end()-1; itr++) {

            double v[3];
            vtkMath::Subtract(itr->pt, (itr+1)->pt, v);

            if (vtkMath::Norm(v) < 1e-5) {
                if (itr->outer != (itr+1)->outer) {
                    if (itr->outer) {
                        omits.push_back(&*itr);
                    } else {
                        omits.push_back(&*(itr+1));
                    }
                } else {
                    // beide
                    omits.push_back(&*itr);
                    omits.push_back(&*(itr+1));
                }

            }
        }

        for (itr = interPts.begin(); itr != interPts.end(); itr++) {
            if (std::count(omits.begin(), omits.end(), &*itr) == 0) {
                interPts2.push_back(*itr);
            } else {
#ifdef DEBUG
                std::cout << "omit " << itr->ind << std::endl;
#endif

            }
        }

    }

    // wenn alles richtig ist, müsste die anzahl gerade sein

    return interPts2;

}

void vtkPolyDataContactFilter::InterPolys (vtkIdType idA, vtkIdType idB) {

#ifdef DEBUG
    std::cout << "InterPolys()" << std::endl;

    std::cout << "idA " << idA << ", idB " << idB << std::endl;
#endif

    vtkIdList *polyA = vtkIdList::New();
    vtkIdList *polyB = vtkIdList::New();

    pdA->GetCellPoints(idA, polyA);
    pdB->GetCellPoints(idB, polyB);

    // ebenen aufstellen

    double nA[3], nB[3], ptA[3], ptB[3], dA, dB;

    GeomHelper::ComputeNormal(pdA->GetPoints(), polyA, nA);
    GeomHelper::ComputeNormal(pdB->GetPoints(), polyB, nB);

    pdA->GetPoint(polyA->GetId(0), ptA);
    pdB->GetPoint(polyB->GetId(0), ptB);

    dA = vtkMath::Dot(nA, ptA);
    dB = vtkMath::Dot(nB, ptB);

    // sind die ebenen parallel?

    double p = std::fabs(vtkMath::Dot(nA, nB));

    if (p < 0.999999) {

        // richtungsvektor

        double r[3];
        vtkMath::Cross(nA, nB, r);
        vtkMath::Normalize(r);

        // lsg. des lin. gls. mittels cramersche regel

        unsigned int i = 0;

         for (unsigned int j = 1; j < 3; j++) {
            if (std::fabs(r[j]) > std::fabs(r[i])) {
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

        // für paraview
        std::cout << "s2 [" << (s[0]+10*r[0])
            << ", " << (s[1]+10*r[1])
            << ", " << (s[2]+10*r[2])
            << "]" << std::endl;


#endif

        InterPtsType intersA(InterPolyLine(pdA->GetPoints(), polyA, r, s));
        InterPtsType intersB(InterPolyLine(pdB->GetPoints(), polyB, r, s));

#ifdef DEBUG
        std::cout << "intersA " << intersA.size()
            << ", intersB " << intersB.size()
            << std::endl;
#endif

        if (intersA.size() != 0 && intersB.size() != 0
            && intersA.size()%2 == 0 && intersB.size()%2 == 0) {

            OverlapsType overlaps(OverlapLines(intersA, intersB));

            OverlapsType::const_iterator itr;

            for (itr = overlaps.begin(); itr != overlaps.end(); itr++) {

#ifdef DEBUG
                std::cout << "first " << itr->first.ind
                    << ", second " << itr->second.ind
                    << std::endl;
#endif

                vtkIdList *linePts = vtkIdList::New();

                linePts->InsertNextId(contPts->InsertNextPoint(itr->first.pt));
                linePts->InsertNextId(contPts->InsertNextPoint(itr->second.pt));

                contLines->InsertNextCell(VTK_LINE, linePts);

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

    OverlapsType ols;

    for (unsigned int i = 0; i < intersA.size(); i += 2) {

        InterPtType &fA = intersA[i];
        InterPtType &sA = intersA[i+1];

        for (unsigned int j = 0; j < intersB.size(); j += 2) {
            InterPtType &fB = intersB[j];
            InterPtType &sB = intersB[j+1];

            std::pair<InterPtType, InterPtType> ol;

            if (fA.t <= fB.t && sA.t > fB.t) {
                ol.first = fB;

                if (sB.t < sA.t) {
                    ol.second = sB;
                } else {
                    ol.second = sA;
                }

                ols.push_back(ol);

            } else if (fB.t <= fA.t && sB.t > fA.t) {
                ol.first = fA;

                if (sA.t < sB.t) {
                    ol.second = sA;
                } else {
                    ol.second = sB;
                }

                ols.push_back(ol);
            }

        }

    }

    return ols;

}

int vtkPolyDataContactFilter::InterOBBNodes (vtkOBBNode *nodeA, vtkOBBNode *nodeB, vtkMatrix4x4 *mat, void *caller) {
    vtkPolyDataContactFilter *self = reinterpret_cast<vtkPolyDataContactFilter*>(caller);

    vtkIdList *cellsA = nodeA->Cells;
    vtkIdList *cellsB = nodeB->Cells;

    for (int i = 0; i < cellsA->GetNumberOfIds(); i++) {
        vtkIdType ci = cellsA->GetId(i);

        for (int j = 0; j < cellsB->GetNumberOfIds(); j++) {
            vtkIdType cj = cellsB->GetId(j);

            self->InterPolys(ci, cj);
        }
    }

    return 0;
}
