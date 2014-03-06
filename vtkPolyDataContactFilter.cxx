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

#include <vtkCellArray.h>

#include "vtkPolyDataContactFilter.h"
#include "GeomHelper.h"

vtkStandardNewMacro(vtkPolyDataContactFilter);

#define CPY(a, b) a[0] = b[0]; a[1] = b[1]; a[2] = b[2];

vtkPolyDataContactFilter::vtkPolyDataContactFilter () {

    contLines = vtkPolyData::New();
    contLines->Allocate(1000);

    contPts = vtkPoints::New();
    contLines->SetPoints(contPts);

    contA = vtkIntArray::New();
    contB = vtkIntArray::New();

    contA->SetName("cA");
    contB->SetName("cB");

    MergeLines = true;

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

        obbA->IntersectWithOBBTree(obbB, mat, IntersectOBBNodes, this);

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

        if (MergeLines) {

            resultA->GetCellData()->RemoveArray("cA");
            resultA->GetCellData()->RemoveArray("cB");

            resultA->BuildLinks();

            vtkIdList *neigs = vtkIdList::New();
            vtkIdList *line = vtkIdList::New();

            toRemove.clear();

            for (int i = 0; i < resultA->GetNumberOfPoints(); i++) {
                resultA->GetPointCells(i, neigs);

                std::vector<std::pair<int, int> > ends;

                for (int j = 0; j < neigs->GetNumberOfIds(); j++) {
                    //if (resultA->GetCellType(neigs->GetId(j)) != VTK_EMPTY_CELL) {

                    if (std::count(toRemove.begin(), toRemove.end(), neigs->GetId(j)) == 0) {

                        resultA->GetCellPoints(neigs->GetId(j), line);

                        ends.push_back(std::make_pair(i == line->GetId(1) ? line->GetId(0) : line->GetId(1), neigs->GetId(j)));

                        line->Reset();
                    }
                }

                std::vector<std::pair<int, int> >::const_iterator itr;

                bool found = false;

                for (itr = ends.begin()+1; itr != ends.end(); itr++) {
                    if (!found && ends.front().first != itr->first) {
                        found = true;
                    } else {
                        //resultA->DeleteCell(itr->second);

                        toRemove.push_back(itr->second);
                    }

                }

                neigs->Reset();

            }

            line->Delete();
            neigs->Delete();

            // das löschen ausführen
            //resultA->RemoveDeletedCells();

            GeomHelper::RemoveCells(resultA, toRemove);


        }

        resultB->DeepCopy(pdA);
        resultC->DeepCopy(pdB);

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

    //pd->RemoveDeletedCells();

    GeomHelper::RemoveCells(pd, toRemove);

}

InterPtType vtkPolyDataContactFilter::IntersectEdgeAndLine (double *edgePtA, double *edgePtB, double *r, double *pt) {

    InterPtType inter;

    // richtungsvektor der kante bestimmen

    double e[3];
    vtkMath::Subtract(edgePtB, edgePtA, e);
    double le = vtkMath::Normalize(e);

    // schnittpunkt ermitteln

    double v[3];
    vtkMath::Cross(r, e, v);

    double nV = vtkMath::Norm(v);

    if (nV*nV > 1e-6) {
        double p[3];

        vtkMath::Subtract(edgePtA, pt, p);

        double t = vtkMath::Determinant3x3(p, e, v)/(nV*nV);
        double s = vtkMath::Determinant3x3(p, r, v)/(nV*nV);

        double sPt[] = {
            pt[0]+t*r[0],
            pt[1]+t*r[1],
            pt[2]+t*r[2]
        };

        CPY(inter.pt, sPt)

        if (std::fabs(s) < 1e-6) {
            inter.onEnd = 0;
            inter.isOnEdge = true;
            inter.t = t;

        } else if (std::fabs(le-s) < 1e-6) {
            inter.onEnd = 1;
            inter.isOnEdge = true;
            inter.t = t;

        } else if (s > 0 && s < le) {

            // der schnittpunkt liegt irgendwo auf der kante

            inter.isOnEdge = true;
            inter.t = t;

        }

    } else {
        // parallel
    }

    return inter;

}

InterPtsType vtkPolyDataContactFilter::IntersectPolyAndLine (vtkPoints *pts, vtkIdList *poly, double *r, double *pt) {

    InterPtsType interPts;

    // durchläuft die kanten und ermittelt die schnittpunkte, die in iterPts gesammelt werden

    int edgeNbr = poly->GetNumberOfIds();

    for (unsigned int i = 0; i < edgeNbr; i++) {
        // kante ermitteln

        double ptA[3], ptB[3];

        pts->GetPoint(poly->GetId(i), ptA);
        pts->GetPoint(poly->GetId((i+1)%poly->GetNumberOfIds()), ptB);

        // schnittpunkt

        InterPtType inter(IntersectEdgeAndLine(ptA, ptB, r, pt));

        if (inter.isOnEdge) {
            if (inter.onEnd == 0 && i == 0) {
                edgeNbr--;
            } else if (inter.onEnd == 1) {
                i++;
            }

            interPts.push_back(inter);

        }

    }

    std::sort(interPts.begin(), interPts.end());

    return interPts;

}

void vtkPolyDataContactFilter::IntersectPolys (vtkIdType idA, vtkIdType idB) {

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

        InterPtsType intersA(IntersectPolyAndLine(pdA->GetPoints(), polyA, r, s));
        InterPtsType intersB(IntersectPolyAndLine(pdB->GetPoints(), polyB, r, s));

        if (intersA.size() > 1 && intersB.size() > 1) {

            OverlapsType overlaps(OverlapLines(intersA, intersB));

            OverlapsType::const_iterator itr;

            for (itr = overlaps.begin(); itr != overlaps.end(); itr++) {

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

    unsigned int i = 0, j;

    while (i < intersA.size()-1) {

        InterPtType &fA = intersA[i];
        InterPtType &sA = intersA[i+1];

        j = 0;

        while (j < intersB.size()-1) {
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


            if (sB.onEnd != -1) {
                j++;
            } else {
                j += 2;
            }

        }

        if (sA.onEnd != -1) {
            i++;
        } else {
            i += 2;
        }

    }

    return ols;

}

int vtkPolyDataContactFilter::IntersectOBBNodes (vtkOBBNode *nodeA, vtkOBBNode *nodeB, vtkMatrix4x4 *mat, void *caller) {
    vtkPolyDataContactFilter *self = reinterpret_cast<vtkPolyDataContactFilter*>(caller);

    vtkIdList *cellsA = nodeA->Cells;
    vtkIdList *cellsB = nodeB->Cells;

    for (int i = 0; i < cellsA->GetNumberOfIds(); i++) {
        vtkIdType ci = cellsA->GetId(i);

        for (int j = 0; j < cellsB->GetNumberOfIds(); j++) {
            vtkIdType cj = cellsB->GetId(j);

            self->IntersectPolys(ci, cj);
        }
    }
}
