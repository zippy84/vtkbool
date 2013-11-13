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

#include <map>
#include <deque>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <functional>

// siehe std::bind
using namespace std::placeholders;

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkObjectFactory.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTrivialProducer.h>
#include <vtkMath.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkAppendPolyData.h>
#include <vtkKdTreePointLocator.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>

#include "vtkPolyDataBooleanFilter.h"
#include "vtkPolyDataContactFilter.h"
#include "GeomHelper.h"

#define CPY(a, b) a[0] = b[0]; a[1] = b[1]; a[2] = b[2];

#ifdef DEBUG
#include <vtkPolyDataWriter.h>
#endif

vtkStandardNewMacro(vtkPolyDataBooleanFilter);

vtkPolyDataBooleanFilter::vtkPolyDataBooleanFilter () {

    SetNumberOfInputPorts(2);
    SetNumberOfOutputPorts(2);

    OperMode = OPER_UNION;

}

vtkPolyDataBooleanFilter::~vtkPolyDataBooleanFilter () {
    // derzeit nix tun
}

int vtkPolyDataBooleanFilter::ProcessRequest(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

    if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA())) {

        vtkInformation *inInfoA = inputVector[0]->GetInformationObject(0);
        vtkInformation *inInfoB = inputVector[1]->GetInformationObject(0);

        vtkPolyData *pdA = vtkPolyData::SafeDownCast(inInfoA->Get(vtkDataObject::DATA_OBJECT()));
        vtkPolyData *pdB = vtkPolyData::SafeDownCast(inInfoB->Get(vtkDataObject::DATA_OBJECT()));

        vtkInformation *outInfoA = outputVector->GetInformationObject(0);
        vtkInformation *outInfoB = outputVector->GetInformationObject(1);

        resultA = vtkPolyData::SafeDownCast(outInfoA->Get(vtkDataObject::DATA_OBJECT()));
        vtkPolyData *resultB = vtkPolyData::SafeDownCast(outInfoB->Get(vtkDataObject::DATA_OBJECT()));

        // eventuell vorhandene regionen vereinen

        vtkCleanPolyData *cleanA = vtkCleanPolyData::New();
#if (VTK_MAJOR_VERSION == 5)
        cleanA->SetInput(pdA);
#else
        cleanA->SetInputData(pdA);
#endif
        cleanA->Update();

        vtkCleanPolyData *cleanB = vtkCleanPolyData::New();
#if (VTK_MAJOR_VERSION == 5)
        cleanB->SetInput(pdB);
#else
        cleanB->SetInputData(pdB);
#endif
        cleanB->Update();

        // ergebnis-vernetzungen

        modPdA = vtkPolyData::New();
        modPdA->DeepCopy(cleanA->GetOutput());

        modPdB = vtkPolyData::New();
        modPdB->DeepCopy(cleanB->GetOutput());

        GeomHelper::PreparePolyData(modPdA);
        GeomHelper::PreparePolyData(modPdB);

#ifdef DEBUG


        vtkPolyDataWriter *w = vtkPolyDataWriter::New();
        w->SetFileName("modPdA.vtk");

#if (VTK_MAJOR_VERSION == 5)
        w->SetInput(modPdA);
#else
        w->SetInputData(modPdA);
#endif

        std::cout << "Exporting modPdA.vtk" << std::endl;
        w->Update();

        w->SetFileName("modPdB.vtk");

#if (VTK_MAJOR_VERSION == 5)
        w->SetInput(modPdB);
#else
        w->SetInputData(modPdB);
#endif

        std::cout << "Exporting modPdB.vtk" << std::endl;
        w->Update();
        w->Delete();

#endif

        // ermittelt kontaktstellen

        vtkTrivialProducer *prA = vtkTrivialProducer::New();
        prA->SetOutput(modPdA);

        vtkTrivialProducer *prB = vtkTrivialProducer::New();
        prB->SetOutput(modPdB);

        vtkPolyDataContactFilter *cl = vtkPolyDataContactFilter::New();
        cl->SetInputConnection(0, prA->GetOutputPort());
        cl->SetInputConnection(1, prB->GetOutputPort());
        cl->MergeLinesOff();
        cl->Update();

        contLines = cl->GetOutput();

        if (contLines->GetNumberOfCells() == 0) {
            vtkErrorMacro("Inputs have no contact.");

            return 1;
        }

        // in den CellDatas steht drin, welche polygone einander schneiden

        vtkIntArray *contsA = reinterpret_cast<vtkIntArray*>(contLines->GetCellData()->GetScalars("cA"));
        vtkIntArray *contsB = reinterpret_cast<vtkIntArray*>(contLines->GetCellData()->GetScalars("cB"));

        PolyStripsType polyStripsA(GetPolyStrips(modPdA, contsA));
        PolyStripsType polyStripsB(GetPolyStrips(modPdB, contsB));

        // polyStripsA/B wird modifiziert
        StripsType allStripsA(GetAllStrips(polyStripsA));
        StripsType allStripsB(GetAllStrips(polyStripsB));

        // trennt die polygone an den linien
        CutCells(modPdA, polyStripsA);
        CutCells(modPdB, polyStripsB);

#ifdef DEBUG
        vtkPolyDataWriter *w2 = vtkPolyDataWriter::New();
        w2->SetFileName("modPdA_2.vtk");

#if (VTK_MAJOR_VERSION == 5)
        w2->SetInput(modPdA);
#else
        w2->SetInputData(modPdA);
#endif

        std::cout << "Exporting modPdA_2.vtk" << std::endl;
        w2->Update();

        w2->SetFileName("modPdB_2.vtk");

#if (VTK_MAJOR_VERSION == 5)
        w2->SetInput(modPdB);
#else
        w2->SetInputData(modPdB);
#endif

        std::cout << "Exporting modPdB_2.vtk" << std::endl;
        w2->Update();
        w2->Delete();

#endif

        AddAdjacentPoints(modPdA, allStripsA);
        AddAdjacentPoints(modPdB, allStripsB);

#ifdef DEBUG
        vtkPolyDataWriter *w3 = vtkPolyDataWriter::New();
        w3->SetFileName("modPdA_3.vtk");

#if (VTK_MAJOR_VERSION == 5)
        w3->SetInput(modPdA);
#else
        w3->SetInputData(modPdA);
#endif

        std::cout << "Exporting modPdA_3.vtk" << std::endl;
        w3->Update();

        w3->SetFileName("modPdB_3.vtk");

#if (VTK_MAJOR_VERSION == 5)
        w3->SetInput(modPdB);
#else
        w3->SetInputData(modPdB);
#endif

        std::cout << "Exporting modPdB_3.vtk" << std::endl;
        w3->Update();
        w3->Delete();

#endif

        DisjoinPolys(modPdA, allStripsA);
        DisjoinPolys(modPdB, allStripsB);

#ifdef DEBUG
        vtkPolyDataWriter *w4 = vtkPolyDataWriter::New();
        w4->SetFileName("modPdA_4.vtk");

#if (VTK_MAJOR_VERSION == 5)
        w4->SetInput(modPdA);
#else
        w4->SetInputData(modPdA);
#endif

        std::cout << "Exporting modPdA_4.vtk" << std::endl;
        w4->Update();

        w4->SetFileName("modPdB_4.vtk");

#if (VTK_MAJOR_VERSION == 5)
        w4->SetInput(modPdB);
#else
        w4->SetInputData(modPdB);
#endif

        std::cout << "Exporting modPdB_4.vtk" << std::endl;
        w4->Update();
        w4->Delete();

#endif


        MergePoints(modPdA, allStripsA);
        MergePoints(modPdB, allStripsB);

        CombineRegions();

        resultB->DeepCopy(contLines);

        /*
        vtkAppendPolyData *app = vtkAppendPolyData::New();
        app->AddInputConnection(prA->GetOutputPort());
        app->AddInputConnection(prB->GetOutputPort());
        app->Update();

        resultA->DeepCopy(app->GetOutput());

        app->Delete();
        */

        cl->Delete();
        prB->Delete();
        prA->Delete();
        modPdB->Delete();
        modPdA->Delete();
        cleanB->Delete();
        cleanA->Delete();

    }

    return 1;

}

StripsType vtkPolyDataBooleanFilter::GetStrips (vtkPolyData *pd, int polyInd, std::deque<int> &lines) {

#ifdef DEBUG
    std::cout << "GetStrips()"
        << std::endl << "polyInd " << polyInd << std::endl;

#endif

    StripsType strips;

    // zuerst wird die lage der enden im polygon analysiert

    std::map<int, StripPtType> stripPts;

    // aus folgenden ecken setzt sich das polygon zusammen

    vtkIdList *polyPts = vtkIdList::New();
    pd->GetCellPoints(polyInd, polyPts);

    RemoveDuplicates(contLines, lines);

    std::deque<int>::const_iterator itr;

    for (itr = lines.begin(); itr != lines.end(); itr++) {
        // def. der linie

        vtkIdList *linePts = vtkIdList::New();

        contLines->GetCellPoints(*itr, linePts);

        // diese punkte durchlaufen

        for (unsigned int i = 0; i < 2; i++) {

            int realInd = linePts->GetId(i);

            if (stripPts.count(realInd) == 0) {
                // lage analysieren

                stripPts[realInd].ind = realInd;

                // die koordinaten
                double pt[3];
                contLines->GetPoint(realInd, pt);

                CPY(stripPts[realInd].pt, pt)

                double nSum = 0;

                // jetzt muss man die kanten durchlaufen
                for (unsigned int j = 0; j < polyPts->GetNumberOfIds(); j++) {
                    // richtungsvektor der kante bestimmen

                    int indA, indB;
                    double a[3], b[3];

                    indA = polyPts->GetId(j);
                    indB = polyPts->GetId((j+1)%polyPts->GetNumberOfIds());

                    pd->GetPoint(indA, a);
                    pd->GetPoint(indB, b);

                    double sA[3], sB[3];
                    vtkMath::Subtract(a, pt, sA);
                    vtkMath::Subtract(b, pt, sB);

                    // richtungsvektor und länge der kante

                    double u[3];

                    vtkMath::Subtract(b, a, u);

                    double n = vtkMath::Norm(u);

                    if (vtkMath::Norm(sA) < 1e-5) {
                        stripPts[realInd].onEdgeVert = true;
                        stripPts[realInd].edgeVert = indA;
                        stripPts[realInd].T = nSum;

                        break;

                    } else if (vtkMath::Norm(sB) < 1e-5) {
                        stripPts[realInd].onEdgeVert = true;
                        stripPts[realInd].edgeVert = indB;
                        stripPts[realInd].T = nSum+n;

                        break;

                    } else {
                        // liegt sie in der mitte der kante?

                        double v[3];
                        vtkMath::Subtract(pt, a, v);

                        double t = vtkMath::Dot(v, u)/(n*n);

                        double w[3];
                        vtkMath::Cross(v, u, w);

                        double d = vtkMath::Norm(w)/n;

                        // t ist bis hier eine fraktion von n

                        t *= n;

                        if (d < 1e-5 && t > 0 && t < n) {
                            stripPts[realInd].edgeVert = indA;
                            stripPts[realInd].secondVert = indB;
                            stripPts[realInd].t = t;
                            stripPts[realInd].T = nSum+t;

                            break;

                        }

                    }

                    nSum += n;

                }

            }

        }

        linePts->Delete();

    }

#ifdef DEBUG
    std::cout << "stripPts:" << std::endl;

    std::map<int, StripPtType>::const_iterator itr2;

    for (itr2 = stripPts.begin(); itr2 != stripPts.end(); itr2++) {

        std::cout << "ind " << itr2->second.ind
                  << ", edgeVert " << itr2->second.edgeVert
                  << ", onEdgeVert " << itr2->second.onEdgeVert
                  << ", T " << itr2->second.T
                  << ", t " << itr2->second.t
                  << std::endl;

    }
#endif

    // nun werden die linien mit hilfe der stripPts zusammengefügt

    StripType strip;

    int i = 0;

    while (lines.size() > 0) {

        vtkIdList *linePts = vtkIdList::New();
        contLines->GetCellPoints(lines[i], linePts);

        int fInd = linePts->GetId(0);
        int sInd = linePts->GetId(1);

        if (strip.empty()) {
            strip.push_back(stripPts[fInd]);
            strip.push_back(stripPts[sInd]);

            lines.erase(lines.begin());
        } else {
            if (strip.back().edgeVert == -1 && strip.back().ind == fInd) {
                strip.push_back(stripPts[sInd]);
                lines.erase(lines.begin()+i);
                i = 0;

            } else if (strip.back().edgeVert == -1 && strip.back().ind == sInd) {
                strip.push_back(stripPts[fInd]);
                lines.erase(lines.begin()+i);
                i = 0;

            } else if (strip.front().edgeVert == -1 && strip.front().ind == fInd) {
                strip.push_front(stripPts[sInd]);
                lines.erase(lines.begin()+i);
                i = 0;

            } else if (strip.front().edgeVert == -1 && strip.front().ind == sInd) {
                strip.push_front(stripPts[fInd]);
                lines.erase(lines.begin()+i);
                i = 0;

            } else {
                i++;

            }
        }

        // ein strip ist abgeschlossen, wenn beide enden auf einer kante liegen
        // oder wenn beide enden gleich sind
        // dann stellt der strip ein eigenes polygon dar

        if ((strip.front().edgeVert > -1 && strip.back().edgeVert > -1)
            || (strip.front().ind == strip.back().ind)
            || (i == lines.size())) {

            // den strip nach dem T sortieren
            if (strip.front().T > strip.back().T) {
                std::reverse(strip.begin(), strip.end());
            }

            // einen neuen strip anlegen
            strips.push_back(strip);

            strip.clear();

            i = 0;
        }

        linePts->Delete();

    }

    polyPts->Delete();

    CompleteStrips(strips);

#ifdef DEBUG
    std::cout << "strips:" << std::endl;

    StripsType::iterator itr3;

    for (itr3 = strips.begin(); itr3 != strips.end(); itr3++) {
        itr3->front().pos = std::distance(strips.begin(), itr3);

        std::cout << "pos " << itr3->front().pos
            << ", front " << itr3->front().ind
            << ", back " << itr3->back().ind
            << ", size " << itr3->size()
            << std::endl;

    }

#endif

    return strips;

}

PolyStripsType vtkPolyDataBooleanFilter::GetPolyStrips (vtkPolyData *pd, vtkIntArray *conts) {

#ifdef DEBUG
    std::cout << "GetPdStrips()" << std::endl;
#endif

    std::map<int, StripsType> strips;

    // zunächst findet die zuordnung der linien zu den polygonen statt

    std::map<int, std::deque<int> > polyLines;

    if (conts != NULL) {

        for (int i = 0; i < conts->GetNumberOfTuples(); i++) {

            int poly = conts->GetValue(i);

            polyLines[poly].push_back(i);
        }
    }

    std::map<int, std::deque<int> >::iterator itr;

    for (itr = polyLines.begin(); itr != polyLines.end(); itr++) {

        strips[itr->first] = GetStrips(pd, itr->first, itr->second);

    }

    return strips;

}

StripsType vtkPolyDataBooleanFilter::GetAllStrips (PolyStripsType &polyStrips) {
    StripsType all;

    PolyStripsType::const_iterator itr;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); itr++) {
        all.insert(all.end(), itr->second.begin(), itr->second.end());
    }

    return all;

}


bool vtkPolyDataBooleanFilter::SortFct (const StripType &stripA, const StripType &stripB, double *n) {

    if (stripA.front().ind != stripB.front().ind) {

        return stripA.front().T < stripB.front().T;

    } else {
        // in umgekehrter reihenfolge

        if (stripA.back().ind != stripB.back().ind) {
            return stripA.back().T > stripB.back().T;

        } else {

            // anfang und ende sind gleich

            // hier muss man mit der normale arbeiten

            vtkPoints *pts = vtkPoints::New();

            vtkIdList *poly = vtkIdList::New();

            StripType::const_iterator itr;

            for (itr = stripA.begin(); itr != stripA.end(); itr++) {
                poly->InsertNextId(pts->InsertNextPoint(itr->pt));
            }

            StripType::const_reverse_iterator itr2;

            for (itr2 = stripB.rbegin()+1; itr2 != stripB.rend()-1; itr2++) {
                poly->InsertNextId(pts->InsertNextPoint(itr2->pt));
            }

            double s[3];

            GeomHelper::ComputeNormal(pts, poly, s);

            // diese normale muss man jetzt mit der normale des polygons vergleichen in dem der strip liegt

            double ang = vtkMath::Dot(n, s);

            poly->Delete();
            pts->Delete();

            return ang < 0;

        }
    }

}



void vtkPolyDataBooleanFilter::CutCells (vtkPolyData *pd, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "CutCells()" << std::endl;
#endif

    vtkPoints *pdPts = pd->GetPoints();

    std::vector<int> toRemove;

    std::map<int, StripsType>::iterator itr;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); itr++) {

        // die zu teilenden polygone
        std::deque<int> polys;
        polys.push_back(itr->first);

        // die normale des polygons wird benötigt

        vtkIdList *orig = vtkIdList::New();
        pd->GetCellPoints(itr->first, orig);

        double n[3];
        GeomHelper::ComputeNormal(pd->GetPoints(), orig, n);

        orig->Delete();

#ifdef DEBUG
        std::cout << "orig " << itr->first << std::endl;
#endif

        StripsType &strips = itr->second;

        // sortieren der strips mit einer speziellen funktion
        // bindet die normale an die SortFct, sodass die funktion die geforderten zwei argument besitzt

        std::function<bool(const StripType&, const StripType&)> _SortFct = std::bind(vtkPolyDataBooleanFilter::SortFct, _1, _2, n);

        std::sort(strips.begin(), strips.end(), _SortFct);

        StripsType::iterator itr2;

        for (itr2 = strips.begin(); itr2 != strips.end(); itr2++) {

#ifdef DEBUG
            std::cout << "strip:" << std::endl
                << "pos " << itr2->front().pos
                << std::endl;
#endif

            // map für die ecken die ersetzt werden müssen

            std::map<int, ReplsType> R;

            int Q = -1;

            if (itr2->front().edgeVert == -1 && itr2->back().edgeVert == -1) {
                // der strip bildet ein polygon im inneren

#ifdef DEBUG
                std::cout << "is a hole" << std::endl;
#endif

                continue;
            }

            int cycle = 0;

            while (true) {

                if (cycle == polys.size()) {
                    // es wurde ein fall nicht berücksichtigt ...

                    break;
                }

                int polyInd = polys.front();
                polys.pop_front();

                vtkIdList *poly = vtkIdList::New();
                pd->GetCellPoints(polyInd, poly);

#ifdef DEBUG
                std::cout << "polyInd " << polyInd << ", poly [";

                for (unsigned int i = 0; i < poly->GetNumberOfIds(); i++) {
                    std::cout << poly->GetId(i) << ", ";
                }

                std::cout << "]" << std::endl;
#endif

                // zwei neue polygone
                vtkIdList *newPolyA = vtkIdList::New();
                vtkIdList *newPolyB = vtkIdList::New();

                if (itr2->front().edgeVert == itr2->back().edgeVert) {

#ifdef DEBUG
                    std::cout << "case 1" << std::endl;
#endif

                    // der strip bildet ein eigenständiges polygon auf einer kante

                    for (unsigned int i = 0; i < poly->GetNumberOfIds(); i++) {
                        int ptInd = poly->GetId(i);

                        if (ptInd == itr2->front().edgeVert) {
                            if (!itr2->front().onEdgeVert) {
                                newPolyA->InsertNextId(ptInd);
                            }

                            // den strip einfügen

                            StripType::iterator itr3;

                            for (itr3 = itr2->begin(); itr3 != itr2->end(); itr3++) {

                                // einen neuen punkt erzeugen
                                double pt[3];
                                contLines->GetPoint(itr3->ind, pt);

                                newPolyA->InsertNextId(pdPts->InsertNextPoint(pt));
                            }

                            // dann das eigenständige polygon bestehend aus dem strip

                            StripType::reverse_iterator itr4;

                            for (itr4 = itr2->rbegin(); itr4 != itr2->rend(); itr4++) {
                                double pt[3];
                                contLines->GetPoint(itr4->ind, pt);

                                newPolyB->InsertNextId(pdPts->InsertNextPoint(pt));
                            }

                            // eine neue referenz hinzufügen
                            R[ptInd].push_back(std::make_pair(itr2->front().t, newPolyB->GetId(newPolyB->GetNumberOfIds()-1)));
                            R[ptInd].push_back(std::make_pair(itr2->back().t, newPolyA->GetId(newPolyA->GetNumberOfIds()-1)));

                        } else {
                            newPolyA->InsertNextId(ptInd);
                        }
                    }

                } else {

#ifdef DEBUG
                    std::cout << "case 2" << std::endl;
#endif

                    vtkIdList *actPoly = newPolyA;

                    // kanten durchlaufen

                    int refs[] = {-1, -1};

                    for (unsigned int i = 0; i < poly->GetNumberOfIds(); i++) {
                        int ptInd = poly->GetId(i);

                        if (ptInd == itr2->front().edgeVert) {
                            // den strip einfügen

                            if (!itr2->front().onEdgeVert) {
                                actPoly->InsertNextId(ptInd);
                            }

                            StripType::iterator itr3;

                            for (itr3 = itr2->begin(); itr3 != itr2->end(); itr3++) {

                                // einen neuen punkt erzeugen
                                double pt[3];
                                contLines->GetPoint(itr3->ind, pt);

                                actPoly->InsertNextId(pdPts->InsertNextPoint(pt));
                            }

                            R[ptInd].push_back(std::make_pair(itr2->front().t, actPoly->GetId(actPoly->GetNumberOfIds()-1)));

                            refs[0] = ptInd;

                            actPoly = (actPoly == newPolyA) ? newPolyB : newPolyA;


                        } else if (ptInd == itr2->back().edgeVert) {

                            if (!itr2->back().onEdgeVert) {
                                actPoly->InsertNextId(ptInd);
                            }

                            StripType::reverse_iterator itr3;

                            for (itr3 = itr2->rbegin(); itr3 != itr2->rend(); itr3++) {

                                double pt[3];
                                contLines->GetPoint(itr3->ind, pt);

                                actPoly->InsertNextId(pdPts->InsertNextPoint(pt));
                            }

                            R[ptInd].push_back(std::make_pair(itr2->back().t, actPoly->GetId(actPoly->GetNumberOfIds()-1)));

                            refs[1] = ptInd;

                            Q = actPoly->GetId(actPoly->GetNumberOfIds()-itr2->size());

                            actPoly = (actPoly == newPolyA) ? newPolyB : newPolyA;

                        } else {

                            actPoly->InsertNextId(ptInd);
                        }

                    }

                    // die ersetzungen müssen getauscht werden

                    if (newPolyB->GetNumberOfIds() > 0) {
                        int tmp = R[refs[0]].front().second;
                        R[refs[0]].front().second = R[refs[1]].front().second;
                        R[refs[1]].front().second = tmp;

                    }

                }

#ifdef DEBUG
                std::cout << "newPolyA [";

                for (unsigned int i = 0; i < newPolyA->GetNumberOfIds(); i++) {
                    std::cout << newPolyA->GetId(i) << ", ";
                }

                std::cout << "]"
                    << std::endl << "newPolyB [";

                for (unsigned int i = 0; i < newPolyB->GetNumberOfIds(); i++) {
                    std::cout << newPolyB->GetId(i) << ", ";
                }

                std::cout << "]" << std::endl;

#endif

                if (newPolyB->GetNumberOfIds() > 0) {

                    // die alte entfernen und die beiden neuen hinzufügen

                    toRemove.push_back(polyInd);

                    // strips können lediglich aus zwei punkten bestehen und auf einer kante des polygons liegen

                    if (newPolyA->GetNumberOfIds() > 2) {
                        polys.push_front(pd->InsertNextCell(VTK_POLYGON, newPolyA));

#ifdef DEBUG
                        std::cout << "newPolyA added" << std::endl;
#endif


                    }

                    if (HasArea(*itr2)) {

                        /*
                        folgendes muss gemacht werden, wenn zwei strips im gleichen punkt enden

                        nachdem der erste strip behandelt wurd, endet der zweite strip auf jeden fall in einer ecke der entstandenen polygone

                        das bedeutet, dass das ende des zweiten strips, ohne einsatz von R, den falschen edgeVert besitzt, nämlich den der kante davor

                        der eigentlich edgeVert wird deshalb dem newPolyB hinzugefügt

                        dieser befindet sich an stelle 0 des neuen polygons, muss deshalb entfernt werden
                        */

                        double eA[3], eB[3];
                        pd->GetPoint(newPolyB->GetId(0), eA);
                        pd->GetPoint(newPolyB->GetId(newPolyB->GetNumberOfIds()-1), eB);

                        // wenn die koordinaten nahezu gleich sind, dann die erste id entfernen

                        double v[3];
                        vtkMath::Subtract(eA, eB, v);

                        if (vtkMath::Norm(v) < 1e-6) {
                            newPolyB->DeleteId(newPolyB->GetId(0));

#ifdef DEBUG
                        std::cout << "first point deleted" << std::endl;
#endif


                        }

                        if (newPolyB->GetNumberOfIds() > 2) {
                            polys.push_front(pd->InsertNextCell(VTK_POLYGON, newPolyB));

#ifdef DEBUG
                            std::cout << "newPolyB added" << std::endl;
#endif

                        }

                    }

                    // die neuen polygone zugänglich machen

                    pd->BuildCells();

                    cycle = 0;

                    // abbruch der while-schleife
                    break;

                } else {
                    // zurückstellen, da der strip nicht in diesem polygon liegt
                    polys.push_back(polyInd);

                    cycle++;
                }

                newPolyB->Delete();
                newPolyA->Delete();

            }

            // enden der verbliebenen strips aktualisieren
            // das T muss nicht aktualisiert werden

            StripsType::iterator itr3;

            for (itr3 = itr2+1; itr3 != strips.end(); itr3++) {
                if (R.count(itr3->front().edgeVert) > 0) {

#ifdef DEBUG
                    std::cout << "pos " << itr3->front().pos
                        << ", back().edgeVert " << itr3->back().edgeVert;
#endif

                    ReplsType &repls = R[itr3->front().edgeVert];

                    ReplsType::const_reverse_iterator itr4;

                    for (itr4 = repls.rbegin(); itr4 != repls.rend(); itr4++) {

                        if (itr3->front().t >= itr4->first) {

                            // wenn sich das ende auf der gleichen kante befindet, müssen auch dessen daten aktualisiert werden
                            // unabhängig von Q

                            if (itr3->back().edgeVert == itr3->front().edgeVert) {
                                itr3->back().edgeVert = itr4->second;
                                itr3->back().t -= itr4->first;
                            }

#ifdef DEBUG
                            std::cout << " -> " << itr3->back().edgeVert
                                << ", front().edgeVert " << itr3->front().edgeVert;
#endif

                            itr3->front().t -= itr4->first;
                            itr3->front().edgeVert = itr4->second;

#ifdef DEBUG
                            std::cout << " -> " << itr3->front().edgeVert << std::endl;
#endif

                            if (std::fabs(itr3->front().t) < 1e-7) {
                                itr3->front().t = 0;
                                itr3->front().onEdgeVert = true;
                            }

                            break;

                        }

                    }

                }

                if (Q > -1) {

                    double ptA[3], ptB[3];
                    pd->GetPoint(Q, ptA);

                    double dV[3], d;

                    contLines->GetPoint(itr3->back().ind, ptB);

                    vtkMath::Subtract(ptA, ptB, dV);

                    d = vtkMath::Norm(dV);

                    if (d < 1e-7) {

#ifdef DEBUG
                        std::cout << "back().edgeVert " << itr3->back().edgeVert << " -> " << Q << std::endl;
#endif

                        itr3->back().edgeVert = Q;

                        itr3->back().onEdgeVert = true;
                        itr3->back().t = 0;

                        Q = -1;

                    }

                }

            }

        }

    }

    std::vector<int>::const_iterator itr2;

    for (itr2 = toRemove.begin(); itr2 != toRemove.end(); itr2++) {
        pd->DeleteCell(*itr2);
    }

    // tatsächliches löschen
    pd->RemoveDeletedCells();

}

void vtkPolyDataBooleanFilter::RemoveDuplicates (vtkPolyData *pd, std::deque<int> &lines) {

    std::deque<int> uniqueLines;

    unsigned int i, j;

    // die indexe der enden auf übereinstimmung prüfen

    vtkIdList *linePtsA = vtkIdList::New();
    vtkIdList *linePtsB = vtkIdList::New();

    for (i = 0; i < lines.size()-1; i++) {

        j = i+1;

        pd->GetCellPoints(lines[i], linePtsA);

        while (j < lines.size()) {

            pd->GetCellPoints(lines[j], linePtsB);

            if ((linePtsA->GetId(0) == linePtsB->GetId(0) && linePtsA->GetId(1) == linePtsB->GetId(1)) ||
                (linePtsA->GetId(0) == linePtsB->GetId(1) && linePtsA->GetId(1) == linePtsB->GetId(0))) {

                // übereinstimmung

                break;
            }

            j++;
        }

        if (j == lines.size()) {
            // keine vorzeitige unterbrechung der while-schleife

            uniqueLines.push_back(lines[i]);
        }
    }

    uniqueLines.push_back(lines.back());

    linePtsA->Delete();
    linePtsB->Delete();

    lines.swap(uniqueLines);

}

void vtkPolyDataBooleanFilter::CompleteStrips (StripsType &strips) {
    for (unsigned int i = 0; i < strips.size(); i++) {
        if (strips[i].front().edgeVert != strips[i].back().edgeVert) {

            if (strips[i].front().edgeVert == -1) {
                strips[i].insert(strips[i].begin(), strips[i].rbegin(), strips[i].rend()-1);
            } else if (strips[i].back().edgeVert == -1) {
                strips[i].insert(strips[i].end(), strips[i].rbegin()+1, strips[i].rend());
            }
        }
    }
}

bool vtkPolyDataBooleanFilter::HasArea (StripType &strip) {
    bool area = true;

    unsigned int n = strip.size();

    if (n%2 == 1) {

        for (unsigned int i = 0; i < (n-1)/2; i++) {
            area = strip[i].ind != strip[n-i-1].ind;
        }
    }

    return area;
}

void vtkPolyDataBooleanFilter::AddAdjacentPoints (vtkPolyData *pd, StripsType &strips) {

#ifdef DEBUG
    std::cout << "AddAdjacentPoints()" << std::endl;
#endif

    std::map<std::pair<int, int>, std::vector<StripPtType> > adjacentPts;

    StripsType::const_iterator itr;

    for (itr = strips.begin(); itr != strips.end(); itr++) {
        // anfang und ende

        if (!itr->front().onEdgeVert && itr->front().edgeVert > -1) {
            adjacentPts[std::make_pair(itr->front().edgeVert, itr->front().secondVert)].push_back(itr->front());
        }

    }

    std::map<std::pair<int, int>, std::vector<StripPtType> >::iterator itr2;
    std::vector<StripPtType>::const_iterator itr3;

    for (itr2 = adjacentPts.begin(); itr2 != adjacentPts.end(); itr2++) {

        // InsertNextLinkedCell funktioniert leider nicht
        pd->BuildLinks();

        // jetzt das polygon ermitteln
        vtkIdList *polysA = vtkIdList::New();
        vtkIdList *polysB = vtkIdList::New();

        pd->GetPointCells(itr2->first.first, polysA);
        pd->GetPointCells(itr2->first.second, polysB);

        // sortieren nach t, allerdings in umgekehrter reihenfolge
        std::sort(itr2->second.rbegin(), itr2->second.rend());

        bool added = false;

        for (unsigned int i = 0; i < polysA->GetNumberOfIds() && !added; i++) {
            for (unsigned int j = 0; j < polysB->GetNumberOfIds() && !added; j++) {

                if (polysA->GetId(i) == polysB->GetId(j)
                    && pd->GetCellType(polysA->GetId(i)) != VTK_EMPTY_CELL) {

                    // es gibt nur eins oder keins

                    vtkIdList *poly = vtkIdList::New();

                    pd->GetCellPoints(polysA->GetId(i), poly);

                    vtkIdList *newPoly = vtkIdList::New();

                    for (unsigned int k = 0; k < poly->GetNumberOfIds(); k++) {

                        newPoly->InsertNextId(poly->GetId(k));

                        if (itr2->second.front().secondVert == poly->GetId(k)
                            && itr2->second.front().edgeVert == poly->GetId((k+1)%poly->GetNumberOfIds())) {

                            // die ursprüngliche kante ist noch vorhanden

                            int lastInd = -1;

                            for (itr3 = itr2->second.begin(); itr3 != itr2->second.end(); itr3++) {

                                if (lastInd < 0 || lastInd != itr3->ind) {

                                    newPoly->InsertNextId(pd->GetPoints()->InsertNextPoint(itr3->pt));

                                    lastInd = itr3->ind;
                                }

                            }

                            added = true;

                        }

                    }

                    // das neue muss jetzt hinzugefügt werden

                    pd->DeleteCell(polysA->GetId(i));
                    pd->InsertNextCell(VTK_POLYGON, newPoly);

                    newPoly->Delete();
                    poly->Delete();

                }

            }
        }

        polysA->Delete();
        polysB->Delete();

    }

    pd->RemoveDeletedCells();

}

void vtkPolyDataBooleanFilter::DisjoinPolys (vtkPolyData *pd, StripsType &strips) {

#ifdef DEBUG
    std::cout << "DisjoinPolys()" << std::endl;
#endif

    pd->BuildLinks();

    std::set<int> touchedPts;

    StripsType::const_iterator itr;

    for (itr = strips.begin(); itr != strips.end(); itr++) {
        if (itr->front().onEdgeVert) {
            touchedPts.insert(itr->front().ind);
        }

        if (itr->back().onEdgeVert) {
            touchedPts.insert(itr->back().ind);
        }
    }

    vtkKdTreePointLocator *pl = vtkKdTreePointLocator::New();
    pl->SetDataSet(pd);
    pl->BuildLocator();

    std::set<int>::const_iterator itr2;

    for (itr2 = touchedPts.begin(); itr2 != touchedPts.end(); itr2++) {
        vtkIdList *pts = vtkIdList::New();

        double pt[3];

        contLines->GetPoint(*itr2, pt);

        pl->FindPointsWithinRadius(1e-7, pt, pts);

        for (unsigned int i = 0; i < pts->GetNumberOfIds(); i++) {
            vtkIdList *polys = vtkIdList::New();

            pd->GetPointCells(pts->GetId(i), polys);

            if (polys->GetNumberOfIds() > 1) {

#ifdef DEBUG
                std::cout << "Point " << pts->GetId(i) << " Polys " << polys->GetNumberOfIds() << std::endl;
#endif

                for (unsigned int j = 0; j < polys->GetNumberOfIds(); j++) {

                    pd->ReplaceCellPoint(polys->GetId(j), pts->GetId(i), pd->GetPoints()->InsertNextPoint(pt));

                }

            }

        }

    }

    pl->FreeSearchStructure();
    pl->Delete();

}

void vtkPolyDataBooleanFilter::MergePoints (vtkPolyData *pd, StripsType &strips) {

#ifdef DEBUG
    std::cout << "MergePoints()" << std::endl;
#endif

    pd->BuildLinks();

    vtkKdTreePointLocator *loc = vtkKdTreePointLocator::New();
    loc->SetDataSet(pd);
    loc->BuildLocator();

    std::map<unsigned int, std::vector<unsigned int> > ends;

    StripsType::const_iterator itr;

    for (itr = strips.begin(); itr != strips.end(); itr++) {

        vtkIdList *neigs = vtkIdList::New();

        loc->FindPointsWithinRadius(1e-7, (itr->begin()+1)->pt, neigs);

        for (unsigned int i = 0; i < neigs->GetNumberOfIds(); i++) {

            ends[itr->front().ind].push_back(neigs->GetId(i));

        }

        loc->FindPointsWithinRadius(1e-7, (itr->rbegin()+1)->pt, neigs);

        for (unsigned int i = 0; i < neigs->GetNumberOfIds(); i++) {

            ends[itr->back().ind].push_back(neigs->GetId(i));

        }

        neigs->Delete();

    }

    std::map<unsigned int, std::vector<unsigned int> >::const_iterator itr2;

    for (itr2 = ends.begin(); itr2 != ends.end(); itr2++) {

        double pt[3];

        contLines->GetPoint(itr2->first, pt);

        std::vector<MergeVertType> mergeVerts;

        vtkIdList *verts = vtkIdList::New();

        loc->FindPointsWithinRadius(1e-7, pt, verts);

#ifdef DEBUG
        std::cout << "ind " << itr2->first << ", neigs [";

        for (int i = 0; i < itr2->second.size(); i++) {
            std::cout << itr2->second[i] << ", ";
        }

        std::cout << "]" << std::endl;

#endif

        vtkIdList *polys = vtkIdList::New();

        for (unsigned int i = 0; i < verts->GetNumberOfIds(); i++) {
            pd->GetPointCells(verts->GetId(i), polys);

            if (polys->GetNumberOfIds() > 0) {
                // dieser vert wird von einem polygon verwendet

                vtkIdList *poly = vtkIdList::New();

                pd->GetCellPoints(polys->GetId(0), poly);

                // rel. index

                unsigned int ptNbr = poly->GetNumberOfIds();

                unsigned int j;

                for (j = 0; j < ptNbr; j++) {
                    if (verts->GetId(i) == poly->GetId(j)) {
                        break;
                    }
                }

                // es können beide zutreffen

                int indA = poly->GetId((j+1)%ptNbr);
                int indB = poly->GetId((j+ptNbr-1)%ptNbr);

                if (std::count(itr2->second.begin(), itr2->second.end(), indA) == 0) {

                    MergeVertType mergeVert;

                    mergeVert.polyInd = polys->GetId(0);
                    mergeVert.vert = poly->GetId(j);
                    mergeVert.vertInd = j;

                    pd->GetPoint(indA, mergeVert.pt);

#ifdef DEBUG
                    mergeVert.ind = indA;
#endif

                    mergeVerts.push_back(mergeVert);
                }

                if (std::count(itr2->second.begin(), itr2->second.end(), indB) == 0) {

                    MergeVertType mergeVert;

                    mergeVert.polyInd = polys->GetId(0);
                    mergeVert.vert = poly->GetId(j);
                    mergeVert.vertInd = j;

                    pd->GetPoint(indB, mergeVert.pt);

#ifdef DEBUG
                    mergeVert.ind = indB;
#endif

                    mergeVerts.push_back(mergeVert);
                }

                poly->Delete();
            }
        }

        polys->Delete();

        verts->Delete();

        // jetzt kommt die paarbildung

        std::vector<MergeVertType>::iterator itr3, itr4;

#ifdef DEBUG
        std::cout << "mergeVerts " << mergeVerts.size() << std::endl;

        for (itr3 = mergeVerts.begin(); itr3 != mergeVerts.end(); itr3++) {

            std::cout << "polyInd " << itr3->polyInd
                << ", vert " << itr3->vert
                << ", ind " << itr3->ind << std::endl;
        }

#endif

        if (mergeVerts.size() > 1) {

            std::deque<std::pair<MergeVertType, MergeVertType> > pairs;

            for (itr3 = mergeVerts.begin(); itr3 != mergeVerts.end()-1; itr3++) {

                if (itr3->used) {
                    continue;
                }

                for (itr4 = itr3+1; itr4 != mergeVerts.end(); itr4++) {

                    if (itr4->used) {
                        continue;
                    }

                    double dV[3], d;

                    vtkMath::Subtract(itr3->pt, itr4->pt, dV);

                    d = vtkMath::Norm(dV);

                    if (d < 1e-7) {
                        // paar gefunden

#ifdef DEBUG
                        std::cout << "pair [" << itr3->vert << ", " << itr4->vert << "]" << std::endl;
#endif

                        pairs.push_back(std::make_pair(*itr3, *itr4));

                        itr3->used = true;
                        itr4->used = true;

                        break;

                    }

                }

            }

            std::deque<MergeVertType> group;

            int i = 0;

            while (pairs.size() > 0) {

                if (group.empty()) {
                    group.push_back(pairs[0].first);
                    group.push_back(pairs[0].second);

                    pairs.pop_front();

                    i = 0;

                } else {
                    if (group.front().vert == pairs[i].first.vert) {
                        group.push_front(pairs[i].second);
                        pairs.erase(pairs.begin()+i);
                        i = 0;

                    } else if (group.front().vert == pairs[i].second.vert) {
                        group.push_front(pairs[i].first);
                        pairs.erase(pairs.begin()+i);
                        i = 0;

                    } else if (group.back().vert == pairs[i].first.vert) {
                        group.push_back(pairs[i].second);
                        pairs.erase(pairs.begin()+i);
                        i = 0;

                    } else if (group.back().vert == pairs[i].second.vert) {
                        group.push_back(pairs[i].first);
                        pairs.erase(pairs.begin()+i);
                        i = 0;

                    } else {
                        i++;
                    }
                }

                if (i == pairs.size()) {

                    std::deque<MergeVertType>::const_iterator itr4;

#ifdef DEBUG
                    std::cout << "group [" << group.front().vert << ", ";;
#endif

                    for (itr4 = group.begin()+1; itr4 != group.end(); itr4++) {
                        pd->ReplaceCellPoint(itr4->polyInd, itr4->vert, group.front().vert);

#ifdef DEBUG
                        std::cout << itr4->vert << ", ";
#endif

                    }

#ifdef DEBUG
                    std::cout << "]" << std::endl;
#endif

                    group.clear();

                }

            }

        }

    }

    loc->FreeSearchStructure();
    loc->Delete();

}


class PolyAtEdge {
    vtkPolyData *pd;
    vtkIdList *poly;

public:
    PolyAtEdge (vtkPolyData *pd, int polyId, int ptIdA, int ptIdB) : pd(pd), polyId(polyId), loc(LOC_NONE), poly(NULL) {
        ptIds[0] = ptIdA;
        ptIds[1] = ptIdB;

        poly = vtkIdList::New();
        pd->GetCellPoints(polyId, poly);

        for (unsigned int i = 0; i < poly->GetNumberOfIds(); i++) {
            if (poly->GetId(i) == ptIds[0]) {
                if (poly->GetId((i+1)%poly->GetNumberOfIds()) != ptIds[1]) {
                    int tmp = ptIds[0];
                    ptIds[0] = ptIds[1];
                    ptIds[1] = tmp;
                }

                break;
            }
        }

        double ptA[3], ptB[3];

        pd->GetPoint(ptIds[0], ptA);
        pd->GetPoint(ptIds[1], ptB);

        vtkMath::Subtract(ptB, ptA, e);
        vtkMath::Normalize(e);

        GeomHelper::ComputeNormal(pd->GetPoints(), poly, n);

        vtkMath::Cross(e, n, r);

    }

    ~PolyAtEdge () {
        poly->Delete();
    }

    int polyId, ptIds[2];
    double n[3], e[3], r[3];

    int loc;

};

class PolyPair {

public:
    double alpha;

    PolyPair (PolyAtEdge *pA, PolyAtEdge *pB) : pA(pA), pB(pB) {
        alpha = GeomHelper::GetAngle(pA->r, pB->r, pA->e);

    }

    ~PolyPair () {
        delete pA;
        delete pB;
    }

    PolyAtEdge *pA, *pB;

    void GetLoc (PolyAtEdge *pT, int mode) {
        double beta = GeomHelper::GetAngle(pA->r, pT->r, pA->e);

#ifdef DEBUG
        std::cout << "GetLoc()" << std::endl;

        std::cout << "polyId " << pT->polyId << ", beta " << (beta*180/M_PI) << std::endl;

#endif

        if (beta < 1e-7 || beta > 2*M_PI-1e-7) {
            // konkruent gegenüber dem polygon hinter pA

            double o = vtkMath::Dot(pA->n, pT->n);

            if (o < .999999) {
                // normalen sind entgegengesetzt gerichtet

                if (mode == OPER_INTERSECTION) {
                    pA->loc = LOC_OUTSIDE;
                    pT->loc = LOC_OUTSIDE;
                } else {
                    pA->loc = LOC_INSIDE;
                    pT->loc = LOC_INSIDE;
                }
            } else if (mode == OPER_UNION || mode == OPER_INTERSECTION) {
                pA->loc = LOC_INSIDE;
                pT->loc = LOC_OUTSIDE;
            }

        } else if (std::fabs(beta-alpha) < 1e-7) {

            double o = vtkMath::Dot(pB->n, pT->n);

            if (o < .999999) {
                // normalen sind entgegengesetzt gerichtet

                if (mode == OPER_INTERSECTION) {
                    pB->loc = LOC_OUTSIDE;
                    pT->loc = LOC_OUTSIDE;
                } else {
                    pB->loc = LOC_INSIDE;
                    pT->loc = LOC_INSIDE;
                }
            } else if (mode == OPER_UNION || mode == OPER_INTERSECTION) {
                pB->loc = LOC_INSIDE;
                pT->loc = LOC_OUTSIDE;
            }

        } else if (beta > alpha) {
            pT->loc = LOC_INSIDE;
        } else {
            pT->loc = LOC_OUTSIDE;
        }
    }

};

PolyPair* GetEdgePolys (vtkPolyData *pd, vtkIdList *ptsA, vtkIdList *ptsB) {

#ifdef DEBUG
    std::cout << "GetEdgePolys()" << std::endl;
#endif

    PolyPair *pp = NULL;

    std::map<int, int> polyPts;

    vtkIdList *polys = vtkIdList::New();

    for (unsigned int i = 0; i < ptsA->GetNumberOfIds(); i++) {
        polys->Reset();

        pd->GetPointCells(ptsA->GetId(i), polys);

        for (unsigned int j = 0; j < polys->GetNumberOfIds(); j++) {

            polyPts[polys->GetId(j)] = ptsA->GetId(i);

        }
    }

    PolyAtEdge *first = NULL;

    for (unsigned int i = 0; i < ptsB->GetNumberOfIds(); i++) {
        polys->Reset();

        pd->GetPointCells(ptsB->GetId(i), polys);

        for (unsigned int j = 0; j < polys->GetNumberOfIds(); j++) {

            if (polyPts.count(polys->GetId(j)) > 0) {

                if (first != NULL) {
                    PolyAtEdge *second = new PolyAtEdge(pd, polys->GetId(j), polyPts[polys->GetId(j)], ptsB->GetId(i));

                    pp = new PolyPair(first, second);

                    break;
                } else {
                    first = new PolyAtEdge(pd, polys->GetId(j), polyPts[polys->GetId(j)], ptsB->GetId(i));
                }

            }

        }

        if (pp != NULL) {
            break;
        }
    }

    polys->Delete();

    // eventuell tritt das hier nie ein
    if (pp == NULL && first != NULL) {
        delete first;
    }

#ifdef DEBUG

    if (pp != NULL) {
        std::cout << "pA: polyInd " << pp->pA->polyId
            << ", ptIdA " << pp->pA->ptIds[0] << ", ptIdB " << pp->pA->ptIds[1] << std::endl;

        std::cout << "pB: polyInd " << pp->pB->polyId
            << ", ptIdA " << pp->pB->ptIds[0] << ", ptIdB " << pp->pB->ptIds[1] << std::endl;

        std::cout << "alpha " << (pp->alpha*180/M_PI) << std::endl;

    }

#endif

    return pp;

}

void vtkPolyDataBooleanFilter::CombineRegions () {

#ifdef DEBUG
    std::cout << "CombineRegions()" << std::endl;
#endif

    // ungenutzte punkte löschen
    vtkCleanPolyData *cleanA = vtkCleanPolyData::New();
    cleanA->PointMergingOff();
#if (VTK_MAJOR_VERSION == 5)
    cleanA->SetInput(modPdA);
#else
    cleanA->SetInputData(modPdA);
#endif

    vtkCleanPolyData *cleanB = vtkCleanPolyData::New();
    cleanB->PointMergingOff();
#if (VTK_MAJOR_VERSION == 5)
    cleanB->SetInput(modPdB);
#else
    cleanB->SetInputData(modPdB);
#endif

    // regionen mit skalaren ausstatten
    vtkPolyDataConnectivityFilter *cfA = vtkPolyDataConnectivityFilter::New();
    cfA->SetExtractionModeToAllRegions();
    cfA->ColorRegionsOn();
    cfA->SetInputConnection(cleanA->GetOutputPort());

    vtkPolyDataConnectivityFilter *cfB = vtkPolyDataConnectivityFilter::New();
    cfB->SetExtractionModeToAllRegions();
    cfB->ColorRegionsOn();
    cfB->SetInputConnection(cleanB->GetOutputPort());

    cfA->Update();
    cfB->Update();

#ifdef DEBUG
    vtkPolyDataWriter *w = vtkPolyDataWriter::New();
    w->SetFileName("modPdA_5.vtk");
    w->SetInputConnection(cfA->GetOutputPort());
    std::cout << "Exporting modPdA_5.vtk" << std::endl;
    w->Update();

    w->SetFileName("modPdB_5.vtk");
    w->SetInputConnection(cfB->GetOutputPort());
    std::cout << "Exporting modPdB_5.vtk" << std::endl;
    w->Update();

    w->Delete();
#endif

    vtkPolyData *_pdA = cfA->GetOutput();
    vtkPolyData *_pdB = cfB->GetOutput();

    // locators erstellen
    vtkKdTreePointLocator *plA = vtkKdTreePointLocator::New();
    plA->SetDataSet(_pdA);
    plA->BuildLocator();

    vtkKdTreePointLocator *plB = vtkKdTreePointLocator::New();
    plB->SetDataSet(_pdB);
    plB->BuildLocator();

    _pdA->BuildLinks();
    _pdB->BuildLinks();

    vtkDataArray *scalarsA = _pdA->GetPointData()->GetScalars();
    vtkDataArray *scalarsB = _pdB->GetPointData()->GetScalars();

    std::map<int, int> locsA, locsB;

    for (unsigned int i = 0; i < contLines->GetNumberOfCells(); i++) {
        vtkIdList *line = vtkIdList::New();

        contLines->GetCellPoints(i, line);

        double ptA[3], ptB[3];

        contLines->GetPoint(line->GetId(0), ptA);
        contLines->GetPoint(line->GetId(1), ptB);

        vtkIdList *fptsA = vtkIdList::New();
        vtkIdList *lptsA = vtkIdList::New();

        plA->FindPointsWithinRadius(1e-7, ptA, fptsA);
        plA->FindPointsWithinRadius(1e-7, ptB, lptsA);

        vtkIdList *fptsB = vtkIdList::New();
        vtkIdList *lptsB = vtkIdList::New();

        plB->FindPointsWithinRadius(1e-7, ptA, fptsB);
        plB->FindPointsWithinRadius(1e-7, ptB, lptsB);

#ifdef DEBUG
        std::cout << "line " << i << std::endl;
#endif

        PolyPair *ppA = GetEdgePolys(_pdA, fptsA, lptsA);
        PolyPair *ppB = GetEdgePolys(_pdB, fptsB, lptsB);

        if (ppA != NULL && ppB != NULL) {

            ppB->GetLoc(ppA->pA, OperMode);
            ppB->GetLoc(ppA->pB, OperMode);

            ppA->GetLoc(ppB->pA, OperMode);
            ppA->GetLoc(ppB->pB, OperMode);

            int fsA = scalarsA->GetTuple1(ppA->pA->ptIds[0]);
            int lsA = scalarsA->GetTuple1(ppA->pB->ptIds[0]);

            int fsB = scalarsB->GetTuple1(ppB->pA->ptIds[0]);
            int lsB = scalarsB->GetTuple1(ppB->pB->ptIds[0]);

#ifdef DEBUG
            std::cout << "polyId " << ppA->pA->polyId << ", sA " << fsA << ", loc " << ppA->pA->loc << std::endl;
            std::cout << "polyId " << ppA->pB->polyId << ", sA " << lsA << ", loc " << ppA->pB->loc << std::endl;
            std::cout << "polyId " << ppB->pA->polyId << ", sB " << fsB << ", loc " << ppB->pA->loc << std::endl;
            std::cout << "polyId " << ppB->pB->polyId << ", sB " << lsB << ", loc " << ppB->pB->loc << std::endl;
#endif

            locsA[fsA] = ppA->pA->loc;
            locsA[lsA] = ppA->pB->loc;

            locsB[fsB] = ppB->pA->loc;
            locsB[lsB] = ppB->pB->loc;

        }

        if (ppA != NULL) {
            delete ppA;
        }

        if (ppB != NULL) {
            delete ppB;
        }

        lptsB->Delete();
        fptsB->Delete();

        lptsA->Delete();
        fptsA->Delete();

        line->Delete();

    }

    // reale kombination der ermittelten regionen

    int comb[] = {LOC_OUTSIDE, LOC_OUTSIDE};

    if (OperMode == OPER_INTERSECTION) {
        comb[0] = LOC_INSIDE;
        comb[1] = LOC_INSIDE;
    } else if (OperMode == OPER_DIFFERENCE) {
        comb[1] = LOC_INSIDE;
    } else if (OperMode == OPER_DIFFERENCE2) {
        comb[0] = LOC_INSIDE;
    }

    cfA->SetExtractionModeToSpecifiedRegions();
    cfB->SetExtractionModeToSpecifiedRegions();

    std::map<int, int>::const_iterator itr;

    for (itr = locsA.begin(); itr != locsA.end(); itr++) {
        if (itr->second == comb[0]) {
            cfA->AddSpecifiedRegion(itr->first);
        }
    }

    for (itr = locsB.begin(); itr != locsB.end(); itr++) {
        if (itr->second == comb[1]) {
            cfB->AddSpecifiedRegion(itr->first);
        }
    }

    // nach innen zeigende normalen umkehren

    cfA->Update();

    vtkPolyData *regsA = cfA->GetOutput();

    if (comb[0] == LOC_INSIDE && OperMode != OPER_INTERSECTION) {
        for (unsigned int i = 0; i < regsA->GetNumberOfCells(); i++) {
            regsA->ReverseCell(i);
        }
    }

    cfB->Update();

    vtkPolyData *regsB = cfB->GetOutput();

    if (comb[1] == LOC_INSIDE && OperMode != OPER_INTERSECTION) {
        for (unsigned int i = 0; i < regsB->GetNumberOfCells(); i++) {
            regsB->ReverseCell(i);
        }
    }

    vtkAppendPolyData *app = vtkAppendPolyData::New();

#if (VTK_MAJOR_VERSION == 5)
    app->AddInput(regsA);
    app->AddInput(regsB);
#else
    app->AddInputData(regsA);
    app->AddInputData(regsB);
#endif

    // entfernt ungenutzte punkte
    vtkCleanPolyData *cleanApp = vtkCleanPolyData::New();
    cleanApp->PointMergingOff();
    cleanApp->SetInputConnection(app->GetOutputPort());

    // färbt die regionen nochmal neu ein, damit mehrere regionen nicht die gleiche farbe haben

    vtkPolyDataConnectivityFilter *cfApp = vtkPolyDataConnectivityFilter::New();
    cfApp->SetExtractionModeToAllRegions();
    cfApp->ColorRegionsOn();
    cfApp->SetInputConnection(cleanApp->GetOutputPort());

    cfApp->Update();

    // resultA ist erster output des filters
    resultA->DeepCopy(cfApp->GetOutput());

    // aufräumen

    cfApp->Delete();
    cleanApp->Delete();
    app->Delete();

    plB->FreeSearchStructure();
    plB->Delete();

    plA->FreeSearchStructure();
    plA->Delete();

    cfB->Delete();
    cfA->Delete();

    cleanB->Delete();
    cleanA->Delete();

}
