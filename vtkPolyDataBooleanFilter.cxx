/*
   Copyright 2012-2016 Ronald Römer

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
#include <queue>
#include <array>

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkObjectFactory.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkMath.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkAppendPolyData.h>
#include <vtkKdTreePointLocator.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>

#include "vtkPolyDataBooleanFilter.h"
#include "vtkPolyDataContactFilter.h"

#include "Utilities.h"
using Utilities::IdsType;
using Utilities::GetNormal;

#include "Decomposer.h"

using Decomposer::Base;
using Decomposer::Transform;
using Decomposer::DecType;
using Decomposer::Decompose;

#ifdef DEBUG
#include <ctime>
#define DIFF(s, e) float(e-s)/CLOCKS_PER_SEC
#include <numeric>
#include <iterator>
#endif

#ifdef _DEBUG
#include <vtkPolyDataWriter.h>
#include <vtkExtractCells.h>
#include <vtkGeometryFilter.h>
#endif

vtkStandardNewMacro(vtkPolyDataBooleanFilter);

vtkPolyDataBooleanFilter::vtkPolyDataBooleanFilter () {

    SetNumberOfInputPorts(2);
    SetNumberOfOutputPorts(2);

    timePdA = 0;
    timePdB = 0;

    contLines = vtkPolyData::New();

    modPdA = vtkPolyData::New();
    modPdB = vtkPolyData::New();

    cellDataA = vtkCellData::New();
    cellDataB = vtkCellData::New();

    OperMode = OPER_UNION;

    MergeAll = false;
    DecPolys = false;

}

vtkPolyDataBooleanFilter::~vtkPolyDataBooleanFilter () {

    cellDataB->Delete();
    cellDataA->Delete();

    modPdB->Delete();
    modPdA->Delete();

    contLines->Delete();

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
        resultB = vtkPolyData::SafeDownCast(outInfoB->Get(vtkDataObject::DATA_OBJECT()));

#ifdef DEBUG
        std::vector<float> times;
        clock_t start;
#endif

        if (pdA->GetMTime() > timePdA || pdB->GetMTime() > timePdB) {

            // eventuell vorhandene regionen vereinen

            vtkCleanPolyData *cleanA = vtkCleanPolyData::New();
            cleanA->SetTolerance(1e-9);
            cleanA->SetInputData(pdA);
            cleanA->Update();

            vtkCleanPolyData *cleanB = vtkCleanPolyData::New();
            cleanB->SetTolerance(1e-9);
            cleanB->SetInputData(pdB);
            cleanB->Update();

#ifdef DEBUG
            std::cout << "Exporting modPdA.vtk" << std::endl;
            Utilities::WriteVTK("modPdA.vtk", cleanA->GetOutput());

            std::cout << "Exporting modPdB.vtk" << std::endl;
            Utilities::WriteVTK("modPdB.vtk", cleanB->GetOutput());
#endif

            // CellData sichern

            cellDataA->DeepCopy(cleanA->GetOutput()->GetCellData());
            cellDataB->DeepCopy(cleanB->GetOutput()->GetCellData());

            // ermittelt kontaktstellen

#ifdef DEBUG

            start = std::clock();
#endif

            vtkPolyDataContactFilter *cl = vtkPolyDataContactFilter::New();
            cl->SetInputConnection(0, cleanA->GetOutputPort());
            cl->SetInputConnection(1, cleanB->GetOutputPort());
            cl->Update();

#ifdef DEBUG
            times.push_back(DIFF(start, std::clock()));
#endif

            contLines->DeepCopy(cl->GetOutput());

#ifdef DEBUG
            std::cout << "Exporting contLines.vtk" << std::endl;
            Utilities::WriteVTK("contLines.vtk", contLines);

            std::cout << "Exporting modPdA_1.vtk" << std::endl;
            Utilities::WriteVTK("modPdA_1.vtk", cl->GetOutput(1));

            std::cout << "Exporting modPdB_1.vtk" << std::endl;
            Utilities::WriteVTK("modPdB_1.vtk", cl->GetOutput(2));
#endif

            modPdA->DeepCopy(cl->GetOutput(1));
            modPdB->DeepCopy(cl->GetOutput(2));

            if (contLines->GetNumberOfCells() == 0) {
                vtkErrorMacro("Inputs have no contact.");

                return 1;
            }

            // in den CellDatas steht drin, welche polygone einander schneiden

            vtkIntArray *contsA = vtkIntArray::SafeDownCast(contLines->GetCellData()->GetScalars("cA"));
            vtkIntArray *contsB = vtkIntArray::SafeDownCast(contLines->GetCellData()->GetScalars("cB"));

            // gültigkeit des schnitts prüfen

            bool valid = true;

            vtkIdList *cells = vtkIdList::New();

            for (int i = 0; i < contLines->GetNumberOfPoints() && valid; i++) {
                cells->Reset();
                contLines->GetPointCells(i, cells);

                // der schnitt endet abrupt
                if (cells->GetNumberOfIds() == 1) {
                    valid = false;
                    break;
                }

                vtkIdList *line = vtkIdList::New();

                std::map<int, std::set<int> > countsA, countsB;

                for (int j = 0; j < cells->GetNumberOfIds(); j++) {
                    contLines->GetCellPoints(cells->GetId(j), line);

                    int target = line->GetId(0) == i ? line->GetId(1) : line->GetId(0);

                    int indA = contsA->GetValue(cells->GetId(j));
                    int indB = contsB->GetValue(cells->GetId(j));

                    // sind am punkt sind mehr als zwei linien beteiligt?

                    countsA[indA].insert(target);

                    if (countsA[indA].size() > 2) {
                        valid = false;
                        break;
                    }

                    countsB[indB].insert(target);

                    if (countsB[indB].size() > 2) {
                        valid = false;
                        break;
                    }

                }

                line->Delete();
            }

            cells->Delete();

            if (!valid) {
                vtkErrorMacro("Contact is ambiguous or incomplete.");

                return 1;
            }

#ifdef DEBUG
            start = std::clock();
#endif

            GetPolyStrips(modPdA, contsA, polyStripsA);
            GetPolyStrips(modPdB, contsB, polyStripsB);

#ifdef DEBUG
            times.push_back(DIFF(start, std::clock()));
#endif

            // trennt die polygone an den linien

#ifdef DEBUG
            start = std::clock();
#endif

            CutCells(modPdA, polyStripsA);
            CutCells(modPdB, polyStripsB);

#ifdef DEBUG
            times.push_back(DIFF(start, std::clock()));
#endif

#ifdef DEBUG
            std::cout << "Exporting modPdA_2.vtk" << std::endl;
            Utilities::WriteVTK("modPdA_2.vtk", modPdA);

            std::cout << "Exporting modPdB_2.vtk" << std::endl;
            Utilities::WriteVTK("modPdB_2.vtk", modPdB);
#endif

#ifdef DEBUG
            start = std::clock();
#endif

            RestoreOrigPoints(modPdA, polyStripsA);
            RestoreOrigPoints(modPdB, polyStripsB);

#ifdef DEBUG
            times.push_back(DIFF(start, std::clock()));
#endif

#ifdef DEBUG
            std::cout << "Exporting modPdA_3.vtk" << std::endl;
            Utilities::WriteVTK("modPdA_3.vtk", modPdA);

            std::cout << "Exporting modPdB_3.vtk" << std::endl;
            Utilities::WriteVTK("modPdB_3.vtk", modPdB);
#endif

#ifdef DEBUG
            start = std::clock();
#endif

            ResolveOverlaps(modPdA, contsA, polyStripsA);
            ResolveOverlaps(modPdB, contsB, polyStripsB);

#ifdef DEBUG
            times.push_back(DIFF(start, std::clock()));
#endif

#ifdef DEBUG
            std::cout << "Exporting modPdA_4.vtk" << std::endl;
            Utilities::WriteVTK("modPdA_4.vtk", modPdA);

            std::cout << "Exporting modPdB_4.vtk" << std::endl;
            Utilities::WriteVTK("modPdB_4.vtk", modPdB);
#endif

#ifdef DEBUG
            start = std::clock();
#endif

            AddAdjacentPoints(modPdA, contsA, polyStripsA);
            AddAdjacentPoints(modPdB, contsB, polyStripsB);

#ifdef DEBUG
            times.push_back(DIFF(start, std::clock()));
#endif

#ifdef DEBUG
            std::cout << "Exporting modPdA_5.vtk" << std::endl;
            Utilities::WriteVTK("modPdA_5.vtk", modPdA);

            std::cout << "Exporting modPdB_5.vtk" << std::endl;
            Utilities::WriteVTK("modPdB_5.vtk", modPdB);
#endif

#ifdef DEBUG
            start = std::clock();
#endif

            DisjoinPolys(modPdA, polyStripsA);
            DisjoinPolys(modPdB, polyStripsB);

#ifdef DEBUG
            times.push_back(DIFF(start, std::clock()));
#endif

#ifdef DEBUG
            std::cout << "Exporting modPdA_6.vtk" << std::endl;
            Utilities::WriteVTK("modPdA_6.vtk", modPdA);

            std::cout << "Exporting modPdB_6.vtk" << std::endl;
            Utilities::WriteVTK("modPdB_6.vtk", modPdB);
#endif

#ifdef DEBUG
            start = std::clock();
#endif

            MergePoints(modPdA, polyStripsA);
            MergePoints(modPdB, polyStripsB);

#ifdef DEBUG
            times.push_back(DIFF(start, std::clock()));
#endif

#ifdef DEBUG
            std::cout << "Exporting modPdA_7.vtk" << std::endl;
            Utilities::WriteVTK("modPdA_7.vtk", modPdA);

            std::cout << "Exporting modPdB_7.vtk" << std::endl;
            Utilities::WriteVTK("modPdB_7.vtk", modPdB);
#endif

            involvedA.clear();
            involvedB.clear();

            int numLines = contLines->GetNumberOfCells();
            for (int i = 0; i < numLines; i++) {
                involvedA.insert(contsA->GetValue(i));
                involvedB.insert(contsB->GetValue(i));
            }

            // aufräumen

            cl->Delete();
            cleanB->Delete();
            cleanA->Delete();

            timePdA = pdA->GetMTime();
            timePdB = pdB->GetMTime();

        }

#ifdef DEBUG
        start = std::clock();
#endif

        if (MergeAll) {
            MergeRegions();
        } else {
            CombineRegions();
        }

#ifdef DEBUG
        times.push_back(DIFF(start, std::clock()));
#endif

        if (DecPolys) {

#ifdef DEBUG
        start = std::clock();
#endif

            DecomposePolys();

#ifdef DEBUG
        times.push_back(DIFF(start, std::clock()));
#endif
        }

#ifdef DEBUG
        // http://www.gnu.org/software/libc/manual/html_node/CPU-Time.html

        float sum = std::accumulate(times.begin(), times.end(), 0.);

        std::vector<float>::const_iterator itr;
        for (itr = times.cbegin(); itr != times.cend(); itr++) {
            std::cout << "Time " << std::distance(times.cbegin(), itr)
                << ": " << *itr << " (" << (*itr/sum)*100 << "%)"
                << std::endl;
        }
#endif

    }

    return 1;

}

void vtkPolyDataBooleanFilter::GetStripPoints (vtkPolyData *pd, PStrips &pStrips, IdsType &lines) {

#ifdef DEBUG
    std::cout << "GetStripPoints()" << std::endl;
#endif

    StripPtsType &pts = pStrips.pts;
    IdsType &poly = pStrips.poly;

    int numPts = poly.size();

    double tol = 1e-5;

    IdsType::const_iterator itr;

    for (itr = lines.begin(); itr != lines.end(); itr++) {

        vtkIdList *linePts = vtkIdList::New();
        contLines->GetCellPoints(*itr, linePts);

        // diese punkte durchlaufen

        for (int i = 0; i < 2; i++) {

            int realInd = linePts->GetId(i);

            if (pts.count(realInd) == 0) {
                // lage analysieren

                pts[realInd].ind = realInd;

                // die koordinaten
                double pt[3];
                contLines->GetPoint(realInd, pt);

                CPY(pts[realInd].pt, pt)

                double lastD = DBL_MAX;

                // jetzt muss man die kanten durchlaufen
                for (int j = 0; j < numPts; j++) {

                    int indA, indB;
                    double a[3], b[3];

                    indA = poly[j];
                    indB = poly[(j+1)%numPts];

                    pd->GetPoint(indA, a);
                    pd->GetPoint(indB, b);

                    double sA[3], sB[3];
                    vtkMath::Subtract(a, pt, sA);
                    vtkMath::Subtract(b, pt, sB);

                    // richtungsvektor und länge der kante

                    double u[3];
                    vtkMath::Subtract(b, a, u);

                    double n = vtkMath::Norm(u);

                    // d und t zur kante

                    double v[3];
                    vtkMath::Subtract(pt, a, v);

                    double t = vtkMath::Dot(v, u)/(n*n);

                    double w[3];
                    vtkMath::Cross(v, u, w);

                    double d = vtkMath::Norm(w)/n;

                    if (d < tol && d < lastD && t > -tol && t < 1+tol) {
                        // liegt im toleranzbereich der kante

                        pts[realInd].edge[0] = indA;
                        pts[realInd].edge[1] = indB;

                        pts[realInd].t = std::min(1., std::max(0., t));

                        if (vtkMath::Norm(sA) < tol) {
                            CPY(pts[realInd].captPt, a)
                            pts[realInd].capt = CAPT_A;

                            break;

                        } else if (vtkMath::Norm(sB) < tol) {
                            CPY(pts[realInd].captPt, b)
                            pts[realInd].capt = CAPT_B;

                            break;

                        } else {
                            // u ist nicht normiert
                            vtkMath::MultiplyScalar(u, t);

                            double x[3];
                            vtkMath::Add(a, u, x);

                            // projektion
                            CPY(pts[realInd].captPt, x)

                            pts[realInd].capt = CAPT_EDGE;

                            lastD = d;
                        }

                    }
                }
            }
        }

        linePts->Delete();

    }

    StripPtsType::iterator itr2;

    for (itr2 = pts.begin(); itr2 != pts.end(); ++itr2) {
        StripPt &sp = itr2->second;

        if (sp.capt != CAPT_NOT) {
            if (sp.capt == CAPT_B) {
                sp.t = 0;

                sp.edge[0] = sp.edge[1];

                int i = std::find(poly.begin(), poly.end(), sp.edge[0])-poly.begin();

                sp.edge[1] = poly[(i+1)%numPts];

                sp.capt = CAPT_A;

            }

            // für den schnitt werden die eingerasteten koordinaten verwendet

            CPY(sp.cutPt, sp.captPt)
        } else {

            CPY(sp.cutPt, sp.pt)
        }

    }

#ifdef DEBUG
    std::cout << "pts: " << std::endl;
    for (itr2 = pts.begin(); itr2 != pts.end(); ++itr2) {
        std::cout << itr2->first << ": " << itr2->second << std::endl;
    }
#endif

}

void vtkPolyDataBooleanFilter::GetPolyStrips (vtkPolyData *pd, vtkIntArray *conts, PolyStripsType &polyStrips) {
#ifdef DEBUG
    std::cout << "GetPolyStrips()" << std::endl;
#endif

    polyStrips.clear();

    std::map<int, IdsType> polyLines;

    for (int i = 0; i < conts->GetNumberOfTuples(); i++) {
        int poly = conts->GetValue(i);
        polyLines[poly].push_back(i);
    }

    std::map<int, IdsType>::iterator itr;

    for (itr = polyLines.begin(); itr != polyLines.end(); ++itr) {
        IdsType &lines = itr->second;
        RemoveDuplicates(lines);

        PStrips &pStrips = polyStrips[itr->first];

        vtkIdList *polyPts = vtkIdList::New();
        pd->GetCellPoints(itr->first, polyPts);

        int numPts = polyPts->GetNumberOfIds();
		std::vector<std::array<double, 3> > pts_(numPts);

        //double pts_[numPts][3];	

        for (int i = 0; i < numPts; i++) {
            pStrips.poly.push_back(polyPts->GetId(i));
            pd->GetPoint(pStrips.poly.back(), pts_[i].data());
        }

        GetNormal(pts_, pStrips.n, numPts);

        GetStripPoints(pd, pStrips, lines);

        StripsType &strips = pStrips.strips;
        StripType strip;

        // zusammensetzen

        std::deque<int> _lines(lines.begin(), lines.end());

        int i = 0;

        while (_lines.size() > 0) {
            vtkIdList *linePts = vtkIdList::New();
            contLines->GetCellPoints(_lines[i], linePts);

            int indA = linePts->GetId(0);
            int indB = linePts->GetId(1);

            if (strip.empty()) {
                strip.push_back(StripPtR(indA));
                strip.push_back(StripPtR(indB));

                _lines.erase(_lines.begin());
            } else {
                StripPt &start = pStrips.pts[strip.front().ind],
                    &end = pStrips.pts[strip.back().ind];

                if (end.capt == CAPT_NOT && end.ind == indA) {
                    strip.push_back(StripPtR(indB));
                    _lines.erase(_lines.begin()+i);
                    i = 0;

                } else if (end.capt == CAPT_NOT && end.ind == indB) {
                    strip.push_back(StripPtR(indA));
                    _lines.erase(_lines.begin()+i);
                    i = 0;

                } else if (start.capt == CAPT_NOT && start.ind == indA) {
                    strip.push_front(StripPtR(indB));
                    _lines.erase(_lines.begin()+i);
                    i = 0;

                } else if (start.capt == CAPT_NOT && start.ind == indB) {
                    strip.push_front(StripPtR(indA));
                    _lines.erase(_lines.begin()+i);
                    i = 0;

                } else {
                    i++;

                }
            }

            StripPt &_start = pStrips.pts[strip.front().ind],
                &_end = pStrips.pts[strip.back().ind];

            if ((_start.capt != CAPT_NOT && _end.capt != CAPT_NOT)
                || (_start.ind == _end.ind)
                || (i == _lines.size())) {

                // einen neuen strip anlegen
                strips.push_back(strip);
                strip.clear();

                i = 0;
            }

            linePts->Delete();

        }

        CompleteStrips(pStrips);

    }

}

void vtkPolyDataBooleanFilter::RemoveDuplicates (IdsType &lines) {

    IdsType unique;

    int i, j;

    // die indexe der enden auf übereinstimmung prüfen

    vtkIdList *linePtsA = vtkIdList::New();
    vtkIdList *linePtsB = vtkIdList::New();

    int numLines = lines.size();

    for (i = 0; i < numLines-1; i++) {
        j = i+1;

        contLines->GetCellPoints(lines[i], linePtsA);

        while (j < lines.size()) {
            contLines->GetCellPoints(lines[j], linePtsB);

            if ((linePtsA->GetId(0) == linePtsB->GetId(0) && linePtsA->GetId(1) == linePtsB->GetId(1)) ||
                (linePtsA->GetId(0) == linePtsB->GetId(1) && linePtsA->GetId(1) == linePtsB->GetId(0))) {
                // stimmen überein
                break;
            }

            j++;
        }

        if (j == numLines) {
            // keine vorzeitige unterbrechung der while-schleife
            unique.push_back(lines[i]);
        }
    }

    unique.push_back(lines.back());

    linePtsA->Delete();
    linePtsB->Delete();

    lines.swap(unique);

}

void vtkPolyDataBooleanFilter::CompleteStrips (PStrips &pStrips) {
    StripsType::iterator itr;

    for (itr = pStrips.strips.begin(); itr != pStrips.strips.end(); ++itr) {
        StripPt &start = pStrips.pts[itr->front().ind],
            &end = pStrips.pts[itr->back().ind];

        if (start.ind != end.ind) {
            if (start.capt == CAPT_NOT) {
                itr->insert(itr->begin(), itr->rbegin(), itr->rend()-1);
            } else if (end.capt == CAPT_NOT) {
                itr->insert(itr->end(), itr->rbegin()+1, itr->rend());
            }
        }
    }
}

bool vtkPolyDataBooleanFilter::HasArea (StripType &strip) {
    bool area = true;

    int n = strip.size();
    if (n%2 == 1) {
        for (int i = 0; i < (n-1)/2; i++) {
            area = strip[i].ind != strip[n-i-1].ind;
        }
    }

    return area;
}

void vtkPolyDataBooleanFilter::CutCells (vtkPolyData *pd, PolyStripsType &polyStrips) {
    vtkPoints *pdPts = pd->GetPoints();

    vtkIntArray *origCellIds = vtkIntArray::SafeDownCast(pd->GetCellData()->GetScalars("OrigCellIds"));

    PolyStripsType::iterator itr;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {

        int polyInd = itr->first;
        PStrips &pStrips = itr->second;

#ifdef DEBUG
        if (polyInd != 204) {
            //continue;
        }
#endif

        StripsType &strips = pStrips.strips;
        StripPtsType &pts = pStrips.pts;

        IdsType &poly = pStrips.poly;

#ifdef DEBUG
        IdsType::iterator itr_;
        std::cout << "polyInd=" << polyInd << ", poly=[";
        for (itr_ = poly.begin(); itr_ != poly.end(); ++itr_) {
            std::cout << *itr_ << ", ";
        }
        std::cout << "]" << std::endl;
#endif

        int numPts = poly.size();

        std::map<int, RefsType> edges;

        StripsType::iterator itr2;
        StripType::iterator itr3;

        // holes löschen
        strips.erase(std::remove_if(strips.begin(), strips.end(), [&](const StripType &s) {
            return pts[s.front().ind].capt == CAPT_NOT && pts[s.back().ind].capt == CAPT_NOT; }), strips.end());

        for (itr2 = strips.begin(); itr2 != strips.end(); ++itr2) {
            StripType &strip = *itr2;

            // init
            if (pts[strip.front().ind].edge[0] == pts[strip.back().ind].edge[0]
                && strip.front().ind != strip.back().ind
                && pts[strip.front().ind].t > pts[strip.back().ind].t) {

                std::reverse(strip.begin(), strip.end());
            }

            StripPt &start = pts[strip.front().ind],
                &end = pts[strip.back().ind];

            strip.front().side = SIDE_START;
            strip.back().side = SIDE_END;

            strip.front().ref = start.edge[0];
            strip.back().ref = end.edge[0];

            int ind = itr2-strips.begin();

            strip.front().strip = strip.back().strip = ind;

            for (itr3 = strip.begin(); itr3 != strip.end(); ++itr3) {
                StripPt &sp = pts[itr3->ind];

                itr3->desc[0] = pdPts->InsertNextPoint(sp.cutPt);
                itr3->desc[1] = pdPts->InsertNextPoint(sp.cutPt);

#ifdef DEBUG
                std::cout << sp << " => " << *itr3 << std::endl;
#endif

            }

            // ordnet zu

            edges[start.edge[0]].push_back(std::ref(strip.front()));
            edges[end.edge[0]].push_back(std::ref(strip.back()));
        }

        // sortiert die punkte auf den kanten

        std::map<int, RefsType>::iterator itr4;

        IdsType::iterator itr5;
        StripType::reverse_iterator itr6;

        RefsType::iterator itr7;

        std::vector<IdsType>::iterator itr8;


        for (itr4 = edges.begin(); itr4 != edges.end(); ++itr4) {
            RefsType &edge = itr4->second;

#ifdef DEBUG
            std::cout << "edge (" << itr4->first << ", _)" << std::endl;
#endif

            std::sort(edge.begin(), edge.end(), [&](const StripPtR &a, const StripPtR &b) {
                StripPt &a_ = pts[a.ind],
                    &b_ = pts[b.ind];

#ifdef DEBUG
                std::cout << "a_: " << a_ << " -> strip " << a.strip << std::endl;
                std::cout << "b_: " << b_ << " -> strip " << b.strip << std::endl;
#endif

                if (a_.ind == b_.ind) {
                    // strips beginnen im gleichen punkt

                    if (a.strip != b.strip) {
                        // gehören nicht dem gleichen strip an

                        StripType &stripA = strips[a.strip],
                            &stripB = strips[b.strip];

                        // andere enden ermitteln

                        StripPtR &eA = (&a == &(stripA.front())) ? stripA.back() : stripA.front(),
                            &eB = (&b == &(stripB.front())) ? stripB.back() : stripB.front();

                        StripPt &eA_ = pts[eA.ind],
                            &eB_ = pts[eB.ind];

#ifdef DEBUG
                        std::cout << "eA_: " << eA_ << std::endl;
                        std::cout << "eB_: " << eA_ << std::endl;
#endif

                        if (eA_.ind != eB_.ind) {
                            int i = std::find(poly.begin(), poly.end(), itr4->first)-poly.begin();

                            int iA = std::find(poly.begin(), poly.end(), eA_.edge[0])-poly.begin(),
                                iB = std::find(poly.begin(), poly.end(), eB_.edge[0])-poly.begin();

                            double dA = Utilities::Mod(iA-i, numPts)+eA_.t,
                                dB = Utilities::Mod(iB-i, numPts)+eB_.t;

                            if (i == iA && a_.t > eA_.t) {
                               dA += numPts;
                            }

                            if (i == iB && b_.t > eB_.t) {
                               dB += numPts;
                            }

#ifdef DEBUG
                            std::cout << "dA=" << dA << ", dB=" << dB << std::endl;
#endif

                            return dB < dA;
                        } else {
                            RefsType poly_;

                            if (a.side == SIDE_START) {
                                poly_.insert(poly_.end(), stripA.begin(), stripA.end());
                            } else {
                                poly_.insert(poly_.end(), stripA.rbegin(), stripA.rend());
                            }

                            if (b.side == SIDE_START) {
                                poly_.insert(poly_.end(), stripB.rbegin()+1, stripB.rend()-1);
                            } else {
                                poly_.insert(poly_.end(), stripB.begin()+1, stripB.end()-1);
                            }

                            int num = poly_.size();
                            //double pts_[num][3];
							std::vector<std::array<double, 3> > pts_(num);
                            for (itr7 = poly_.begin(); itr7 != poly_.end(); ++itr7) {
                                StripPtR &sp = *itr7;
                                int i = itr7-poly_.begin();

                                CPY(pts_[i], pts[sp.ind].cutPt)
                            }

                            double n[3];
                            GetNormal(pts_, n, num);

                            double ang = vtkMath::Dot(pStrips.n, n);

#ifdef DEBUG
                            std::cout << "ang=" << ang*180/pi << std::endl;
#endif

                            return ang < .999999;

                        }
                    } else {
                        // gleicher strip

                        StripType &strip = strips[a.strip];

                        if (HasArea(strip)) {
                            RefsType poly_(strip.begin(), strip.end()-1);

                            int num = poly_.size();
                            //double pts_[num][3];
							std::vector<std::array<double, 3> > pts_(num);

                            for (itr7 = poly_.begin(); itr7 != poly_.end(); ++itr7) {
                                StripPtR &sp = *itr7;
                                int i = itr7-poly_.begin();

                                CPY(pts_[i], pts[sp.ind].cutPt)
                            }

                            double n[3];
                            GetNormal(pts_, n, num);

                            double ang = vtkMath::Dot(pStrips.n, n);

                            if (ang > .999999) {
                                std::reverse(strip.begin(), strip.end());
                                return true;
                            } else {
                                return false;
                            }

                        } else {
                            // reihenfolge von a und b bereits richtig
                            return false;
                        }
                    }

                } else {
                    return a_.t < b_.t;
                }
            });

#ifdef DEBUG
            std::cout << "after sort" << std::endl;
            for (itr7 = edge.begin(); itr7 != edge.end(); ++itr7) {
                std::cout << itr7->get().strip << ": " << pts[itr7->get().ind] << std::endl;
            }
#endif

        }

        // baut die strips ein

        std::deque<IdsType> polys;
        polys.push_back(pStrips.poly);

        int i = 0;

        std::deque<IdsType>::iterator itr9;

        for (itr2 = strips.begin(); itr2 != strips.end(); ++itr2) {
            StripType &strip = *itr2;

            StripPtR &start = strip.front(),
                &end = strip.back();

#ifdef DEBUG
            std::cout << "strip " << start.strip
                << " refs=[" << start.ref << ", " << end.ref << "]"
                << std::endl;
#endif

            int cycle = 0;

            while (true) {

                if (cycle == polys.size()) {
                    break;
                }

                IdsType p = polys.front();
                polys.pop_front();

                std::vector<IdsType> divided(2);

                if (std::find(p.begin(), p.end(), start.ref) != p.end()) {
                    if (start.ref == end.ref) {
                        for (itr5 = p.begin(); itr5 != p.end(); ++itr5) {
                            divided[0].push_back(*itr5);

#ifdef DEBUG
                            std::cout << "add d[0] " << *itr5 << std::endl;
#endif

                            if (*itr5 == start.ref) {
                                for (itr3 = strip.begin(); itr3 != strip.end(); ++itr3) {
                                    divided[0].push_back(itr3->desc[0]);

#ifdef DEBUG
                                    std::cout << "add d[0] " << itr3->desc[0] << std::endl;
#endif

                                }
                            }
                        }

                        // strip selbst ist ein polygon

                        for (itr6 = strip.rbegin(); itr6 != strip.rend(); ++itr6) {
                            divided[1].push_back(itr6->desc[1]);

#ifdef DEBUG
                            std::cout << "add d[1] " << itr6->desc[1] << std::endl;
#endif

                        }

                    } else {
                        int d = 0;
                        for (itr5 = p.begin(); itr5 != p.end(); ++itr5) {
                            divided[d].push_back(*itr5);

#ifdef DEBUG
                            std::cout << "add d[" << d << "] " << *itr5 << std::endl;
#endif

                            if (*itr5 == start.ref) {
                                for (itr3 = strip.begin(); itr3 != strip.end(); ++itr3) {
                                    divided[d].push_back(itr3->desc[0]);

#ifdef DEBUG
                                    std::cout << "add d[" << d << "] " << itr3->desc[0] << std::endl;
#endif

                                }
                                d++;
                            } else if (*itr5 == end.ref) {
                                for (itr6 = strip.rbegin(); itr6 != strip.rend(); ++itr6) {
                                    divided[d].push_back(itr6->desc[1]);

#ifdef DEBUG
                                    std::cout << "add d[" << d << "] " << itr6->desc[1] << std::endl;
#endif

                                }
                                d++;
                            }

                            d %= 2;
                        }
                    }
                }

                if (divided[1].size() > 0) {

                    // refs aktualisieren

                    int stripInd = itr2-strips.begin();

                    for (itr4 = edges.begin(); itr4 != edges.end(); ++itr4) {
                        RefsType &edge = itr4->second;

                        for (itr7 = edge.begin()+1; itr7 != edge.end(); ++itr7) {
                            StripPtR &sp = *itr7;

                            RefsType::reverse_iterator itr7_(itr7);

                            for (; itr7_ != edge.rend(); ++itr7_) {
                                // pt davor
                                StripPtR &pre = *itr7_;

                                if (pre.strip != sp.strip) {
                                    if (pre.strip <= stripInd) {
                                        // der davor wurde bereits verwendet

                                        if (pre.side == SIDE_END) {
                                            sp.ref = pre.desc[0];
                                        } else {
                                            sp.ref = pre.desc[1];
                                        }

                                        break;
                                    }
                                } else {
                                    // gehört dem gleichen strip an
                                    sp.ref = pre.ref;
                                    break;
                                }
                            }


                        }

                        // sonderfall
                        if (edge.size() > 1) {
                            StripPtR &a = edge.front(),
                                &b = *(edge.begin()+1);

                            if (a.ind == b.ind
                                && b.strip == stripInd
                                && pts[a.ind].capt == CAPT_A) {

                                if (b.side == SIDE_START) {
                                    a.ref = b.desc[0];
                                } else {
                                    a.ref = b.desc[1];
                                }

                            }
                        }
                    }

                    // doppelte punkte entfernen

                    for (itr8 = divided.begin(); itr8 != divided.end(); ++itr8) {
#ifdef DEBUG
                        std::cout << "d[" << itr8-divided.begin() << "]" << std::endl;
#endif

                        IdsType &div = *itr8;

                        IdsType _div;

                        int num = div.size();

                        for (int i = 0; i < num; i++) {
                            double ptA[3], ptB[3];
                            pd->GetPoint(div[i], ptA);
                            pd->GetPoint(div[(i+1)%num], ptB);

                            double d = Utilities::GetD3(ptA, ptB);

#ifdef DEBUG
                            std::cout << "edge (" << div[i] << ", "
                                      << div[(i+1)%num] << "), d=" << d
                                      << std::endl;
#endif

                            if (d > 1e-6) {
                                _div.push_back(div[i]);
                            } else {
#ifdef DEBUG
                                std::cout << "rm " << div[i] << std::endl;
#endif
                            }
                        }

                        div = std::move(_div);
                    }

                    // prüfen ob erstellte polygone gültig

                    if (divided[0].size() > 2) {
                        polys.push_back(divided[0]);
                    }

                    if (HasArea(strip) && divided[1].size() > 2) {
                        polys.push_back(divided[1]);
                    }

                    cycle = 0;

                    break;

                } else {
                    polys.push_back(p);

                    cycle++;
                }

            }

            i++;

            if (i == 4) {
                //break;
            }

        }

        // erzeugte polys hinzufügen

        for (itr9 = polys.begin(); itr9 != polys.end(); ++itr9) {
            IdsType &p = *itr9;

            int num = p.size();

            vtkIdList *cell = vtkIdList::New();
            cell->SetNumberOfIds(num);

            for (int i = 0; i < num; i++) {
                cell->SetId(i, p[i]);
            }

            pd->InsertNextCell(VTK_POLYGON, cell);

            origCellIds->InsertNextValue(origCellIds->GetValue(polyInd));

            cell->Delete();
        }

        pd->DeleteCell(polyInd);

    }

    pd->RemoveDeletedCells();
    pd->BuildCells();

}

void vtkPolyDataBooleanFilter::RestoreOrigPoints (vtkPolyData *pd, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "RestoreOrigPoints()" << std::endl;
#endif

    std::vector<StripPtL> ends;

    PolyStripsType::iterator itr;
    StripPtsType::iterator itr2;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
        PStrips &pStrips = itr->second;

        for (itr2 = pStrips.pts.begin(); itr2 != pStrips.pts.end(); ++itr2) {
            StripPt &sp = itr2->second;
            ends.push_back(sp);
        }
    }

    vtkKdTreePointLocator *loc = vtkKdTreePointLocator::New();
    loc->SetDataSet(pd);
    loc->BuildLocator();

    vtkIdList *pts = vtkIdList::New();

    std::vector<StripPtL>::const_iterator itr3;

    for (itr3 = ends.begin(); itr3 != ends.end(); ++itr3) {
        Utilities::FindPoints(loc, itr3->cutPt, pts);
        int numPts = pts->GetNumberOfIds();

        for (int i = 0; i < numPts; i++) {
            pd->GetPoints()->SetPoint(pts->GetId(i), itr3->pt);
        }
    }

    pts->Delete();

    loc->FreeSearchStructure();
    loc->Delete();

}

void vtkPolyDataBooleanFilter::DisjoinPolys (vtkPolyData *pd, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "DisjoinPolys()" << std::endl;
#endif

    pd->BuildLinks();

    std::set<StripPtL> ends;

    PolyStripsType::iterator itr;
    StripPtsType::iterator itr2;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
        PStrips &pStrips = itr->second;

        for (itr2 = pStrips.pts.begin(); itr2 != pStrips.pts.end(); ++itr2) {
            StripPt &sp = itr2->second;
            if (sp.capt == CAPT_A) {
                ends.insert(sp);
            }
        }
    }

    vtkKdTreePointLocator *loc = vtkKdTreePointLocator::New();
    loc->SetDataSet(pd);
    loc->BuildLocator();

    vtkIdList *pts = vtkIdList::New();
    vtkIdList *cells = vtkIdList::New();

    std::set<StripPtL>::const_iterator itr3;

    for (itr3 = ends.begin(); itr3 != ends.end(); ++itr3) {
        Utilities::FindPoints(loc, itr3->pt, pts);
        int numPts = pts->GetNumberOfIds();

        for (int i = 0; i < numPts; i++) {
            pd->GetPointCells(pts->GetId(i), cells);
            int numCells = cells->GetNumberOfIds();

            if (numCells > 1) {
                for (int j = 0; j < numCells; j++) {
                    pd->ReplaceCellPoint(cells->GetId(j), pts->GetId(i), pd->GetPoints()->InsertNextPoint(itr3->pt));
                }
            }
        }
    }

    cells->Delete();
    pts->Delete();

    loc->FreeSearchStructure();
    loc->Delete();

}

void vtkPolyDataBooleanFilter::ResolveOverlaps (vtkPolyData *pd, vtkIntArray *conts, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "ResolveOverlaps()" << std::endl;
#endif

    pd->BuildCells();
    pd->BuildLinks();

    contLines->BuildLinks();

    std::vector<StripPtL2> ends;

    PolyStripsType::iterator itr;
    StripPtsType::iterator itr2;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
        PStrips &pStrips = itr->second;

        for (itr2 = pStrips.pts.begin(); itr2 != pStrips.pts.end(); ++itr2) {
            StripPt &sp = itr2->second;
            if (sp.capt == CAPT_EDGE) {
                ends.push_back(sp);
            }
        }
    }

    vtkIdList *links = vtkIdList::New();

    vtkIdList *cellsA = vtkIdList::New(),
        *cellsB = vtkIdList::New();

    int numCellsA, numCellsB;

    typedef std::map<int, int> CountsType;
    std::map<Pair, CountsType> skipped;

    std::vector<StripPtL2>::const_iterator itr3;

    for (itr3 = ends.begin(); itr3 != ends.end(); ++itr3) {
        contLines->GetPointCells(itr3->ind, links);

        if (links->GetNumberOfIds() == 2
            && conts->GetValue(links->GetId(0)) != conts->GetValue(links->GetId(1))) {

            // kein berührender schnitt

            pd->GetPointCells(itr3->edge[0], cellsA);
            pd->GetPointCells(itr3->edge[1], cellsB);

            numCellsA = cellsA->GetNumberOfIds();
            numCellsB = cellsB->GetNumberOfIds();

            for (int i = 0; i < numCellsA; i++) {
                for (int j = 0; j < numCellsB; j++) {
                    if (cellsA->GetId(i) == cellsB->GetId(j)) {
                        // kante existiert noch

                        vtkIdList *poly = vtkIdList::New();
                        pd->GetCellPoints(cellsA->GetId(i), poly);

                        int numPts = poly->GetNumberOfIds();

                        for (int k = 0; k < numPts; k++) {
                            if (poly->GetId(k) == itr3->edge[1]
                                && poly->GetId((k+1)%numPts) == itr3->edge[0]) {

                                CountsType &c = skipped[Pair(itr3->ind, cellsA->GetId(i))];

                                c[itr3->edge[0]]++;
                                c[itr3->edge[1]]++;
                            }
                        }

                        poly->Delete();

                    }
                }
            }

        }

        links->Reset();
    }

    cellsB->Delete();
    cellsA->Delete();

    links->Delete();

    std::map<Pair, CountsType>::iterator itr4;
    CountsType::iterator itr5;

    for (itr4 = skipped.begin(); itr4 != skipped.end(); ++itr4) {
        const Pair &pair = itr4->first;
        CountsType &c = itr4->second;

#ifdef DEBUG
        std::cout << pair << ": [";
        for (itr5 = c.begin(); itr5 != c.end(); ++itr5) {
            std::cout << "(" << itr5->first << ", " << itr5->second << "), ";
        }
        std::cout << "]" << std::endl;
#endif

        double pt[3];
        contLines->GetPoint(pair.f, pt);

        for (itr5 = c.begin(); itr5 != c.end(); ++itr5) {
            if (itr5->second == 2) {
                int i = pd->GetPoints()->InsertNextPoint(pt);

#ifdef DEBUG
                std::cout << "repl " << itr5->first << " -> " << i << std::endl;
#endif

                pd->ReplaceCellPoint(pair.g, itr5->first, i);
            }
        }
    }
}

void vtkPolyDataBooleanFilter::AddAdjacentPoints (vtkPolyData *pd, vtkIntArray *conts, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "AddAdjacentPoints()" << std::endl;
#endif

    PolyStripsType::iterator itr;
    StripPtsType::iterator itr2;

    std::set<StripPtL2> ends;
    std::set<StripPtL2>::const_iterator itr3;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
        PStrips &pStrips = itr->second;

        for (itr2 = pStrips.pts.begin(); itr2 != pStrips.pts.end(); ++itr2) {
            StripPt &sp = itr2->second;

            if (sp.capt == CAPT_EDGE) {
                ends.insert(sp);
            }
        }
    }

    vtkIntArray *origCellIds = vtkIntArray::SafeDownCast(pd->GetCellData()->GetScalars("OrigCellIds"));

    typedef std::vector<StripPtL3> AType;
    typedef std::map<Pair, AType> BType;

    BType edges;

    for (itr3 = ends.begin(); itr3 != ends.end(); ++itr3) {
        edges[Pair(itr3->edge[0], itr3->edge[1])].push_back(StripPtL3(itr3->pt, itr3->t, itr3->ind));
    }

    AType::iterator itr4;
    BType::iterator itr5;

    pd->BuildLinks();

    vtkKdTreePointLocator *loc = vtkKdTreePointLocator::New();
    loc->SetDataSet(pd);
    loc->BuildLocator();

    typedef std::vector<Pair> CType;

    IdsType::const_iterator itr6;
    CType::const_iterator itr7, itr8;

    vtkIdList *cells = vtkIdList::New();

    for (itr5 = edges.begin(); itr5 != edges.end(); ++itr5) {
        const Pair &pair = itr5->first;
        AType &pts = itr5->second;

        double ptA[3], ptB[3];

        pd->GetPoint(pair.f, ptA);
        pd->GetPoint(pair.g, ptB);

        pts.push_back(StripPtL3(ptA, 0.));
        pts.push_back(StripPtL3(ptB, 1.));

        std::sort(pts.rbegin(), pts.rend());

#ifdef DEBUG
        if (std::find_if(pts.begin(), pts.end(), [](const StripPtL3 &s) {
            return s.ind == 123;
        }) == pts.end()) {
            //continue;
        }

        for (itr4 = pts.begin(); itr4 != pts.end(); ++itr4) {
            std::cout << *itr4 << std::endl;
        }
#endif

        IdsType voids;

        voids.push_back(0);

        for (itr4 = pts.begin()+1; itr4 != pts.end()-1; ++itr4) {
            contLines->GetPointCells(itr4->ind, cells);
            int numCells = cells->GetNumberOfIds();

            std::set<int> seeds;

            for (int i = 0; i < numCells; i++) {
                seeds.insert(cells->GetId(i));
            }

            if (seeds.size() == 1) {
                voids.push_back(itr4-pts.begin());
            }
        }

        voids.push_back(pts.size()-1);

#ifdef DEBUG
        std::cout << "voids=[";
        for (itr6 = voids.begin(); itr6 != voids.end(); ++itr6) {
            std::cout << *itr6 << ", ";
        }
        std::cout << "]" << std::endl;
#endif

        for (itr6 = voids.begin(); itr6 != voids.end()-1; ++itr6) {
            AType pts_(pts.begin()+(*itr6), pts.begin()+(*(itr6+1))+1);

            if (pts_.size() > 2) {

                vtkIdList *ptsA = vtkIdList::New();
                vtkIdList *ptsB = vtkIdList::New();

                Utilities::FindPoints(loc, pts_.front().pt, ptsA);
                Utilities::FindPoints(loc, pts_.back().pt, ptsB);

                int numPtsA = ptsA->GetNumberOfIds(),
                    numPtsB = ptsB->GetNumberOfIds();

                CType polysA, polysB;

                for (int j = 0; j < numPtsA; j++) {
                    pd->GetPointCells(ptsA->GetId(j), cells);
                    int numCells = cells->GetNumberOfIds();
                    for (int k = 0; k < numCells; k++) {
                        polysA.push_back(Pair(cells->GetId(k), ptsA->GetId(j)));
                    }
                }

                for (int j = 0; j < numPtsB; j++) {
                    pd->GetPointCells(ptsB->GetId(j), cells);
                    int numCells = cells->GetNumberOfIds();
                    for (int k = 0; k < numCells; k++) {
                        polysB.push_back(Pair(cells->GetId(k), ptsB->GetId(j)));
                    }
                }

                for (itr7 = polysA.begin(); itr7 != polysA.end(); ++itr7) {
                    for (itr8 = polysB.begin(); itr8 != polysB.end(); ++itr8) {

                        if (itr7->f == itr8->f
                            && pd->GetCellType(itr7->f) != VTK_EMPTY_CELL) {

                            vtkIdList *poly = vtkIdList::New();
                            pd->GetCellPoints(itr7->f, poly);

                            int numPts = poly->GetNumberOfIds();

                            vtkIdList *poly_ = vtkIdList::New();

                            for (int j = 0; j < numPts; j++) {
                                poly_->InsertNextId(poly->GetId(j));

                                if (itr7->g == poly->GetId(j)
                                    && itr8->g == poly->GetId((j+1)%numPts)) {

                                    // ursprüngliche kante

                                    for (itr4 = pts_.begin()+1; itr4 != pts_.end()-1; ++itr4) {
                                        poly_->InsertNextId(pd->InsertNextLinkedPoint(itr4->pt, 1));
                                    }

                                }

                                pd->RemoveReferenceToCell(poly->GetId(j), itr7->f);
                            }

                            pd->DeleteCell(itr7->f);

                            pd->InsertNextLinkedCell(VTK_POLYGON, poly_->GetNumberOfIds(), poly_->GetPointer(0));

                            origCellIds->InsertNextValue(origCellIds->GetValue(itr7->f));

                            poly->Delete();
                            poly_->Delete();

                        }
                    }
                }

                ptsA->Delete();
                ptsB->Delete();

            }

        }

    }

    cells->Delete();

    loc->FreeSearchStructure();
    loc->Delete();

    pd->RemoveDeletedCells();

}

void vtkPolyDataBooleanFilter::MergePoints (vtkPolyData *pd, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "MergePoints()" << std::endl;
#endif

    vtkKdTreePointLocator *loc = vtkKdTreePointLocator::New();
    loc->SetDataSet(pd);
    loc->BuildLocator();

    // essenziell
    pd->BuildLinks();

    PolyStripsType::iterator itr;
    StripsType::iterator itr2;

    typedef std::set<int> AType;
    typedef std::map<int, AType> BType;

    BType inds;

    vtkIdList *pts = vtkIdList::New();

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
        PStrips &pStrips = itr->second;

        for (itr2 = pStrips.strips.begin(); itr2 != pStrips.strips.end(); ++itr2) {
            StripType &strip = *itr2;

            StripPtR &s = strip.front(),
                &e = strip.back();

            Utilities::FindPoints(loc, pStrips.pts[(strip.begin()+1)->ind].pt, pts);
            int numPts = pts->GetNumberOfIds();

            for (int i = 0; i < numPts; i++) {
                inds[s.ind].insert(pts->GetId(i));
            }

            Utilities::FindPoints(loc, pStrips.pts[(strip.end()-2)->ind].pt, pts);
            numPts = pts->GetNumberOfIds();

            for (int i = 0; i < numPts; i++) {
                inds[e.ind].insert(pts->GetId(i));
            }
        }
    }

    vtkIdList *poly = vtkIdList::New(),
        *polys = vtkIdList::New();

    typedef std::vector<MergePt> CType;
    typedef std::vector<Pair> DType;
    typedef std::deque<int> EType;

    BType::iterator itr3;
    CType::iterator itr4, itr5;
    EType::iterator itr6;

#ifdef DEBUG
    AType::const_iterator _itr;
#endif

    for (itr3 = inds.begin(); itr3 != inds.end(); ++itr3) {
        AType &bads = itr3->second;

#ifdef DEBUG
        std::cout << "ind=" << itr3->first << ", bads=[";
        for (_itr = bads.begin(); _itr != bads.end(); ++_itr) {
            std::cout << *_itr << ", ";
        }
        std::cout << "]" << std::endl;
#endif

        CType mergePts;

        double pt[3];
        contLines->GetPoint(itr3->first, pt);

        Utilities::FindPoints(loc, pt, pts);
        int numPts = pts->GetNumberOfIds();

        assert(numPts > 0);

#ifdef DEBUG
        std::cout << "pts=[";
        for (int i = 0; i < numPts; i++) {
            int ind = pts->GetId(i);
            pd->GetPointCells(ind, polys);
            if (polys->GetNumberOfIds() > 0) {
                std::cout << ind << " -> " << polys->GetId(0) << ", ";
            }
        }
        std::cout << "]" << std::endl;
#endif

        for (int i = 0; i < numPts; i++) {
            int ind = pts->GetId(i);

            pd->GetPointCells(ind, polys);

            if (polys->GetNumberOfIds() > 0) {
                // sollte nur von einem verwendet werden

                pd->GetCellPoints(polys->GetId(0), poly);
                int numPts_ = poly->GetNumberOfIds();

                int j = 0;
                for (; j < numPts_; j++) {
                    if (poly->GetId(j) == ind) {
                        break;
                    }
                }

                // davor und danach
                int before = poly->GetId((j+numPts_-1)%numPts_),
                    after = poly->GetId((j+1)%numPts_);

                double pt_[3];

                if (std::count(bads.begin(), bads.end(), before) == 0) {
                    pd->GetPoint(before, pt_);
                    mergePts.push_back(MergePt(polys->GetId(0), ind, pt_));
                }

                if (std::count(bads.begin(), bads.end(), after) == 0) {
                    pd->GetPoint(after, pt_);
                    mergePts.push_back(MergePt(polys->GetId(0), ind, pt_));
                }

            }
        }

#ifdef DEBUG
        for (itr4 = mergePts.begin(); itr4 != mergePts.end(); ++itr4) {
            std::cout << (itr4-mergePts.begin())
                      << ": (" << itr4->polyInd << ", " << itr4->ind
                      << ", [" << itr4->pt[0] << ", " << itr4->pt[1] << ", " << itr4->pt[2]
                      << "])" << std::endl;
        }
#endif

        // benachbarte polygone finden

        DType pairs;

        for (itr4 = mergePts.begin(); itr4 != mergePts.end(); ++itr4) {
            for (itr5 = itr4+1; itr5 != mergePts.end(); ++itr5) {
                if (Utilities::GetD3(itr4->pt, itr5->pt) < 1e-6) {
                    // paar gefunden
                    pairs.push_back(Pair(itr4-mergePts.begin(), itr5-mergePts.begin()));
                }
            }
        }

#ifdef DEBUG
        DType::const_iterator _itr2;
        for (_itr2 = pairs.begin(); _itr2 != pairs.end(); ++_itr2) {
            std::cout << *_itr2 << std::endl;
        }
#endif

        // gruppiert anhand von aneinandergrenzung

        EType group;

        int i = 0;

        while (pairs.size() > 0) {
            Pair &pair = pairs[i];

            if (group.empty()) {
                group.push_back(pair.f);
                group.push_back(pair.g);

                pairs.erase(pairs.begin());
            } else {
                MergePt &a = mergePts[pair.f],
                    &b = mergePts[pair.g];

                if (mergePts[group.front()].ind == a.ind) {
                    group.push_front(pair.g);
                    pairs.erase(pairs.begin()+i);
                    i = 0;
                } else if (mergePts[group.front()].ind == b.ind) {
                    group.push_front(pair.f);
                    pairs.erase(pairs.begin()+i);
                    i = 0;
                } else if (mergePts[group.back()].ind == a.ind) {
                    group.push_back(pair.g);
                    pairs.erase(pairs.begin()+i);
                    i = 0;
                } else if (mergePts[group.back()].ind == b.ind) {
                    group.push_back(pair.f);
                    pairs.erase(pairs.begin()+i);
                    i = 0;
                } else {
                    i++;
                }

            }

            if (i == pairs.size()) {
                for (itr6 = group.begin()+1; itr6 != group.end(); ++itr6) {
#ifdef DEBUG
                    std::cout << "repl " <<  mergePts[*itr6].ind << " -> " << mergePts[group.front()].ind << std::endl;
#endif

                    pd->ReplaceCellPoint(mergePts[*itr6].polyInd, mergePts[*itr6].ind, mergePts[group.front()].ind);
                }

                group.clear();
                i = 0;
            }

        }

    }

    pts->Delete();
    poly->Delete();
    polys->Delete();

    loc->FreeSearchStructure();
    loc->Delete();

}

void vtkPolyDataBooleanFilter::DecomposePolys () {

#ifdef DEBUG
    std::cout << "DecomposePolys()" << std::endl;
#endif

    vtkPolyData *pd = vtkPolyData::New();
    pd->DeepCopy(resultA);

    vtkCellData *cellData = pd->GetCellData();

    vtkIntArray *origCellIdsA = vtkIntArray::SafeDownCast(cellData->GetScalars("OrigCellIdsA")),
        *origCellIdsB = vtkIntArray::SafeDownCast(cellData->GetScalars("OrigCellIdsB"));

    int numCells = pd->GetNumberOfCells();

    vtkIdList *cells = vtkIdList::New();
    for (int i = 0; i < numCells; i++) {
        if (involvedA.count(origCellIdsA->GetValue(i)) > 0
            || involvedB.count(origCellIdsB->GetValue(i)) > 0) {
            cells->InsertNextId(i);
        }
    }

#ifdef _DEBUG

    vtkExtractCells *extract = vtkExtractCells::New();
    extract->SetInputData(pd);
    extract->SetCellList(cells);

    vtkGeometryFilter *gf = vtkGeometryFilter::New();
    gf->SetInputConnection(extract->GetOutputPort());

    vtkCleanPolyData *clean = vtkCleanPolyData::New();
    clean->PointMergingOff();
    clean->SetInputConnection(gf->GetOutputPort());
    clean->Update();

    Utilities::WriteVTK("extracted.vtk", clean->GetOutput());

    clean->Delete();
    gf->Delete();
    extract->Delete();
    cells->Delete();

#else

    for (int i = 0; i < cells->GetNumberOfIds(); i++) {
        vtkIdList *poly = vtkIdList::New();
        pd->GetCellPoints(cells->GetId(i), poly);

        int numPts = poly->GetNumberOfIds();

        /*double pts[numPts][3],
            pts2[numPts][2];*/

		std::vector<std::array<double, 3> > pts(numPts);
		std::vector<std::array<double, 2> > pts2(numPts);

        for (int j = 0; j < numPts; j++) {
            pd->GetPoint(poly->GetId(j), pts[j].data());
        }

        // transformiert in 2d
        Base base(pts, numPts);
        Transform(pts, pts2, numPts, base);

        // zerlegt
        DecType dec;
        Decompose(pts2, numPts).collect(dec);

        if (dec.size() > 1) {
            vtkIdList *poly2 = vtkIdList::New();

            DecType::const_iterator itr;
            IdsType::const_iterator itr2;
            for (itr = dec.begin(); itr != dec.end(); ++itr) {
                poly2->Reset();
                for (itr2 = itr->begin(); itr2 != itr->end(); ++itr2) {
                    poly2->InsertNextId(poly->GetId(*itr2));
                }

                pd->InsertNextCell(VTK_POLYGON, poly2);

                origCellIdsA->InsertNextValue(origCellIdsA->GetValue(cells->GetId(i)));
                origCellIdsB->InsertNextValue(origCellIdsB->GetValue(cells->GetId(i)));
            }

            poly2->Delete();

            // löscht das nicht-konvexe polygon
            pd->DeleteCell(cells->GetId(i));
        }

        poly->Delete();

    }

    pd->RemoveDeletedCells();

#endif

    resultA->ShallowCopy(pd);

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

        int numPts = poly->GetNumberOfIds();

        for (int i = 0; i < numPts; i++) {
            if (poly->GetId(i) == ptIds[0]) {
                if (poly->GetId((i+1)%numPts) != ptIds[1]) {
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

        Utilities::ComputeNormal(pd->GetPoints(), poly, n);

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
        alpha = Utilities::GetAngle(pA->r, pB->r, pA->e);
    }

    ~PolyPair () {
        delete pA;
        delete pB;
    }

    PolyAtEdge *pA, *pB;

    void GetLoc (PolyAtEdge *pT, int mode) {
        double beta = Utilities::GetAngle(pA->r, pT->r, pA->e);

#ifdef DEBUG
        std::cout << "GetLoc() -> polyId "
                  << pT->polyId << ", beta " << (beta*180/pi) << std::endl;
#endif

        if (beta < 1e-7 || beta > 2*pi-1e-7) {
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

        } else if (std::abs(beta-alpha) < 1e-7) {

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

    int numPtsA = ptsA->GetNumberOfIds(),
        numPtsB = ptsB->GetNumberOfIds();

    vtkIdList *polys = vtkIdList::New();

    for (int i = 0; i < numPtsA; i++) {
        polys->Reset();

        pd->GetPointCells(ptsA->GetId(i), polys);
        int numCells = polys->GetNumberOfIds();

        for (int j = 0; j < numCells; j++) {
            polyPts[polys->GetId(j)] = ptsA->GetId(i);
        }
    }

    PolyAtEdge *first = NULL;

    for (int i = 0; i < numPtsB; i++) {
        polys->Reset();

        pd->GetPointCells(ptsB->GetId(i), polys);
        int numCells = polys->GetNumberOfIds();

        for (int j = 0; j < numCells; j++) {

            if (polyPts.count(polys->GetId(j)) > 0) {

                // die kante muss auch existieren

                bool adj = false;

                vtkIdList *poly = vtkIdList::New();
                pd->GetCellPoints(polys->GetId(j), poly);

                int numPts = poly->GetNumberOfIds();

                for (int k = 0; k < numPts; k++) {
                    if (poly->GetId(k) == ptsB->GetId(i)
                        && (polyPts[polys->GetId(j)] == poly->GetId((k+1)%numPts)
                            || polyPts[polys->GetId(j)] == poly->GetId((k+numPts-1)%numPts))) {
                        adj = true;

                        break;
                    }
                }

                poly->Delete();

                if (adj) {
                    if (first != NULL) {
                        PolyAtEdge *second = new PolyAtEdge(pd, polys->GetId(j), polyPts[polys->GetId(j)], ptsB->GetId(i));

                        pp = new PolyPair(first, second);

                        break;
                    } else {
                        first = new PolyAtEdge(pd, polys->GetId(j), polyPts[polys->GetId(j)], ptsB->GetId(i));
                    }
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

        std::cout << "alpha " << (pp->alpha*180/pi) << std::endl;

    }
#endif

    return pp;

}


void vtkPolyDataBooleanFilter::CombineRegions () {

#ifdef DEBUG
    std::cout << "CombineRegions()" << std::endl;
#endif

    // double-Koordinaten sichern
    modPdA->GetPointData()->AddArray(modPdA->GetPoints()->GetData());
    modPdB->GetPointData()->AddArray(modPdB->GetPoints()->GetData());

    // ungenutzte punkte löschen
    vtkCleanPolyData *cleanA = vtkCleanPolyData::New();
    cleanA->PointMergingOff();
    cleanA->SetInputData(modPdA);

    vtkCleanPolyData *cleanB = vtkCleanPolyData::New();
    cleanB->PointMergingOff();
    cleanB->SetInputData(modPdB);

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

    vtkPolyData *pdA = cfA->GetOutput();
    vtkPolyData *pdB = cfB->GetOutput();

    // double-Pts wiederherstellen
    pdA->GetPoints()->SetData(pdA->GetPointData()->GetArray(0));
    pdB->GetPoints()->SetData(pdB->GetPointData()->GetArray(0));

#ifdef DEBUG
    std::cout << "Exporting modPdA_8.vtk" << std::endl;
    Utilities::WriteVTK("modPdA_8.vtk", cfA->GetOutput());

    std::cout << "Exporting modPdB_8.vtk" << std::endl;
    Utilities::WriteVTK("modPdB_8.vtk", cfB->GetOutput());
#endif

    // locators erstellen
    vtkKdTreePointLocator *plA = vtkKdTreePointLocator::New();
    plA->SetDataSet(pdA);
    plA->BuildLocator();

    vtkKdTreePointLocator *plB = vtkKdTreePointLocator::New();
    plB->SetDataSet(pdB);
    plB->BuildLocator();

    pdA->BuildLinks();
    pdB->BuildLinks();

    vtkDataArray *scalarsA = pdA->GetPointData()->GetScalars();
    vtkDataArray *scalarsB = pdB->GetPointData()->GetScalars();

    vtkIdList *line = vtkIdList::New();

    double ptA[3], ptB[3];

    vtkIdList *fptsA = vtkIdList::New();
    vtkIdList *lptsA = vtkIdList::New();

    vtkIdList *fptsB = vtkIdList::New();
    vtkIdList *lptsB = vtkIdList::New();

    std::map<int, int> locsA, locsB;

    for (int i = 0; i < contLines->GetNumberOfCells(); i++) {

        contLines->GetCellPoints(i, line);

        contLines->GetPoint(line->GetId(0), ptA);
        contLines->GetPoint(line->GetId(1), ptB);

        Utilities::FindPoints(plA, ptA, fptsA);
        Utilities::FindPoints(plB, ptA, fptsB);

#ifdef DEBUG
        std::cout << "line " << i << std::endl;
#else

        // bereits behandelte regionen werden nicht noch einmal untersucht

        int notLocated = 0;

        for (int j = 0; j < fptsA->GetNumberOfIds(); j++) {
            if (locsA.count(scalarsA->GetTuple1(fptsA->GetId(j))) == 0) {
                notLocated++;
            }
        }

        for (int j = 0; j < fptsB->GetNumberOfIds(); j++) {
            if (locsB.count(scalarsB->GetTuple1(fptsB->GetId(j))) == 0) {
                notLocated++;
            }
        }

        if (notLocated == 0) {
            continue;
        }

#endif

        Utilities::FindPoints(plA, ptB, lptsA);
        Utilities::FindPoints(plB, ptB, lptsB);

        PolyPair *ppA = GetEdgePolys(pdA, fptsA, lptsA);
        PolyPair *ppB = GetEdgePolys(pdB, fptsB, lptsB);

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

            if (locsA.count(fsA) > 0 && locsA[fsA] != ppA->pA->loc) {
                std::cout << "sA " << fsA << ": " << locsA[fsA] << " -> " << ppA->pA->loc << std::endl;
            }

            if (locsA.count(lsA) > 0 && locsA[lsA] != ppA->pB->loc) {
                std::cout << "sA " << lsA << ": " << locsA[lsA] << " -> " << ppA->pB->loc << std::endl;
            }

            if (locsB.count(fsB) > 0 && locsB[fsB] != ppB->pA->loc) {
                std::cout << "sB " << fsB << ": " << locsB[fsB] << " -> " << ppB->pA->loc << std::endl;
            }

            if (locsB.count(lsB) > 0 && locsB[lsB] != ppB->pB->loc) {
                std::cout << "sB " << lsB << ": " << locsB[lsB] << " -> " << ppB->pB->loc << std::endl;
            }

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

    }

    lptsB->Delete();
    fptsB->Delete();

    lptsA->Delete();
    fptsA->Delete();

    line->Delete();

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
        for (int i = 0; i < regsA->GetNumberOfCells(); i++) {
            regsA->ReverseCell(i);
        }
    }

    cfB->Update();

    vtkPolyData *regsB = cfB->GetOutput();

    if (comb[1] == LOC_INSIDE && OperMode != OPER_INTERSECTION) {
        for (int i = 0; i < regsB->GetNumberOfCells(); i++) {
            regsB->ReverseCell(i);
        }
    }

    // OrigCellIds und CellData

    vtkIntArray *origCellIdsA = vtkIntArray::SafeDownCast(regsA->GetCellData()->GetScalars("OrigCellIds"));
    vtkIntArray *origCellIdsB = vtkIntArray::SafeDownCast(regsB->GetCellData()->GetScalars("OrigCellIds"));

    vtkIntArray *newOrigCellIdsA = vtkIntArray::New();
    newOrigCellIdsA->SetName("OrigCellIdsA");

    vtkIntArray *newOrigCellIdsB = vtkIntArray::New();
    newOrigCellIdsB->SetName("OrigCellIdsB");

    vtkCellData *newCellDataA = vtkCellData::New();
    newCellDataA->CopyAllocate(cellDataA);

    vtkCellData *newCellDataB = vtkCellData::New();
    newCellDataB->CopyAllocate(cellDataB);

    for (int i = 0; i < regsA->GetNumberOfCells(); i++) {
        newOrigCellIdsA->InsertNextValue(origCellIdsA->GetValue(i));
        newOrigCellIdsB->InsertNextValue(-1);

        newCellDataA->CopyData(cellDataA, origCellIdsA->GetValue(i), i);
    }

    for (int i = 0; i < regsB->GetNumberOfCells(); i++) {
        newOrigCellIdsB->InsertNextValue(origCellIdsB->GetValue(i));
        newOrigCellIdsA->InsertNextValue(-1);

        newCellDataB->CopyData(cellDataB, origCellIdsB->GetValue(i), i);
    }

    regsA->GetCellData()->Initialize();
    regsB->GetCellData()->Initialize();

    regsA->GetCellData()->ShallowCopy(newCellDataA);
    regsB->GetCellData()->ShallowCopy(newCellDataB);

    newCellDataA->Delete();
    newCellDataB->Delete();

    // zusammenführung

    vtkAppendPolyData *app = vtkAppendPolyData::New();
    app->AddInputData(regsA);
    app->AddInputData(regsB);

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

    vtkPolyData *cfPd = cfApp->GetOutput();

    // double-Koordinaten wiederherstellen
    if (cfPd->GetPoints() != NULL) {
        cfPd->GetPoints()->SetData(cfPd->GetPointData()->GetArray(0));
        cfPd->GetPointData()->RemoveArray(0);
    }

    // resultA ist erster output des filters
    resultA->ShallowCopy(cfPd);

    resultA->GetCellData()->AddArray(newOrigCellIdsA);
    resultA->GetCellData()->AddArray(newOrigCellIdsB);

    resultB->ShallowCopy(contLines);

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


void vtkPolyDataBooleanFilter::MergeRegions () {

#ifdef DEBUG
    std::cout << "MergeRegions()" << std::endl;
#endif

    vtkIntArray *origCellIdsA = vtkIntArray::New();
    origCellIdsA->DeepCopy(modPdA->GetCellData()->GetScalars("OrigCellIds"));
    origCellIdsA->SetName("OrigCellIdsA");

    vtkIntArray *origCellIdsB = vtkIntArray::New();
    origCellIdsB->DeepCopy(modPdB->GetCellData()->GetScalars("OrigCellIds"));
    origCellIdsB->SetName("OrigCellIdsB");

    vtkPolyData *pdA = vtkPolyData::New();
    pdA->CopyStructure(modPdA);

    vtkPolyData *pdB = vtkPolyData::New();
    pdB->CopyStructure(modPdB);

    pdA->GetCellData()->AddArray(origCellIdsA);
    pdB->GetCellData()->AddArray(origCellIdsB);

    origCellIdsA->Delete();
    origCellIdsB->Delete();

    vtkIntArray *padIdsA = vtkIntArray::New();
    vtkIntArray *padIdsB = vtkIntArray::New();

    for (int i = 0; i < pdA->GetNumberOfCells(); i++) {
        padIdsA->InsertNextValue(-1);
    }

    for (int i = 0; i < pdB->GetNumberOfCells(); i++) {
        padIdsB->InsertNextValue(-1);
    }

    padIdsA->SetName("OrigCellIdsB");
    padIdsB->SetName("OrigCellIdsA");

    pdA->GetCellData()->AddArray(padIdsA);
    pdB->GetCellData()->AddArray(padIdsB);

    padIdsA->Delete();
    padIdsB->Delete();

    vtkAppendPolyData *app = vtkAppendPolyData::New();
    app->AddInputData(pdA);
    app->AddInputData(pdB);

    app->Update();

    vtkPolyData *appPd = app->GetOutput();

    appPd->GetPointData()->AddArray(appPd->GetPoints()->GetData());

    vtkCleanPolyData *clean = vtkCleanPolyData::New();
    clean->PointMergingOff();
    clean->SetInputData(appPd);

    clean->Update();

    vtkPolyData *cleanPd = clean->GetOutput();

    cleanPd->GetPoints()->SetData(cleanPd->GetPointData()->GetArray(0));
    cleanPd->GetPointData()->RemoveArray(0);

    resultA->ShallowCopy(cleanPd);
    resultB->ShallowCopy(contLines);

    clean->Delete();
    app->Delete();

    pdB->Delete();
    pdA->Delete();

}
