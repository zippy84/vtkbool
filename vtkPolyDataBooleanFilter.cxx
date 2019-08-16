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

#include <map>
#include <deque>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <functional>
#include <queue>
#include <sstream>

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
#include <vtkSmartPointer.h>

#include "vtkPolyDataBooleanFilter.h"
#include "vtkPolyDataContactFilter.h"

#include "Utilities.h"

#include "Merger.h"
#include "Decomposer.h"
#include "AABB.h"

#ifdef DEBUG
#include <chrono>
#include <numeric>
#include <iterator>
#endif

// #include <csignal>

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

    cellIdsA = vtkIntArray::New();
    cellIdsB = vtkIntArray::New();

    OperMode = OPER_UNION;

    MergeRegs = false;
    DecPolys = true;

}

vtkPolyDataBooleanFilter::~vtkPolyDataBooleanFilter () {

    cellIdsA->Delete();
    cellIdsB->Delete();

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
        using clock = std::chrono::steady_clock;
        std::vector<clock::duration> times;
        clock::time_point start;
#endif

        if (pdA->GetMTime() > timePdA || pdB->GetMTime() > timePdB) {

            // eventuell vorhandene regionen vereinen

            vtkSmartPointer<vtkCleanPolyData> cleanA = vtkSmartPointer<vtkCleanPolyData>::New();
            cleanA->SetOutputPointsPrecision(DOUBLE_PRECISION);
            cleanA->SetTolerance(1e-6);
            cleanA->SetInputData(pdA);
            cleanA->Update();

            vtkSmartPointer<vtkCleanPolyData> cleanB = vtkSmartPointer<vtkCleanPolyData>::New();
            cleanB->SetOutputPointsPrecision(DOUBLE_PRECISION);
            cleanB->SetTolerance(1e-6);
            cleanB->SetInputData(pdB);
            cleanB->Update();

#ifdef DEBUG
            std::cout << "Exporting modPdA.vtk" << std::endl;
            WriteVTK("modPdA.vtk", cleanA->GetOutput());

            std::cout << "Exporting modPdB.vtk" << std::endl;
            WriteVTK("modPdB.vtk", cleanB->GetOutput());
#endif

            // CellData sichern

            cellDataA->DeepCopy(cleanA->GetOutput()->GetCellData());
            cellDataB->DeepCopy(cleanB->GetOutput()->GetCellData());

            // ermittelt kontaktstellen

#ifdef DEBUG

            start = clock::now();
#endif

            vtkSmartPointer<vtkPolyDataContactFilter> cl = vtkSmartPointer<vtkPolyDataContactFilter>::New();
            cl->SetInputConnection(0, cleanA->GetOutputPort());
            cl->SetInputConnection(1, cleanB->GetOutputPort());
            cl->Update();

#ifdef DEBUG
            times.push_back(clock::now()-start);
#endif

            contLines->DeepCopy(cl->GetOutput());

#ifdef DEBUG
            std::cout << "Exporting contLines.vtk" << std::endl;
            WriteVTK("contLines.vtk", contLines);

            std::cout << "Exporting modPdA_1.vtk" << std::endl;
            WriteVTK("modPdA_1.vtk", cl->GetOutput(1));

            std::cout << "Exporting modPdB_1.vtk" << std::endl;
            WriteVTK("modPdB_1.vtk", cl->GetOutput(2));
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

            vtkIntArray *sourcesA = vtkIntArray::SafeDownCast(contLines->GetCellData()->GetScalars("sourcesA"));
            vtkIntArray *sourcesB = vtkIntArray::SafeDownCast(contLines->GetCellData()->GetScalars("sourcesB"));

            int i, numPts = contLines->GetNumberOfPoints();

            vtkIdList *cells = vtkIdList::New();

            for (i = 0; i < numPts; i++) {
                contLines->GetPointCells(i, cells);

                if (cells->GetNumberOfIds() == 1) {
                    break;
                }
            }

            cells->Delete();

            if (i < numPts) {
                vtkErrorMacro("Contact ends suddenly at point " << i << ".");

                return 1;
            }

            // sichert die OrigCellIds

            vtkIntArray *origCellIdsA = vtkIntArray::SafeDownCast(modPdA->GetCellData()->GetScalars("OrigCellIds"));
            vtkIntArray *origCellIdsB = vtkIntArray::SafeDownCast(modPdB->GetCellData()->GetScalars("OrigCellIds"));

            cellIdsA->DeepCopy(origCellIdsA);
            cellIdsB->DeepCopy(origCellIdsB);

            for (int i = 0; i < modPdA->GetNumberOfCells(); i++) {
                origCellIdsA->SetValue(i, i);
            }

            for (int i = 0; i < modPdB->GetNumberOfCells(); i++) {
                origCellIdsB->SetValue(i, i);
            }


#ifdef DEBUG
            start = clock::now();
#endif

            if (GetPolyStrips(modPdA, contsA, sourcesA, polyStripsA) ||
                GetPolyStrips(modPdB, contsB, sourcesB, polyStripsB)) {

                vtkErrorMacro("Strips are invalid.");

                return 1;

            }

#ifdef DEBUG
            times.push_back(clock::now()-start);
#endif

            // löst ein sehr spezielles problem

#ifdef DEBUG
            start = clock::now();
#endif

            CollapseCaptPoints(modPdA, polyStripsA);
            CollapseCaptPoints(modPdB, polyStripsB);

#ifdef DEBUG
            times.push_back(clock::now()-start);
#endif

            // trennt die polygone an den linien

#ifdef DEBUG
            start = clock::now();
#endif

            CutCells(modPdA, polyStripsA);
            CutCells(modPdB, polyStripsB);

#ifdef DEBUG
            times.push_back(clock::now()-start);
#endif

#ifdef DEBUG
            std::cout << "Exporting modPdA_2.vtk" << std::endl;
            WriteVTK("modPdA_2.vtk", modPdA);

            std::cout << "Exporting modPdB_2.vtk" << std::endl;
            WriteVTK("modPdB_2.vtk", modPdB);
#endif

#ifdef DEBUG
            start = clock::now();
#endif

            RestoreOrigPoints(modPdA, polyStripsA);
            RestoreOrigPoints(modPdB, polyStripsB);

#ifdef DEBUG
            times.push_back(clock::now()-start);
#endif

#ifdef DEBUG
            std::cout << "Exporting modPdA_3.vtk" << std::endl;
            WriteVTK("modPdA_3.vtk", modPdA);

            std::cout << "Exporting modPdB_3.vtk" << std::endl;
            WriteVTK("modPdB_3.vtk", modPdB);
#endif

#ifdef DEBUG
            start = clock::now();
#endif

            ResolveOverlaps(modPdA, contsA, polyStripsA);
            ResolveOverlaps(modPdB, contsB, polyStripsB);

#ifdef DEBUG
            times.push_back(clock::now()-start);
#endif

#ifdef DEBUG
            std::cout << "Exporting modPdA_4.vtk" << std::endl;
            WriteVTK("modPdA_4.vtk", modPdA);

            std::cout << "Exporting modPdB_4.vtk" << std::endl;
            WriteVTK("modPdB_4.vtk", modPdB);
#endif

#ifdef DEBUG
            start = clock::now();
#endif

            AddAdjacentPoints(modPdA, contsA, polyStripsA);
            AddAdjacentPoints(modPdB, contsB, polyStripsB);

#ifdef DEBUG
            times.push_back(clock::now()-start);
#endif

#ifdef DEBUG
            std::cout << "Exporting modPdA_5.vtk" << std::endl;
            WriteVTK("modPdA_5.vtk", modPdA);

            std::cout << "Exporting modPdB_5.vtk" << std::endl;
            WriteVTK("modPdB_5.vtk", modPdB);
#endif

#ifdef DEBUG
            start = clock::now();
#endif

            DisjoinPolys(modPdA, polyStripsA);
            DisjoinPolys(modPdB, polyStripsB);

#ifdef DEBUG
            times.push_back(clock::now()-start);
#endif

#ifdef DEBUG
            std::cout << "Exporting modPdA_6.vtk" << std::endl;
            WriteVTK("modPdA_6.vtk", modPdA);

            std::cout << "Exporting modPdB_6.vtk" << std::endl;
            WriteVTK("modPdB_6.vtk", modPdB);
#endif

#ifdef DEBUG
            start = clock::now();
#endif

            MergePoints(modPdA, polyStripsA);
            MergePoints(modPdB, polyStripsB);

#ifdef DEBUG
            times.push_back(clock::now()-start);
#endif

#ifdef DEBUG
            std::cout << "Exporting modPdA_7.vtk" << std::endl;
            WriteVTK("modPdA_7.vtk", modPdA);

            std::cout << "Exporting modPdB_7.vtk" << std::endl;
            WriteVTK("modPdB_7.vtk", modPdB);
#endif

            involvedA.clear();
            involvedB.clear();

            int numLines = contLines->GetNumberOfCells();

            for (int i = 0; i < numLines; i++) {
                involvedA.insert(contsA->GetValue(i));
                involvedB.insert(contsB->GetValue(i));
            }

            relsA.clear();
            relsB.clear();

            // aufräumen

            /*cl->Delete();
            cleanB->Delete();
            cleanA->Delete();*/

            timePdA = pdA->GetMTime();
            timePdB = pdB->GetMTime();

        }

#ifdef DEBUG
        times.push_back(clock::now()-start);
#endif

        DecPolys_(modPdA, involvedA, relsA);
        DecPolys_(modPdB, involvedB, relsB);

#ifdef DEBUG
        times.push_back(clock::now()-start);
#endif

#ifdef DEBUG
        std::cout << "Exporting modPdA_8.vtk" << std::endl;
        WriteVTK("modPdA_8.vtk", modPdA);

        std::cout << "Exporting modPdB_8.vtk" << std::endl;
        WriteVTK("modPdB_8.vtk", modPdB);
#endif

#ifdef DEBUG
        start = clock::now();
#endif

        if (MergeRegs) {
            MergeRegions();
        } else {
            CombineRegions();
        }

#ifdef DEBUG
        times.push_back(clock::now()-start);
#endif


#ifdef DEBUG
        double sum = std::chrono::duration_cast<std::chrono::duration<double>>(std::accumulate(times.begin(), times.end(), clock::duration())).count();

        std::vector<clock::duration>::const_iterator itr;
        for (itr = times.begin(); itr != times.end(); itr++) {
            double time = std::chrono::duration_cast<std::chrono::duration<double>>(*itr).count();

            std::cout << "Time " << (itr-times.begin())
                << ": " << time << "s (" << (time/sum) << "%)"
                << std::endl;
        }
#endif

    }

    return 1;

}

void vtkPolyDataBooleanFilter::GetStripPoints (vtkPolyData *pd, vtkIntArray *sources, PStrips &pStrips, IdsType &lines) {

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

                Cpy(pts[realInd].pt, pt, 3);

                int src = sources->GetComponent(*itr, i);

                pts[realInd].src = src;

                double lastD = DBL_MAX;

                // jetzt muss man die kanten durchlaufen
                for (int j = 0; j < numPts; j++) {

                    if (src > -1 && !(j == src || (j+1)%numPts == src)) {
                        continue;
                    }

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
                            Cpy(pts[realInd].captPt, a, 3);
                            pts[realInd].capt = CAPT_A;

                            break;

                        } else if (vtkMath::Norm(sB) < tol) {
                            Cpy(pts[realInd].captPt, b, 3);
                            pts[realInd].capt = CAPT_B;

                            break;

                        } else {
                            // u ist nicht normiert
                            vtkMath::MultiplyScalar(u, t);

                            double x[3];
                            vtkMath::Add(a, u, x);

                            // projektion
                            Cpy(pts[realInd].captPt, x, 3);

                            pts[realInd].capt = CAPT_EDGE;

                            lastD = d;
                        }

                    }

                }

                if (src != NO_USE && pts[realInd].edge[0] == NO_USE) {
                    pts[realInd].catched = false;
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

            Cpy(sp.cutPt, sp.captPt, 3);
        } else {

            Cpy(sp.cutPt, sp.pt, 3);
        }

        sp.history.push_back({sp.edge[0], sp.edge[1]});

    }

#ifdef DEBUG
    std::cout << "pts: " << std::endl;
    for (itr2 = pts.begin(); itr2 != pts.end(); ++itr2) {
        std::cout << itr2->first << ": " << itr2->second << std::endl;
    }
#endif

}

bool vtkPolyDataBooleanFilter::GetPolyStrips (vtkPolyData *pd, vtkIntArray *conts, vtkIntArray *sources, PolyStripsType &polyStrips) {
#ifdef DEBUG
    std::cout << "GetPolyStrips()" << std::endl;
#endif

    polyStrips.clear();

    std::map<int, IdsType> polyLines;

    for (int i = 0; i < conts->GetNumberOfTuples(); i++) {
        int poly = conts->GetValue(i);
        polyLines[poly].push_back(i);
    }

    std::vector<std::reference_wrapper<StripPt>> notCatched;

    std::map<int, IdsType>::iterator itr;

    for (itr = polyLines.begin(); itr != polyLines.end(); ++itr) {
        IdsType &lines = itr->second;
        RemoveDuplicates(lines);

        PStrips &pStrips = polyStrips[itr->first];

        vtkIdList *polyPts = vtkIdList::New();
        pd->GetCellPoints(itr->first, polyPts);

        int numPts = polyPts->GetNumberOfIds();

        for (int i = 0; i < numPts; i++) {
            pStrips.poly.push_back(polyPts->GetId(i));
        }

        ComputeNormal(pd->GetPoints(), pStrips.n, polyPts);

        GetStripPoints(pd, sources, pStrips, lines);

        for (auto &sp : pStrips.pts) {
            sp.second.polyId = itr->first;

            if (!sp.second.catched) {
                notCatched.push_back(sp.second);
            }
        }

        polyPts->Delete();
    }

    auto Next = [](const IdsType &ids, int id) -> int {
        IdsType::const_iterator itr;

        itr = std::find(ids.begin(), ids.end(), id);

        if (++itr == ids.end()) {
            itr = ids.begin();
        }

        return *itr;
    };

    for (StripPt &sp : notCatched) {
        for (itr = polyLines.begin(); itr != polyLines.end(); ++itr) {
            const PStrips &pStrips = polyStrips[itr->first];

            try {
                const StripPt &corr = pStrips.pts.at(sp.ind);

                if (&corr != &sp) {
                    if (corr.capt == CAPT_A) {
                        sp.capt = CAPT_A;
                        sp.edge[0] = corr.edge[0];
                        sp.edge[1] = Next(polyStrips[sp.polyId].poly, sp.edge[0]);

                        sp.t = 0;

                        Cpy(sp.captPt, corr.captPt, 3);
                        Cpy(sp.cutPt, sp.captPt, 3);

                        sp.history.push_back({sp.edge[0], sp.edge[1]});

                        sp.catched = true;

                    }
                }
            } catch (...) {}

        }

        assert(sp.catched);

    }

    StripPtsType::const_iterator itr5;

    for (itr = polyLines.begin(); itr != polyLines.end(); ++itr) {
        PStrips &pStrips = polyStrips[itr->first];

        const IdsType &lines = itr->second;
        const StripPtsType &pts = pStrips.pts;

        for (itr5 = pts.begin(); itr5 != pts.end(); ++itr5) {
            // ist der punkt auf einer kante, dürfen von ihm mehr als 2 linien ausgehen

            const StripPt &pt = itr5->second;

            if (pt.capt == CAPT_NOT) {

                vtkIdList *cells = vtkIdList::New(),
                    *line = vtkIdList::New();

                contLines->GetPointCells(pt.ind, cells);

                int numCells = cells->GetNumberOfIds();

                std::set<int> ends;

                for (int i = 0; i < numCells; i++) {
                    contLines->GetCellPoints(cells->GetId(i), line);

                    ends.insert(pt.ind == line->GetId(0) ? line->GetId(1) : line->GetId(0));
                }

                line->Delete();
                cells->Delete();

                if (ends.size() > 2) {
                    return true;
                }
            }
        }

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
                const StripPt &start = pts.at(strip.front().ind),
                    &end = pts.at(strip.back().ind);

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

            const StripPt &_start = pts.at(strip.front().ind),
                &_end = pts.at(strip.back().ind);

            if ((_start.capt != CAPT_NOT && _end.capt != CAPT_NOT)
                || (_start.ind == _end.ind)
                || (i == _lines.size())) {

                // einen neuen strip anlegen
                pStrips.strips.push_back(strip);
                strip.clear();

                i = 0;
            }

            linePts->Delete();

        }

        CompleteStrips(pStrips);

    }

    // validierung
    // sucht nach schnitten zw. den strips

    /*
    PolyStripsType::const_iterator itr2;
    StripsType::const_iterator itr3, itr4;

    for (itr2 = polyStrips.begin(); itr2 != polyStrips.end(); ++itr2) {
        const PStrips &pStrips = itr2->second;

        const StripsType &strips = pStrips.strips;
        const StripPtsType &pts = pStrips.pts;

        const IdsType &poly = pStrips.poly;

        auto Coord = [&poly](const StripPt &a, const StripPt &b) -> double {
            double t = 1-a.t+b.t;

            auto iA = std::find(poly.begin(), poly.end(), a.edge[0]),
                iB = std::find(poly.begin(), poly.end(), b.edge[0]);

            if (iA == iB && a.t < b.t) {
                return b.t-a.t;
            }

            for (;;) {
                if (++iA == poly.end()) {
                    iA = poly.begin();
                }

                if (iA == iB) {
                    break;
                }

                t += 1;
            }

            return t;
        };

        if (!strips.empty()) {
            for (itr3 = strips.begin(); itr3 != strips.end()-1; ++itr3) {
                const StripType &stripA = *itr3;

                const StripPt &pA = pts.at(stripA.front().ind),
                    &pB = pts.at(stripA.back().ind);

                if (pA.capt != CAPT_NOT && pB.capt != CAPT_NOT) {
                    double c1 = Coord(pA, pB);

                    for (itr4 = itr3+1; itr4 != strips.end(); ++itr4) {
                        const StripType &stripB = *itr4;

                        const StripPt &pC = pts.at(stripB.front().ind),
                            &pD = pts.at(stripB.back().ind);

                        if (pC.capt != CAPT_NOT && pD.capt != CAPT_NOT) {

                            if (pC.ind != pA.ind && pC.ind != pB.ind
                                && pD.ind != pA.ind && pD.ind != pB.ind) {

                                double c2 = Coord(pA, pC),
                                    c3 = Coord(pA, pD);

                                if ((c2 < c1) != (c3 < c1)) {
#ifdef DEBUG
                                    std::cout << "c1=" << c1 << ", c2=" << c2 << ", c3=" << c3
                                        << ", poly=" << itr2->first
                                        << ", " << ((pd == modPdA) ? "A" : "B")
                                        << std::endl;
#endif

                                    return true;
                                }
                            }
                        }

                    }

                }

            }
        }

    }
    */

    PolyStripsType::const_iterator itr2;
    StripsType::const_iterator itr3;
    StripType::const_iterator itr4;

    for (itr2 = polyStrips.begin(); itr2 != polyStrips.end(); ++itr2) {
        const PStrips &pStrips = itr2->second;

        const StripsType &strips = pStrips.strips;
        const StripPtsType &pts = pStrips.pts;

        vtkIdList *cell = vtkIdList::New();
        pd->GetCellPoints(itr2->first, cell);

        Base base(pd->GetPoints(), cell);

        cell->Delete();

        AABB tree;

        std::vector<std::shared_ptr<Line>> lines;

        for (itr3 = strips.begin(); itr3 != strips.end(); ++itr3) {
            const StripType &strip = *itr3;

            for (itr4 = strip.begin(); itr4 != strip.end()-1; ++itr4) {
                const StripPt &spA = pts.at(itr4->ind),
                    &spB = pts.at((itr4+1)->ind);

                double ptA[2], ptB[2];

                Transform(spA.pt, ptA, base);
                Transform(spB.pt, ptB, base);

                int grp = itr3-strips.begin();

                lines.emplace_back(new Line({ptA, itr4->ind}, {ptB, (itr4+1)->ind}, grp));
            }

        }

        for (auto &line : lines) {
            // std::cout << *line << std::endl;

            tree.InsertObj(line);
        }

        Bnds bnds(-E, E, -E, E);

        std::vector<std::shared_ptr<Line>>::const_iterator itr5;
        std::vector<std::shared_ptr<Obj>>::const_iterator itr6;

        for (itr5 = lines.begin(); itr5 != lines.end(); ++itr5) {
            auto found = tree.Search(*itr5);

            const Line &lA = **itr5;

            for (itr6 = found.begin(); itr6 != found.end(); ++itr6) {
                // std::cout << itr6->use_count() << std::endl;

                const Line &lB = dynamic_cast<Line&>(**itr6);

                // die linien dürfen nicht zum gleichen strip gehören und sich nicht an den enden berühren

                if (lA.grp != lB.grp
                    && lA.pA.id != lB.pA.id
                    && lA.pA.id != lB.pB.id
                    && lA.pB.id != lB.pA.id
                    && lA.pB.id != lB.pB.id
                    && Intersect2(lA.pA.pt, lA.pB.pt, lB.pA.pt, lB.pB.pt, bnds)) {
                    return true;
                }
            }
        }

    }

    return false;

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
                StripType s(itr->rbegin(), itr->rend()-1);
                itr->insert(itr->begin(), s.begin(), s.end());

            } else if (end.capt == CAPT_NOT) {
                StripType s(itr->rbegin()+1, itr->rend());
                itr->insert(itr->end(), s.begin(), s.end());

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
#ifdef DEBUG
    std::cout << "CutCells()" << std::endl;
#endif

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

        int origId = origCellIds->GetValue(polyInd);

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

        // holes sammeln

        HolesType holes;

        for (itr2 = strips.begin(); itr2 != strips.end(); ++itr2) {
            StripType &s = *itr2;
            if (pts[s.front().ind].capt == CAPT_NOT && pts[s.back().ind].capt == CAPT_NOT) {
                IdsType hole;

                for (auto& sp : s) {
                    hole.push_back(pdPts->InsertNextPoint(pts[sp.ind].pt));
                }

                // anfang und ende sind ja gleich
                hole.pop_back();

                holes.push_back(std::move(hole));

            }
        }

        // holes löschen
        strips.erase(std::remove_if(strips.begin(), strips.end(), [&](const StripType &s) {
            return pts[s.front().ind].capt == CAPT_NOT && pts[s.back().ind].capt == CAPT_NOT; }), strips.end());

        for (itr2 = strips.begin(); itr2 != strips.end(); ++itr2) {
            StripType &strip = *itr2;

#ifdef DEBUG
            std::cout << "strip [";
            for (auto &s : strip) {
                std::cout << s.ind << ", ";
            }
            std::cout << "]" << std::endl;
#endif

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

                            double dA = Mod(iA-i, numPts)+eA_.t,
                                dB = Mod(iB-i, numPts)+eB_.t;

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

                            vtkPoints *pts_ = vtkPoints::New();
                            pts_->SetNumberOfPoints(num);

                            for (itr7 = poly_.begin(); itr7 != poly_.end(); ++itr7) {
                                StripPtR &sp = *itr7;
                                int i = itr7-poly_.begin();

                                pts_->SetPoint(i, pts[sp.ind].cutPt);
                            }

                            double n[3];
                            ComputeNormal(pts_, n);

                            pts_->Delete();

                            double ang = vtkMath::Dot(pStrips.n, n);

#ifdef DEBUG
                            std::cout << "ang=" << ang*180/PI << std::endl;
#endif

                            return ang < .999999;

                        }
                    } else {
                        // gleicher strip

                        StripType &strip = strips[a.strip];

                        if (HasArea(strip)) {
                            RefsType poly_(strip.begin(), strip.end()-1);

                            int num = poly_.size();

                            vtkPoints *pts_ = vtkPoints::New();
                            pts_->SetNumberOfPoints(num);

                            for (itr7 = poly_.begin(); itr7 != poly_.end(); ++itr7) {
                                StripPtR &sp = *itr7;
                                int i = itr7-poly_.begin();

                                pts_->SetPoint(i, pts[sp.ind].cutPt);
                            }

                            double n[3];
                            ComputeNormal(pts_, n);

                            pts_->Delete();

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

                        // erstellt die history

                        auto _s = std::find_if(edge.begin(), edge.end(), [&start](const StripPtR &r) {
                            return &start == &r;
                        });

                        if (_s != edge.end()) {
                            for (itr7 = edge.begin(); itr7 != edge.end(); ++itr7) {
                                StripPtR &sp = *itr7;
                                StripPt &_sp = pts[sp.ind];
                                auto &history = _sp.history;

                                if (&sp != &start) {
                                    if (_sp.t < pts[start.ind].t) {
                                        history.push_back({history.back().f, start.desc[0]});
                                    } else {
                                        history.push_back({start.desc[1], history.back().g});
                                    }

                                }

                            }
                        }

                        auto _e = std::find_if(edge.begin(), edge.end(), [&end](const StripPtR &r) {
                            return &end == &r;
                        });

                        if (_e != edge.end()) {
                            for (itr7 = edge.begin(); itr7 != edge.end(); ++itr7) {
                                StripPtR &sp = *itr7;
                                StripPt &_sp = pts[sp.ind];
                                auto &history = _sp.history;

                                if (&sp != &end) {
                                    if (_sp.t < pts[end.ind].t) {
                                        history.push_back({history.back().f, end.desc[1]});
                                    } else {
                                        history.push_back({end.desc[0], history.back().g});
                                    }

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

                            double d = GetD(ptA, ptB);

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

        IdsType descIds;

        for (itr9 = polys.begin(); itr9 != polys.end(); ++itr9) {
            IdsType &p = *itr9;

            int num = p.size();

            vtkIdList *cell = vtkIdList::New();
            cell->SetNumberOfIds(num);

            for (int i = 0; i < num; i++) {
                cell->SetId(i, p[i]);
            }

            descIds.push_back(pd->InsertNextCell(VTK_POLYGON, cell));

            origCellIds->InsertNextValue(origId);

            cell->Delete();
        }

        pd->DeleteCell(polyInd);

        // holes verarbeiten

        if (!holes.empty()) {
            _Wrapper w(pd, descIds, origId);

            for (auto& hole : holes) {
                w.Add(hole);
            }

            w.MergeAll();
        }

    }

    pd->RemoveDeletedCells();
    pd->BuildCells();

}

void vtkPolyDataBooleanFilter::CollapseCaptPoints (vtkPolyData *pd, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "CollapseCaptPoints()" << std::endl;
#endif

    struct Cmp {
        bool operator() (const StripPt &l, const StripPt &r) const {
            const int x1 = static_cast<int>(l.cutPt[0]*1e5),
                y1 = static_cast<int>(l.cutPt[1]*1e5),
                z1 = static_cast<int>(l.cutPt[2]*1e5),
                x2 = static_cast<int>(r.cutPt[0]*1e5),
                y2 = static_cast<int>(r.cutPt[1]*1e5),
                z2 = static_cast<int>(r.cutPt[2]*1e5);

            return std::tie(x1, y1, z1) < std::tie(x2, y2, z2);
        }
    };

    struct Cmp2 {
        bool operator() (const StripPt &l, const StripPt &r) const {
            return l.ind < r.ind;
        }
    };

    std::map<StripPt, std::set<std::reference_wrapper<StripPt>, Cmp2>, Cmp> test;

    PolyStripsType::iterator itr;
    StripPtsType::iterator itr2;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
        PStrips &pStrips = itr->second;

        for (itr2 = pStrips.pts.begin(); itr2 != pStrips.pts.end(); ++itr2) {
            StripPt &sp = itr2->second;

            test[sp].insert(sp);
        }
    }

    vtkIdList *cells = vtkIdList::New();

    for (auto &s : test) {
        auto &pts = s.second;

        if (pts.size() > 1) {
            StripPt &_a = *(pts.begin()),
                &_b = *(std::next(pts.begin()));

            std::reference_wrapper<StripPt> __a(_a), __b(_b);

            if (_a.capt == CAPT_EDGE && _b.capt == CAPT_A) {
                // so ist a niemals CAPT_EDGE, wenn b CAPT_A ist
                __a = _b;
                __b = _a;
            }

            StripPt &a = __a,
                &b = __b;

            int indA = a.ind,
                indB = b.ind;

#ifdef DEBUG
            std::cout << "collapsing " << b.ind << " -> " << a.ind << std::endl;
#endif

            contLines->GetPointCells(b.ind, cells);

            for (int i = 0; i < cells->GetNumberOfIds(); i++) {
                contLines->ReplaceCellPoint(cells->GetId(i), indB, indA);

                // aktualisiert die links
                contLines->RemoveReferenceToCell(indB, cells->GetId(i));
                contLines->ResizeCellList(indA, 1);
                contLines->AddReferenceToCell(indA, cells->GetId(i));
            }

            std::set<Pair> pairs;

            if (a.capt == CAPT_A && b.capt == CAPT_A) {
                std::map<int, IdsType> shared;

                vtkIdList *lines = vtkIdList::New(),
                    *cell = vtkIdList::New();
                contLines->GetPointCells(indA, lines);

                for (int i = 0; i < lines->GetNumberOfIds(); i++) {
                    contLines->GetCellPoints(lines->GetId(i), cell);

                    int pA = cell->GetId(0),
                        pB = cell->GetId(1);

                    shared[pA == indA ? pB : pA].push_back(lines->GetId(i));
                }

                cell->Delete();
                lines->Delete();

                for (auto &s : shared) {
                    if (s.second.size() > 1) {
                        assert(s.second.size() == 2);

                        for (int l : s.second) {
                            contLines->DeleteCell(l);
                            contLines->RemoveCellReference(l);

                            pairs.insert({indA, s.first});
                            pairs.insert({s.first, indA});
                        }
                    }
                }

            }

            if (a.capt == CAPT_A && b.capt == CAPT_EDGE) {
                vtkIdList *lines = vtkIdList::New(),
                    *cell = vtkIdList::New();
                contLines->GetPointCells(indA, lines);

                for (int i = 0; i < lines->GetNumberOfIds(); i++) {
                    if (contLines->GetCellType(lines->GetId(i)) != VTK_EMPTY_CELL) {
                        contLines->GetCellPoints(lines->GetId(i), cell);

                        int pA = cell->GetId(0),
                            pB = cell->GetId(1);

                        if (pA == indA && pB == indA) {
                            contLines->DeleteCell(lines->GetId(i));
                            contLines->RemoveCellReference(lines->GetId(i));
                        }
                    }

                }

                cell->Delete();
                lines->Delete();
            }

            auto Fct = [&](PolyStripsType &polyStrips_) -> void {
                for (auto &ps : polyStrips_) {
                    StripsType &strips = ps.second.strips;
                    StripPtsType &pts = ps.second.pts;

                    if (pts.count(indB) == 1) {
                        if (pts.count(indA) == 0) {
                            if (a.capt == pts[indB].capt) {
                                pts[indA] = pts[indB];
                                pts[indA].ind = indA;
                                Cpy(pts[indA].pt, a.pt, 3);
                            } else /*if (a.capt == CAPT_A && pts[indB].capt == CAPT_EDGE)*/ {
                                pts[indA] = a;
                            }
                        }

                        for (auto &strip : strips) {
                            if (strip.front().ind == indB) {
                                strip.front().ind = indA;
                            }

                            if (strip.back().ind == indB) {
                                strip.back().ind = indA;
                            }

                        }

                        pts.erase(indB);
                    }

                    if (pairs.size() > 0) {
                        strips.erase(std::remove_if(strips.begin(), strips.end(), [&pairs](const StripType &strip) {
                            return pairs.count({strip.front().ind, strip.back().ind}) == 1;
                        }), strips.end());
                    }
                }
            };

            Fct(polyStripsA);
            Fct(polyStripsB);

        }
    }

    cells->Delete();

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
        FindPoints(loc, itr3->cutPt, pts);
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
        FindPoints(loc, itr3->pt, pts);
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

    vtkKdTreePointLocator *loc = vtkKdTreePointLocator::New();
    loc->SetDataSet(pd);
    loc->BuildLocator();

    vtkIdList *ptsA = vtkIdList::New();
    vtkIdList *ptsB = vtkIdList::New();

    vtkIdList *cells = vtkIdList::New();

    vtkIdList *links = vtkIdList::New();

    typedef std::map<int, int> CountsType;
    std::map<Pair, CountsType> skipped;

    std::vector<StripPtL2>::const_iterator itr3;

    for (itr3 = ends.begin(); itr3 != ends.end(); ++itr3) {
        contLines->GetPointCells(itr3->ind, links);

        if (links->GetNumberOfIds() == 2
            && conts->GetValue(links->GetId(0)) != conts->GetValue(links->GetId(1))) {

            // kein berührender schnitt

            auto &history = itr3->history;

#ifdef DEBUG
            std::cout << "ind: " << itr3->ind << ", history=[";

            for (auto& h : history) {
                std::cout << h << ", ";
            }

            std::cout << "]" << std::endl;
#endif

            std::vector<Pair>::const_reverse_iterator itr4;

            for (itr4 = history.rbegin(); itr4 != history.rend(); ++itr4) {

                double ptA[3], ptB[3];
                pd->GetPoint(itr4->f, ptA);
                pd->GetPoint(itr4->g, ptB);

                FindPoints(loc, ptA, ptsA);
                FindPoints(loc, ptB, ptsB);

                int numPtsA = ptsA->GetNumberOfIds();
                int numPtsB = ptsB->GetNumberOfIds();

                std::vector<Pair> cellsA, cellsB;

                cells->Reset();

                for (int i = 0; i < numPtsA; i++) {
                    pd->GetPointCells(ptsA->GetId(i), cells);
                    for (int j = 0; j < cells->GetNumberOfIds(); j++) {
                        cellsA.push_back({static_cast<int>(ptsA->GetId(i)), static_cast<int>(cells->GetId(j))});
                    }
                }

                cells->Reset();

                for (int i = 0; i < numPtsB; i++) {
                    pd->GetPointCells(ptsB->GetId(i), cells);
                    for (int j = 0; j < cells->GetNumberOfIds(); j++) {
                        cellsB.push_back({static_cast<int>(ptsB->GetId(i)), static_cast<int>(cells->GetId(j))});
                    }
                }

                for (auto &a : cellsA) {
                    for (auto &b : cellsB) {
                        if (a.g == b.g) {
                            // kante existiert noch

#ifdef DEBUG
                            std::cout << "poly: " << a.g
                                << ", edge: (" << a.f << ", " << b.f << ")"
                                << std::endl;
#endif

                            vtkIdList *poly = vtkIdList::New();
                            pd->GetCellPoints(a.g, poly);

                            int numPts = poly->GetNumberOfIds();

                            for (int k = 0; k < numPts; k++) {
                                if (poly->GetId(k) == b.f
                                    && poly->GetId((k+1)%numPts) == a.f) {

                                    CountsType &c = skipped[Pair(itr3->ind, a.g)];

                                    c[b.f]++;
                                    c[a.f]++;
                                }
                            }

                            poly->Delete();

                        }
                    }
                }



            }

        }

        links->Reset();
    }

    links->Delete();
    cells->Delete();

    ptsB->Delete();
    ptsA->Delete();

    loc->FreeSearchStructure();
    loc->Delete();

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

    vtkIntArray *origCellIds = vtkIntArray::SafeDownCast(pd->GetCellData()->GetScalars("OrigCellIds"));

    struct Cmp {
        bool operator() (const StripPtL3 &l, const StripPtL3 &r) const {
            return l.ind < r.ind;
        }
    };

    typedef std::set<StripPtL3, Cmp> AType;
    typedef std::map<Pair, AType> BType;
    typedef std::vector<StripPtL3> CType;

    BType edges;

    PolyStripsType::iterator itr;
    StripPtsType::iterator itr2;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
        PStrips &pStrips = itr->second;

        for (itr2 = pStrips.pts.begin(); itr2 != pStrips.pts.end(); ++itr2) {
            StripPt &sp = itr2->second;

            if (sp.capt == CAPT_EDGE) {
                edges[Pair(sp.edge[0], sp.edge[1])].insert(StripPtL3(sp.pt, sp.t, sp.ind));
            }
        }
    }

    BType::iterator itr4;
    CType::iterator itr5;

    pd->BuildLinks();

    vtkKdTreePointLocator *loc = vtkKdTreePointLocator::New();
    loc->SetDataSet(pd);
    loc->BuildLocator();

    typedef std::vector<Pair> DType;

    IdsType::const_iterator itr6;
    DType::const_iterator itr7, itr8;

    vtkIdList *cells = vtkIdList::New();

    for (itr4 = edges.begin(); itr4 != edges.end(); ++itr4) {
        const Pair &pair = itr4->first;
        const AType &ends = itr4->second;

        CType pts(ends.begin(), ends.end());

        double ptA[3], ptB[3];

        pd->GetPoint(pair.f, ptA);
        pd->GetPoint(pair.g, ptB);

        pts.push_back(StripPtL3(ptA, 0.));
        pts.push_back(StripPtL3(ptB, 1.));

        std::sort(pts.rbegin(), pts.rend());

#ifdef DEBUG
        std::cout << "edge=" << pair << std::endl;

        for (itr5 = pts.begin(); itr5 != pts.end(); ++itr5) {
            std::cout << *itr5 << std::endl;
        }
#endif

        IdsType voids;

        voids.push_back(0);

        for (itr5 = pts.begin()+1; itr5 != pts.end()-1; ++itr5) {
            contLines->GetPointCells(itr5->ind, cells);
            int numCells = cells->GetNumberOfIds();

            std::set<int> seeds;

            for (int i = 0; i < numCells; i++) {
                seeds.insert(conts->GetValue(cells->GetId(i)));
            }

#ifdef DEBUG
            std::cout << itr5->ind << " -> seeds=[";
            for (auto s : seeds) {
                std::cout << s << ", ";
            }
            std::cout << "]" << std::endl;
#endif

            if (seeds.size() != 1) {
                voids.push_back(itr5-pts.begin());
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
            CType pts_(pts.begin()+(*itr6), pts.begin()+(*(itr6+1))+1);

            if (pts_.size() > 2) {

                vtkIdList *ptsA = vtkIdList::New();
                vtkIdList *ptsB = vtkIdList::New();

                FindPoints(loc, pts_.front().pt, ptsA);
                FindPoints(loc, pts_.back().pt, ptsB);

                int numPtsA = ptsA->GetNumberOfIds(),
                    numPtsB = ptsB->GetNumberOfIds();

                DType polysA, polysB;

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

                                    for (itr5 = pts_.begin()+1; itr5 != pts_.end()-1; ++itr5) {
                                        poly_->InsertNextId(pd->InsertNextLinkedPoint(itr5->pt, 1));
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

            FindPoints(loc, pStrips.pts[(strip.begin()+1)->ind].pt, pts);
            int numPts = pts->GetNumberOfIds();

            for (int i = 0; i < numPts; i++) {
                inds[s.ind].insert(pts->GetId(i));
            }

            FindPoints(loc, pStrips.pts[(strip.end()-2)->ind].pt, pts);
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

        FindPoints(loc, pt, pts);
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
                if (GetD(itr4->pt, itr5->pt) < 1e-6) {
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


class PolyAtEdge {
    vtkPolyData *pd;

public:
    PolyAtEdge (vtkPolyData *_pd, int _polyId, int _ptIdA, int _ptIdB) : pd(_pd), polyId(_polyId), ptIdA(_ptIdA), ptIdB(_ptIdB), loc(LOC_NONE) {
        vtkIdList *poly = vtkIdList::New();
        pd->GetCellPoints(polyId, poly);

        double ptA[3], ptB[3];

        pd->GetPoint(ptIdA, ptA);
        pd->GetPoint(ptIdB, ptB);

        vtkMath::Subtract(ptB, ptA, e);
        vtkMath::Normalize(e);

        ComputeNormal(pd->GetPoints(), n, poly);

        vtkMath::Cross(e, n, r);

        poly->Delete();

    }

    int polyId, ptIdA, ptIdB;
    double n[3], e[3], r[3];

    int loc;

    friend std::ostream& operator<< (std::ostream &out, const PolyAtEdge &p) {
        out << "polyId " << p.polyId << ", ptIdA " << p.ptIdA << ", ptIdB " << p.ptIdB;
        return out;
    }

};


class PolyPair {
public:
    double alpha;

    PolyPair (PolyAtEdge _pA, PolyAtEdge _pB) : pA(_pA), pB(_pB) {
        alpha = GetAngle(pA.r, pB.r, pA.e);
    }

    PolyAtEdge pA, pB;

    void GetLoc (PolyAtEdge &pT, int mode) {
        double beta = GetAngle(pA.r, pT.r, pA.e);

#ifdef DEBUG
        std::cout << "GetLoc() -> polyId "
                  << pT.polyId << ", beta " << (beta*180/PI) << std::endl;
#endif

        if (beta < 1e-7 || beta > 2*PI-1e-7) {
            // konkruent gegenüber dem polygon hinter pA

            double o = vtkMath::Dot(pA.n, pT.n);

            if (o < .999999) {
                // normalen sind entgegengesetzt gerichtet

                if (mode == OPER_INTERSECTION) {
                    pA.loc = LOC_OUTSIDE;
                    pT.loc = LOC_OUTSIDE;
                } else {
                    pA.loc = LOC_INSIDE;
                    pT.loc = LOC_INSIDE;
                }
            } else if (mode == OPER_UNION || mode == OPER_INTERSECTION) {
                pA.loc = LOC_INSIDE;
                pT.loc = LOC_OUTSIDE;
            }

        } else if (std::abs(beta-alpha) < 1e-7) {

            double o = vtkMath::Dot(pB.n, pT.n);

            if (o < .999999) {
                // normalen sind entgegengesetzt gerichtet

                if (mode == OPER_INTERSECTION) {
                    pB.loc = LOC_OUTSIDE;
                    pT.loc = LOC_OUTSIDE;
                } else {
                    pB.loc = LOC_INSIDE;
                    pT.loc = LOC_INSIDE;
                }
            } else if (mode == OPER_UNION || mode == OPER_INTERSECTION) {
                pB.loc = LOC_INSIDE;
                pT.loc = LOC_OUTSIDE;
            }

        } else if (beta > alpha) {
            pT.loc = LOC_INSIDE;
        } else {
            pT.loc = LOC_OUTSIDE;
        }
    }

};


PolyPair GetEdgePolys (vtkPolyData *pd, vtkIdList *ptsA, vtkIdList *ptsB) {

#ifdef DEBUG
    std::cout << "GetEdgePolys()" << std::endl;
#endif

    std::vector<Pair> p;

    int numPtsA = ptsA->GetNumberOfIds(),
        numPtsB = ptsB->GetNumberOfIds();

    vtkIdList *polys = vtkIdList::New();

    for (int i = 0; i < numPtsA; i++) {
        polys->Reset();

        pd->GetPointCells(ptsA->GetId(i), polys);
        int numCells = polys->GetNumberOfIds();

        for (int j = 0; j < numCells; j++) {
            p.push_back(Pair(ptsA->GetId(i), polys->GetId(j)));
        }
    }

    for (int i = 0; i < numPtsB; i++) {
        polys->Reset();

        pd->GetPointCells(ptsB->GetId(i), polys);
        int numCells = polys->GetNumberOfIds();

        for (int j = 0; j < numCells; j++) {
            p.push_back(Pair(ptsB->GetId(i), polys->GetId(j)));
        }
    }

    polys->Delete();

    std::map<int, IdsType> pEdges;

    std::vector<Pair>::const_iterator itr;
    for (itr = p.begin(); itr != p.end(); ++itr) {
        pEdges[itr->g].push_back(itr->f);
    }

    std::vector<PolyAtEdge> opp;

    std::map<int, IdsType>::const_iterator itr2;

    for (itr2 = pEdges.begin(); itr2 != pEdges.end(); ++itr2) {
        const std::vector<int> &pts = itr2->second;

        if (pts.size() > 1) {
            vtkIdList *poly = vtkIdList::New();
            pd->GetCellPoints(itr2->first, poly);

            int l = poly->GetNumberOfIds();
            for (int i = 0; i < l; i++) {
                int a = poly->GetId(i),
                    b = poly->GetId((i+1)%l);

                if (std::find(pts.begin(), pts.end(), a) != pts.end()
                    && std::find(pts.begin(), pts.end(), b) != pts.end()) {
                    //std::cout << "(" << a << ", " << b << ") in " << itr2->first << std::endl;

                    opp.push_back(PolyAtEdge(pd, itr2->first, a, b));
                }
            }

            poly->Delete();
        }
    }

    assert(opp.size() == 2);

    PolyPair pp(opp[0], opp[1]);

#ifdef DEBUG
    std::cout << pp.pA << std::endl;
    std::cout << pp.pB << std::endl;

    std::cout << "alpha " << (pp.alpha*180/PI) << std::endl;
#endif

    return pp;

}


void vtkPolyDataBooleanFilter::CombineRegions () {

#ifdef DEBUG
    std::cout << "CombineRegions()" << std::endl;
#endif

    vtkPolyData *filterdA = vtkPolyData::New();
    filterdA->DeepCopy(modPdA);

    vtkPolyData *filterdB = vtkPolyData::New();
    filterdB->DeepCopy(modPdB);

    auto FilterCells = [&](vtkPolyData *pd, RelationsType &rels) -> void {
        RelationsType::const_iterator itr;
        for (itr = rels.begin(); itr != rels.end(); ++itr) {
            if (itr->second == (DecPolys ? Rel::ORIG : Rel::DEC)) {
                pd->DeleteCell(itr->first);
            }
        }

        pd->RemoveDeletedCells();
    };

    FilterCells(filterdA, relsA);
    FilterCells(filterdB, relsB);

    // ungenutzte punkte löschen
    vtkCleanPolyData *cleanA = vtkCleanPolyData::New();
    cleanA->PointMergingOff();
    cleanA->SetInputData(filterdA);

    vtkCleanPolyData *cleanB = vtkCleanPolyData::New();
    cleanB->PointMergingOff();
    cleanB->SetInputData(filterdB);

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

#ifdef DEBUG
    std::cout << "Exporting modPdA_9.vtk" << std::endl;
    WriteVTK("modPdA_9.vtk", cfA->GetOutput());

    std::cout << "Exporting modPdB_9.vtk" << std::endl;
    WriteVTK("modPdB_9.vtk", cfB->GetOutput());
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

        if (contLines->GetCellType(i) == VTK_EMPTY_CELL) {
            continue;
        }

        contLines->GetCellPoints(i, line);

        contLines->GetPoint(line->GetId(0), ptA);
        contLines->GetPoint(line->GetId(1), ptB);

        FindPoints(plA, ptA, fptsA);
        FindPoints(plB, ptA, fptsB);

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

        FindPoints(plA, ptB, lptsA);
        FindPoints(plB, ptB, lptsB);

        PolyPair ppA = GetEdgePolys(pdA, fptsA, lptsA);
        PolyPair ppB = GetEdgePolys(pdB, fptsB, lptsB);

        ppB.GetLoc(ppA.pA, OperMode);
        ppB.GetLoc(ppA.pB, OperMode);

        ppA.GetLoc(ppB.pA, OperMode);
        ppA.GetLoc(ppB.pB, OperMode);

        int fsA = scalarsA->GetTuple1(ppA.pA.ptIdA);
        int lsA = scalarsA->GetTuple1(ppA.pB.ptIdA);

        int fsB = scalarsB->GetTuple1(ppB.pA.ptIdA);
        int lsB = scalarsB->GetTuple1(ppB.pB.ptIdA);

#ifdef DEBUG
        std::cout << "polyId " << ppA.pA.polyId << ", sA " << fsA << ", loc " << ppA.pA.loc << std::endl;
        std::cout << "polyId " << ppA.pB.polyId << ", sA " << lsA << ", loc " << ppA.pB.loc << std::endl;
        std::cout << "polyId " << ppB.pA.polyId << ", sB " << fsB << ", loc " << ppB.pA.loc << std::endl;
        std::cout << "polyId " << ppB.pB.polyId << ", sB " << lsB << ", loc " << ppB.pB.loc << std::endl;

        if (locsA.count(fsA) > 0 && locsA[fsA] != ppA.pA.loc) {
            std::cout << "sA " << fsA << ": " << locsA[fsA] << " -> " << ppA.pA.loc << std::endl;
        }

        if (locsA.count(lsA) > 0 && locsA[lsA] != ppA.pB.loc) {
            std::cout << "sA " << lsA << ": " << locsA[lsA] << " -> " << ppA.pB.loc << std::endl;
        }

        if (locsB.count(fsB) > 0 && locsB[fsB] != ppB.pA.loc) {
            std::cout << "sB " << fsB << ": " << locsB[fsB] << " -> " << ppB.pA.loc << std::endl;
        }

        if (locsB.count(lsB) > 0 && locsB[lsB] != ppB.pB.loc) {
            std::cout << "sB " << lsB << ": " << locsB[lsB] << " -> " << ppB.pB.loc << std::endl;
        }

#endif

        locsA[fsA] = ppA.pA.loc;
        locsA[lsA] = ppA.pB.loc;

        locsB[fsB] = ppB.pA.loc;
        locsB[lsB] = ppB.pB.loc;

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
        int cellId = cellIdsA->GetValue(origCellIdsA->GetValue(i));

        newOrigCellIdsA->InsertNextValue(cellId);
        newOrigCellIdsB->InsertNextValue(-1);

        newCellDataA->CopyData(cellDataA, cellId, i);
    }

    for (int i = 0; i < regsB->GetNumberOfCells(); i++) {
        int cellId = cellIdsB->GetValue(origCellIdsB->GetValue(i));

        newOrigCellIdsB->InsertNextValue(cellId);
        newOrigCellIdsA->InsertNextValue(-1);

        newCellDataB->CopyData(cellDataB, cellId, i);
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

    // resultA ist erster output des filters
    resultA->ShallowCopy(cfPd);

    resultA->GetCellData()->AddArray(newOrigCellIdsA);
    resultA->GetCellData()->AddArray(newOrigCellIdsB);

    resultB->ShallowCopy(contLines);

    // aufräumen

    cfApp->Delete();
    cleanApp->Delete();
    app->Delete();

    newOrigCellIdsB->Delete();
    newOrigCellIdsA->Delete();

    plB->FreeSearchStructure();
    plB->Delete();

    plA->FreeSearchStructure();
    plA->Delete();

    cfB->Delete();
    cfA->Delete();

    cleanB->Delete();
    cleanA->Delete();

    filterdB->Delete();
    filterdA->Delete();

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
    pdA->DeepCopy(modPdA);

    vtkPolyData *pdB = vtkPolyData::New();
    pdB->DeepCopy(modPdB);

    pdA->GetCellData()->AddArray(origCellIdsA);
    pdB->GetCellData()->AddArray(origCellIdsB);

    origCellIdsA->Delete();
    origCellIdsB->Delete();

    auto FilterCells = [&](vtkPolyData *pd, RelationsType &rels) -> void {
        RelationsType::const_iterator itr;
        for (itr = rels.begin(); itr != rels.end(); ++itr) {
            if (itr->second == (DecPolys ? Rel::ORIG : Rel::DEC)) {
                pd->DeleteCell(itr->first);
            }
        }

        pd->RemoveDeletedCells();
    };

    FilterCells(pdA, relsA);
    FilterCells(pdB, relsB);

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
    appPd->GetCellData()->RemoveArray("OrigCellIds");

    vtkCleanPolyData *clean = vtkCleanPolyData::New();
    clean->PointMergingOff();
    clean->SetInputData(appPd);

    clean->Update();

    vtkPolyData *cleanPd = clean->GetOutput();

    resultA->ShallowCopy(cleanPd);
    resultB->ShallowCopy(contLines);

    clean->Delete();
    app->Delete();

    pdB->Delete();
    pdA->Delete();

}

void _Wrapper::MergeAll () {
    vtkPoints *pdPts = pd->GetPoints();

    vtkIntArray *origCellIds = vtkIntArray::SafeDownCast(pd->GetCellData()->GetScalars("OrigCellIds"));

    vtkIdList *cell = vtkIdList::New();

    // descendants in holes einfügen

    std::set<int> outerIds;

    for (auto& id : descIds) {
        pd->GetCellPoints(id, cell);
        IdsType hole;

        for (int i = 0; i < cell->GetNumberOfIds(); i++) {
            hole.push_back(cell->GetId(i));
        }

        outerIds.insert(hole.begin(), hole.end());

        holes.push_back(std::move(hole));

        // löscht
        pd->DeleteCell(id);
    }

    base = Base(pdPts, cell);

    Merger m;

    for (auto& hole : holes) {
        PolyType casted;

        for (int id : hole) {
            double pt[3];
            pd->GetPoint(id, pt);

            double _pt[2];
            Transform(pt, _pt, base);

            casted.push_back({_pt, id});
        }

        if (TestCW(casted)) {
            std::reverse(casted.begin(), casted.end());
        }

        m.AddPoly(casted);
    }

    PolysType merged;
    m.GetMerged(merged);

    std::set<int> usedIds;

    for (auto& poly : merged) {
        cell->Reset();

        // poly ist immer ccw

        std::map<int, int> repl;

        for (auto& p : poly) {
            if (outerIds.count(p.id) == 1 || usedIds.count(p.id) == 0) {
                repl[p.id] = p.id;
            } else {
                double _pt[3];
                //BackTransform(p.pt, _pt, base);

                pd->GetPoint(p.id, _pt);

                repl[p.id] = pdPts->InsertNextPoint(_pt);

            }
        }


        for (auto& p : poly) {
            cell->InsertNextId(repl[p.id]);
            usedIds.insert(p.id);

        }

        pd->InsertNextCell(VTK_POLYGON, cell);
        origCellIds->InsertNextValue(origId);

    }

    cell->Delete();

}

void vtkPolyDataBooleanFilter::DecPolys_ (vtkPolyData *pd, InvolvedType &involved, RelationsType &rels) {

#ifdef DEBUG
    std::cout << "DecPolys_()" << std::endl;
#endif

    if (!DecPolys || !rels.empty()) {
        return;
    }

    vtkPoints *pdPts = pd->GetPoints();

    vtkIntArray *origCellIds = vtkIntArray::SafeDownCast(pd->GetCellData()->GetScalars("OrigCellIds"));

    int numCells = pd->GetNumberOfCells();

    vtkIdList *cells = vtkIdList::New();

    for (int i = 0; i < numCells; i++) {
        if (involved.count(origCellIds->GetTuple1(i)) == 1) {
            cells->InsertNextId(i);
        }
    }

    for (int i = 0; i < cells->GetNumberOfIds(); i++) {

        int cellId = cells->GetId(i),
            origId = origCellIds->GetValue(cellId);

#ifdef DEBUG
        std::cout << "cellId " << cellId << std::endl;
#endif

        // if (cellId != 12) {
        //     continue;
        // }

        auto _ps = (pd == modPdA ? polyStripsA : polyStripsB).at(origId);

        vtkIdList *cell = vtkIdList::New();

        pd->GetCellPoints(cellId, cell);

        Base base(pdPts, cell);

        int numPts = cell->GetNumberOfIds();

        if (numPts > 3) {

            IdsType ptIds;

            for (int k = 0; k < numPts; k++) {
                ptIds.push_back(cell->GetId(k));
            }

            std::reverse(ptIds.begin(), ptIds.end());

            PolyType poly;

            int j = 0;

            for (int id : ptIds) {
                double pt[3],
                    _pt[2];

                pd->GetPoint(id, pt);
                Transform(pt, _pt, base);

                poly.push_back({_pt, j++});
            }

            assert(TestCW(poly));

            vtkIdList *newCell = vtkIdList::New();

            try {

                Decomposer d(poly);

                DecResType decs;
                d.GetDecomposed(decs);

                for (auto& dec : decs) {
                    newCell->Reset();

                    std::reverse(dec.begin(), dec.end());

                    for (int id : dec) {
                        newCell->InsertNextId(ptIds[id]);
                    }

                    int newId = pd->InsertNextCell(VTK_POLYGON, newCell);
                    origCellIds->InsertNextValue(origId);

                    rels[newId] = Rel::DEC;

                }

                rels[cellId] = Rel::ORIG;

                // std::raise(SIGSEGV);

            } catch (const std::exception &e) {
                std::stringstream ss;
                ss << "Exception on " << GetAbsolutePath(poly);

                std::cerr << ss.str()
                    << ", " << e.what() << std::endl;
            }

            newCell->Delete();

        }

        cell->Delete();
    }

    cells->Delete();

}
