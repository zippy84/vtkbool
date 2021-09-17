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

#include <map>
#include <deque>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <functional>
#include <queue>
#include <memory>

#include <chrono>
#include <numeric>
#include <iterator>

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

vtkStandardNewMacro(vtkPolyDataBooleanFilter);

vtkPolyDataBooleanFilter::vtkPolyDataBooleanFilter () {

    SetNumberOfInputPorts(2);
    SetNumberOfOutputPorts(3);

    timePdA = 0;
    timePdB = 0;

    contLines = vtkSmartPointer<vtkPolyData>::New();

    modPdA = vtkSmartPointer<vtkPolyData>::New();
    modPdB = vtkSmartPointer<vtkPolyData>::New();

    cellDataA = vtkSmartPointer<vtkCellData>::New();
    cellDataB = vtkSmartPointer<vtkCellData>::New();

    cellIdsA = vtkSmartPointer<vtkIdTypeArray>::New();
    cellIdsB = vtkSmartPointer<vtkIdTypeArray>::New();

    OperMode = OPER_UNION;

}

vtkPolyDataBooleanFilter::~vtkPolyDataBooleanFilter () {
    // nix mehr
}

int vtkPolyDataBooleanFilter::ProcessRequest(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

    if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA())) {

        vtkInformation *inInfoA = inputVector[0]->GetInformationObject(0);
        vtkInformation *inInfoB = inputVector[1]->GetInformationObject(0);

        vtkPolyData *pdA = vtkPolyData::SafeDownCast(inInfoA->Get(vtkDataObject::DATA_OBJECT()));
        vtkPolyData *pdB = vtkPolyData::SafeDownCast(inInfoB->Get(vtkDataObject::DATA_OBJECT()));

        vtkInformation *outInfoA = outputVector->GetInformationObject(0);
        vtkInformation *outInfoB = outputVector->GetInformationObject(1);
        vtkInformation *outInfoC = outputVector->GetInformationObject(2);

        resultA = vtkPolyData::SafeDownCast(outInfoA->Get(vtkDataObject::DATA_OBJECT()));
        resultB = vtkPolyData::SafeDownCast(outInfoB->Get(vtkDataObject::DATA_OBJECT()));
        resultC = vtkPolyData::SafeDownCast(outInfoC->Get(vtkDataObject::DATA_OBJECT()));

        using clock = std::chrono::steady_clock;
        std::vector<clock::duration> times;
        clock::time_point start;

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

            start = clock::now();

            vtkSmartPointer<vtkPolyDataContactFilter> cl = vtkSmartPointer<vtkPolyDataContactFilter>::New();
            cl->SetInputConnection(0, cleanA->GetOutputPort());
            cl->SetInputConnection(1, cleanB->GetOutputPort());
            cl->Update();

            times.push_back(clock::now()-start);

            contLines->DeepCopy(cl->GetOutput());

            modPdA->DeepCopy(cl->GetOutput(1));
            modPdB->DeepCopy(cl->GetOutput(2));

#ifdef DEBUG
            std::cout << "Exporting contLines.vtk" << std::endl;
            WriteVTK("contLines.vtk", contLines);

            std::cout << "Exporting modPdA_1.vtk" << std::endl;
            WriteVTK("modPdA_1.vtk", modPdA);

            std::cout << "Exporting modPdB_1.vtk" << std::endl;
            WriteVTK("modPdB_1.vtk", modPdB);
#endif

            if (contLines->GetNumberOfCells() == 0) {
                vtkErrorMacro("Inputs have no contact.");

                return 1;
            }

            // in den CellDatas steht drin, welche polygone einander schneiden

            vtkIdTypeArray *contsA = vtkIdTypeArray::SafeDownCast(contLines->GetCellData()->GetScalars("cA"));
            vtkIdTypeArray *contsB = vtkIdTypeArray::SafeDownCast(contLines->GetCellData()->GetScalars("cB"));

            vtkIdTypeArray *sourcesA = vtkIdTypeArray::SafeDownCast(contLines->GetCellData()->GetScalars("sourcesA"));
            vtkIdTypeArray *sourcesB = vtkIdTypeArray::SafeDownCast(contLines->GetCellData()->GetScalars("sourcesB"));

            vtkIdType i, numPts = contLines->GetNumberOfPoints();

            vtkIdList *cells = vtkIdList::New();

            for (i = 0; i < numPts; i++) {
                contLines->GetPointCells(i, cells);

                if (cells->GetNumberOfIds() == 1) {
                    std::cout << "Contact ends at " << i << " (line " << cells->GetId(0) << ")" << std::endl;

                    break;
                }
            }

            cells->Delete();

            if (i < numPts) {
                vtkErrorMacro("Contact ends suddenly at point " << i << ".");

                return 1;
            }

            // sichert die OrigCellIds

            vtkIdTypeArray *origCellIdsA = vtkIdTypeArray::SafeDownCast(modPdA->GetCellData()->GetScalars("OrigCellIds"));
            vtkIdTypeArray *origCellIdsB = vtkIdTypeArray::SafeDownCast(modPdB->GetCellData()->GetScalars("OrigCellIds"));

            cellIdsA->DeepCopy(origCellIdsA);
            cellIdsB->DeepCopy(origCellIdsB);

            vtkIdType numCellsA = modPdA->GetNumberOfCells(),
                numCellsB = modPdB->GetNumberOfCells();

            for (i = 0; i < numCellsA; i++) {
                origCellIdsA->SetValue(i, i);
            }

            for (i = 0; i < numCellsB; i++) {
                origCellIdsB->SetValue(i, i);
            }

            start = clock::now();

            if (GetPolyStrips(modPdA, contsA, sourcesA, polyStripsA) ||
                GetPolyStrips(modPdB, contsB, sourcesB, polyStripsB)) {

                vtkErrorMacro("Strips are invalid.");

                return 1;

            }

            times.push_back(clock::now()-start);

            // trennt die polygone an den linien

            start = clock::now();

            CutCells(modPdA, polyStripsA);
            CutCells(modPdB, polyStripsB);

            times.push_back(clock::now()-start);

#ifdef DEBUG
            std::cout << "Exporting modPdA_2.vtk" << std::endl;
            WriteVTK("modPdA_2.vtk", modPdA);

            std::cout << "Exporting modPdB_2.vtk" << std::endl;
            WriteVTK("modPdB_2.vtk", modPdB);
#endif

            start = clock::now();

            RestoreOrigPoints(modPdA, polyStripsA);
            RestoreOrigPoints(modPdB, polyStripsB);

            times.push_back(clock::now()-start);

#ifdef DEBUG
            std::cout << "Exporting modPdA_3.vtk" << std::endl;
            WriteVTK("modPdA_3.vtk", modPdA);

            std::cout << "Exporting modPdB_3.vtk" << std::endl;
            WriteVTK("modPdB_3.vtk", modPdB);
#endif

            start = clock::now();

            ResolveOverlaps(modPdA, contsA, polyStripsA);
            ResolveOverlaps(modPdB, contsB, polyStripsB);

            times.push_back(clock::now()-start);

#ifdef DEBUG
            std::cout << "Exporting modPdA_4.vtk" << std::endl;
            WriteVTK("modPdA_4.vtk", modPdA);

            std::cout << "Exporting modPdB_4.vtk" << std::endl;
            WriteVTK("modPdB_4.vtk", modPdB);
#endif

            start = clock::now();

            AddAdjacentPoints(modPdA, contsA, polyStripsA);
            AddAdjacentPoints(modPdB, contsB, polyStripsB);

            times.push_back(clock::now()-start);

#ifdef DEBUG
            std::cout << "Exporting modPdA_5.vtk" << std::endl;
            WriteVTK("modPdA_5.vtk", modPdA);

            std::cout << "Exporting modPdB_5.vtk" << std::endl;
            WriteVTK("modPdB_5.vtk", modPdB);
#endif

            start = clock::now();

            DisjoinPolys(modPdA, polyStripsA);
            DisjoinPolys(modPdB, polyStripsB);

            times.push_back(clock::now()-start);

#ifdef DEBUG
            std::cout << "Exporting modPdA_6.vtk" << std::endl;
            WriteVTK("modPdA_6.vtk", modPdA);

            std::cout << "Exporting modPdB_6.vtk" << std::endl;
            WriteVTK("modPdB_6.vtk", modPdB);
#endif

            start = clock::now();

            MergePoints(modPdA, polyStripsA);
            MergePoints(modPdB, polyStripsB);

            times.push_back(clock::now()-start);

#ifdef DEBUG
            std::cout << "Exporting modPdA_7.vtk" << std::endl;
            WriteVTK("modPdA_7.vtk", modPdA);

            std::cout << "Exporting modPdB_7.vtk" << std::endl;
            WriteVTK("modPdB_7.vtk", modPdB);
#endif

            timePdA = pdA->GetMTime();
            timePdB = pdB->GetMTime();

        }

        start = clock::now();

        if (CombineRegions()) {
            vtkErrorMacro("Boolean operation failed.");

            return 1;
        }

        times.push_back(clock::now()-start);

// #ifdef DEBUG
        double sum = std::chrono::duration_cast<std::chrono::duration<double>>(std::accumulate(times.begin(), times.end(), clock::duration())).count();

        std::vector<clock::duration>::const_iterator itr;
        for (itr = times.begin(); itr != times.end(); itr++) {
            double time = std::chrono::duration_cast<std::chrono::duration<double>>(*itr).count();

            std::cout << "Time " << (itr-times.begin())
                << ": " << time << "s (" << (time/sum*100) << "%)"
                << std::endl;
        }
// #endif

    }

    return 1;

}

void vtkPolyDataBooleanFilter::GetStripPoints (vtkPolyData *pd, vtkIdTypeArray *sources, PStrips &pStrips, IdsType &lines) {

#ifdef DEBUG
    std::cout << "GetStripPoints()" << std::endl;
#endif

    StripPtsType &pts = pStrips.pts;
    const IdsType &poly = pStrips.poly;

    std::size_t i, j, numPts = poly.size();

    double tol = 1e-5;

    vtkIdList *linePts = vtkIdList::New();

    IdsType::const_iterator itr;

    for (itr = lines.begin(); itr != lines.end(); itr++) {
        contLines->GetCellPoints(*itr, linePts);

        // diese punkte durchlaufen

        for (int _i : {0, 1}) {

            vtkIdType realInd = linePts->GetId(_i);

            if (pts.count(realInd) == 0) {
                // lage analysieren

                StripPt sp;

                sp.ind = realInd;

                // die koordinaten
                double pt[3];
                contLines->GetPoint(realInd, pt);

                Cpy(sp.pt, pt, 3);

                vtkIdType src = sources->GetTypedComponent(*itr, _i);

                sp.src = src;

                double lastD = DBL_MAX;

                // jetzt muss man die kanten durchlaufen
                for (i = 0; i < numPts; i++) {
                    j = i+1 == numPts ? 0 : i+1;

                    if (src != NO_USE && !(i == src || j == src)) {
                        continue;
                    }

                    vtkIdType indA, indB;
                    double a[3], b[3];

                    indA = poly[i];
                    indB = poly[j];

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

                        sp.edge[0] = indA;
                        sp.edge[1] = indB;

                        sp.t = std::min(1., std::max(0., t));

                        if (vtkMath::Norm(sA) < tol) {
                            Cpy(sp.captPt, a, 3);
                            sp.capt = Capt::A;

                            break;

                        } else if (vtkMath::Norm(sB) < tol) {
                            Cpy(sp.captPt, b, 3);
                            sp.capt = Capt::B;

                            break;

                        } else {
                            // u ist nicht normiert
                            vtkMath::MultiplyScalar(u, t);

                            double x[3];
                            vtkMath::Add(a, u, x);

                            // projektion
                            Cpy(sp.captPt, x, 3);

                            sp.capt = Capt::EDGE;

                            lastD = d;
                        }

                    }

                }

                if (src != NO_USE && sp.edge[0] == NO_USE) {
                    sp.catched = false;
                }

                pts.emplace(realInd, std::move(sp));

            }
        }
    }

    linePts->Delete();

    StripPtsType::iterator itr2;

    for (itr2 = pts.begin(); itr2 != pts.end(); ++itr2) {
        StripPt &sp = itr2->second;

        if (sp.capt != Capt::NOT) {
            if (sp.capt == Capt::B) {
                sp.t = 0;

                sp.edge[0] = sp.edge[1];

                vtkIdType i = std::find(poly.begin(), poly.end(), sp.edge[0])-poly.begin(); // falscher typ!

                sp.edge[1] = poly[i+1 == numPts ? 0 : i+1];

                sp.capt = Capt::A;

            }

            // für den schnitt werden die eingerasteten koordinaten verwendet

            Cpy(sp.cutPt, sp.captPt, 3);
        } else {

            Cpy(sp.cutPt, sp.pt, 3);
        }

        sp.history.emplace_back(sp.edge[0], sp.edge[1]);

    }

#ifdef DEBUG
    std::cout << "pts: " << std::endl;
    for (itr2 = pts.begin(); itr2 != pts.end(); ++itr2) {
        std::cout << itr2->first << ": " << itr2->second << std::endl;
    }
#endif

}

bool vtkPolyDataBooleanFilter::GetPolyStrips (vtkPolyData *pd, vtkIdTypeArray *conts, vtkIdTypeArray *sources, PolyStripsType &polyStrips) {
#ifdef DEBUG
    std::cout << "GetPolyStrips()" << std::endl;
#endif

    polyStrips.clear();

    std::map<vtkIdType, IdsType> polyLines;

    for (vtkIdType i = 0; i < conts->GetNumberOfTuples(); i++) {
        vtkIdType poly = conts->GetValue(i);

        /*if (poly != 1641) {
            continue;
        }*/

        polyLines[poly].push_back(i);
    }

    std::vector<std::reference_wrapper<StripPt>> notCatched;

    std::map<vtkIdType, IdsType>::iterator itr;

    for (itr = polyLines.begin(); itr != polyLines.end(); ++itr) {

        IdsType &lines = itr->second;
        RemoveDuplicates(lines);

        PStrips &pStrips = polyStrips[itr->first];

        const vtkIdType *polyPts;

        vtkIdType i, numPts;
        pd->GetCellPoints(itr->first, numPts, polyPts);

        for (i = 0; i < numPts; i++) {
            pStrips.poly.push_back(polyPts[i]);
        }

        ComputeNormal(pd->GetPoints(), pStrips.n, numPts, polyPts);

        GetStripPoints(pd, sources, pStrips, lines);

        for (auto &sp : pStrips.pts) {
            sp.second.polyId = itr->first;

            if (!sp.second.catched) {
                notCatched.push_back(sp.second);
            }
        }

    }

    auto Next = [](const IdsType &ids, vtkIdType id) -> vtkIdType {
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
                    if (corr.capt == Capt::A) {
                        sp.capt = Capt::A;
                        sp.edge[0] = corr.edge[0];
                        sp.edge[1] = Next(polyStrips[sp.polyId].poly, sp.edge[0]);

                        sp.t = 0;

                        Cpy(sp.captPt, corr.captPt, 3);
                        Cpy(sp.cutPt, sp.captPt, 3);

                        sp.history.emplace_back(sp.edge[0], sp.edge[1]);

                        sp.catched = true;

                    }
                }
            } catch (...) {}

        }

        if (!sp.catched) {
            std::cout << sp << std::endl;
        }

        assert(sp.catched);

    }

    // einrasten von unterschiedlichen ind auf einen punkt auf der grenze des polygons

    {
        struct Cmp {
            bool operator() (const StripPt &l, const StripPt &r) const {
                return l.ind < r.ind;
            }
        };

        std::map<Point3d, std::set<std::reference_wrapper<const StripPt>, Cmp>> equalPts;

        PolyStripsType::const_iterator itr;
        StripPtsType::const_iterator itr2;

        for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
            const PStrips &pStrips = itr->second;

            for (itr2 = pStrips.pts.begin(); itr2 != pStrips.pts.end(); ++itr2) {
                const StripPt &sp = itr2->second;

                if (sp.capt != Capt::NOT) {
                    equalPts[{sp.cutPt[0], sp.cutPt[1], sp.cutPt[2]}].emplace(sp);
                }
            }
        }

        vtkSmartPointer<vtkIdList> lines = vtkSmartPointer<vtkIdList>::New();

        vtkIdType i, numLines;

        decltype(equalPts)::const_iterator itr3;

        for (itr3 = equalPts.begin(); itr3 != equalPts.end(); ++itr3) {
            auto &pts = itr3->second;

            if (pts.size() > 1) {
                assert(pts.size() == 2);

                const StripPt &spA = *(pts.begin()),
                    &spB = *(std::next(pts.begin()));

                // spA auf spB kopieren und die line dazwischen löschen

#ifdef DEBUG
                std::cout << "collapsing " << spB.ind << " -> " << spA.ind << std::endl;
#endif

                contLines->GetPointCells(spB.ind, lines);

                numLines = lines->GetNumberOfIds();

                for (i = 0; i < numLines; i++) {
                    contLines->ReplaceCellPoint(lines->GetId(i), spB.ind, spA.ind);

                    // aktualisiert die links
                    contLines->RemoveReferenceToCell(spB.ind, lines->GetId(i));
                    contLines->ResizeCellList(spA.ind, 1);
                    contLines->AddReferenceToCell(spA.ind, lines->GetId(i));
                }

            }
        }

    }

    vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New(),
        line = vtkSmartPointer<vtkIdList>::New();

    vtkIdType i, numCells;

    StripPtsType::const_iterator itr5;

    for (itr = polyLines.begin(); itr != polyLines.end(); ++itr) {
        PStrips &pStrips = polyStrips[itr->first];

        const IdsType &lines = itr->second;
        const StripPtsType &pts = pStrips.pts;

        for (itr5 = pts.begin(); itr5 != pts.end(); ++itr5) {
            // ist der punkt auf einer kante, dürfen von ihm mehr als 2 linien ausgehen

            const StripPt &pt = itr5->second;

            if (pt.capt == Capt::NOT) {
                contLines->GetPointCells(pt.ind, cells);

                numCells = cells->GetNumberOfIds();

                std::set<vtkIdType> ends;

                for (i = 0; i < numCells; i++) {
                    contLines->GetCellPoints(cells->GetId(i), line);

                    ends.insert(pt.ind == line->GetId(0) ? line->GetId(1) : line->GetId(0));
                }

                if (ends.size() > 2) {
                    return true;
                }
            }
        }

        StripType strip;

        // zusammensetzen

        std::deque<vtkIdType> _lines(lines.begin(), lines.end());

        std::size_t i = 0;

        while (_lines.size() > 0) {
            vtkIdList *linePts = vtkIdList::New();
            contLines->GetCellPoints(_lines[i], linePts);

            vtkIdType indA = linePts->GetId(0),
                indB = linePts->GetId(1);

            if (strip.empty()) {
                strip.emplace_back(indA);
                strip.emplace_back(indB);

                _lines.erase(_lines.begin());
            } else {
                const StripPt &start = pts.at(strip.front().ind),
                    &end = pts.at(strip.back().ind);

                if (end.capt == Capt::NOT && end.ind == indA) {
                    strip.emplace_back(indB);
                    _lines.erase(_lines.begin()+i);
                    i = 0;

                } else if (end.capt == Capt::NOT && end.ind == indB) {
                    strip.emplace_back(indA);
                    _lines.erase(_lines.begin()+i);
                    i = 0;

                } else if (start.capt == Capt::NOT && start.ind == indA) {
                    strip.emplace_front(indB);
                    _lines.erase(_lines.begin()+i);
                    i = 0;

                } else if (start.capt == Capt::NOT && start.ind == indB) {
                    strip.emplace_front(indA);
                    _lines.erase(_lines.begin()+i);
                    i = 0;

                } else {
                    i++;

                }
            }

            const StripPt &_start = pts.at(strip.front().ind),
                &_end = pts.at(strip.back().ind);

            if ((_start.capt != Capt::NOT && _end.capt != Capt::NOT)
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

    // PolyStripsType::const_iterator itr2;
    // StripsType::const_iterator itr3;
    // StripType::const_iterator itr4;

    // for (itr2 = polyStrips.begin(); itr2 != polyStrips.end(); ++itr2) {
    //     const PStrips &pStrips = itr2->second;

    //     const StripsType &strips = pStrips.strips;
    //     const StripPtsType &pts = pStrips.pts;

    //     vtkIdList *cell = vtkIdList::New();
    //     pd->GetCellPoints(itr2->first, cell);

    //     Base base(pd->GetPoints(), cell);

    //     cell->Delete();

    //     AABB tree;

    //     std::vector<std::shared_ptr<Line>> lines;

    //     for (itr3 = strips.begin(); itr3 != strips.end(); ++itr3) {
    //         const StripType &strip = *itr3;

    //         for (itr4 = strip.begin(); itr4 != strip.end()-1; ++itr4) {
    //             const StripPt &spA = pts.at(itr4->ind),
    //                 &spB = pts.at((itr4+1)->ind);

    //             double ptA[2], ptB[2];

    //             Transform(spA.pt, ptA, base);
    //             Transform(spB.pt, ptB, base);

    //             vtkIdType grp = itr3-strips.begin(); // falscher typ!

    //             lines.emplace_back(new Line({ptA, itr4->ind}, {ptB, (itr4+1)->ind}, grp));
    //         }

    //     }

    //     for (auto &line : lines) {
    //         // std::cout << *line << std::endl;

    //         tree.InsertObj(line);
    //     }

    //     Bnds bnds(-E, E, -E, E);

    //     std::vector<std::shared_ptr<Line>>::const_iterator itr5;
    //     std::vector<std::shared_ptr<Obj>>::const_iterator itr6;

    //     for (itr5 = lines.begin(); itr5 != lines.end(); ++itr5) {
    //         auto found = tree.Search(*itr5);

    //         const Line &lA = **itr5;

    //         for (itr6 = found.begin(); itr6 != found.end(); ++itr6) {
    //             // std::cout << itr6->use_count() << std::endl;

    //             const Line &lB = dynamic_cast<Line&>(**itr6);

    //             // die linien dürfen nicht zum gleichen strip gehören und sich nicht an den enden berühren

    //             if (lA.grp != lB.grp
    //                 && lA.pA.id != lB.pA.id
    //                 && lA.pA.id != lB.pB.id
    //                 && lA.pB.id != lB.pA.id
    //                 && lA.pB.id != lB.pB.id
    //                 && Intersect2(lA.pA.pt, lA.pB.pt, lB.pA.pt, lB.pB.pt, bnds)) {
    //                 return true;
    //             }
    //         }
    //     }

    // }

    return false;

}

void vtkPolyDataBooleanFilter::RemoveDuplicates (IdsType &lines) {

    IdsType unique;

    // die indexe der enden auf übereinstimmung prüfen

    vtkIdList *linePtsA = vtkIdList::New();
    vtkIdList *linePtsB = vtkIdList::New();

    std::size_t i, j, numLines = lines.size();

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
        const StripPt &start = pStrips.pts[itr->front().ind],
            &end = pStrips.pts[itr->back().ind];

        if (start.ind != end.ind) {
            if (start.capt == Capt::NOT) {
                StripType s(itr->rbegin(), itr->rend()-1);
                itr->insert(itr->begin(), s.begin(), s.end());

            } else if (end.capt == Capt::NOT) {
                StripType s(itr->rbegin()+1, itr->rend());
                itr->insert(itr->end(), s.begin(), s.end());

            }
        }
    }
}

bool vtkPolyDataBooleanFilter::HasArea (const StripType &strip) const {
    bool area = true;

    std::size_t i, n = strip.size();

    if (n%2 == 1) {
        for (i = 0; i < (n-1)/2; i++) {
            area = strip[i].ind != strip[n-i-1].ind;
        }
    }

    return area;
}

void ComputeNormal (const StripPtsType &pts, const RefsType &poly, double *n) {
    n[0] = 0; n[1] = 0; n[2] = 0;

    RefsType::const_iterator itrA, itrB;

    for (itrA = poly.begin(); itrA != poly.end(); ++itrA) {
        itrB = itrA+1;

        if (itrB == poly.end()) {
            itrB = poly.begin();
        }

        const StripPtR &spA = *itrA,
            &spB = *itrB;

        auto pA = pts.find(spA.ind);
        auto pB = pts.find(spB.ind);

        const double *ptA = pA->second.cutPt,
            *ptB = pB->second.cutPt;

        n[0] += (ptA[1]-ptB[1])*(ptA[2]+ptB[2]);
        n[1] += (ptA[2]-ptB[2])*(ptA[0]+ptB[0]);
        n[2] += (ptA[0]-ptB[0])*(ptA[1]+ptB[1]);
    }

    vtkMath::Normalize(n);
}

void vtkPolyDataBooleanFilter::CutCells (vtkPolyData *pd, PolyStripsType &polyStrips) {
#ifdef DEBUG
    std::cout << "CutCells()" << std::endl;
#endif

    vtkPoints *pdPts = pd->GetPoints();

    vtkIdTypeArray *origCellIds = vtkIdTypeArray::SafeDownCast(pd->GetCellData()->GetScalars("OrigCellIds"));

    PolyStripsType::iterator itr;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {

        vtkIdType polyInd = itr->first;
        PStrips &pStrips = itr->second;

        /*if (polyInd != 1641) {
            continue;
        }*/

        StripsType &strips = pStrips.strips;
        StripPtsType &pts = pStrips.pts;

        IdsType &poly = pStrips.poly;

        vtkIdType origId = origCellIds->GetValue(polyInd);

#ifdef DEBUG
        std::cout << "polyInd " << polyInd << ", poly [";
        for (auto &p : poly) {
            std::cout << p << ", ";
        }
        std::cout << "]" << std::endl;
#endif

        std::size_t numPts = poly.size();

        std::map<vtkIdType, RefsType> edges;

        StripsType::iterator itr2;
        StripType::iterator itr3;

        // holes sammeln

        // HolesType holes;

        // for (itr2 = strips.begin(); itr2 != strips.end(); ++itr2) {
        //     StripType &s = *itr2;
        //     if (pts[s.front().ind].capt == Capt::NOT && pts[s.back().ind].capt == Capt::NOT) {
        //         IdsType hole;

        //         for (auto& sp : s) {
        //             hole.push_back(pdPts->InsertNextPoint(pts[sp.ind].pt));
        //         }

        //         // anfang und ende sind ja gleich
        //         hole.pop_back();

        //         holes.push_back(std::move(hole));

        //     }
        // }

        // holes löschen
        strips.erase(std::remove_if(strips.begin(), strips.end(), [&](const StripType &s) {
            return pts[s.front().ind].capt == Capt::NOT && pts[s.back().ind].capt == Capt::NOT; }), strips.end());

        for (itr2 = strips.begin(); itr2 != strips.end(); ++itr2) {
            StripType &strip = *itr2;

#ifdef DEBUG
            std::cout << (itr2-strips.begin()) << ". strip [";
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

            strip.front().side = Side::START;
            strip.back().side = Side::END;

            strip.front().ref = start.edge[0];
            strip.back().ref = end.edge[0];

            vtkIdType ind = itr2-strips.begin(); // falscher typ!

            strip.front().strip = strip.back().strip = ind;

            // nachfolgend könnte man dann anfang und ende weglassen

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

        std::map<vtkIdType, RefsType>::iterator itr4;

        IdsType::iterator itr5;
        StripType::reverse_iterator itr6;

        RefsType::iterator itr7;

        for (itr4 = edges.begin(); itr4 != edges.end(); ++itr4) {
            RefsType &edge = itr4->second;

#ifdef DEBUG
            std::cout << "edge (" << itr4->first << ", _)" << std::endl;
#endif

            std::sort(edge.begin(), edge.end(), [&](const StripPtR &a, const StripPtR &b) {
                const StripPt &a_ = pts[a.ind],
                    &b_ = pts[b.ind];

#ifdef DEBUG
                // std::cout << "a_: " << a_ << " -> strip " << a.strip << std::endl;
                // std::cout << "b_: " << b_ << " -> strip " << b.strip << std::endl;
#endif

                if (a_.ind == b_.ind) {
                    // strips beginnen im gleichen punkt

                    if (a.strip != b.strip) {
                        // gehören nicht dem gleichen strip an

                        StripType &stripA = strips[a.strip],
                            &stripB = strips[b.strip];

                        // andere enden ermitteln

                        const StripPtR &eA = (&a == &(stripA.front())) ? stripA.back() : stripA.front(),
                            &eB = (&b == &(stripB.front())) ? stripB.back() : stripB.front();

                        const StripPt &eA_ = pts[eA.ind],
                            &eB_ = pts[eB.ind];

#ifdef DEBUG
                        // std::cout << "eA_: " << eA_ << std::endl;
                        // std::cout << "eB_: " << eA_ << std::endl;
#endif

                        if (eA_.ind != eB_.ind) {
                            auto i = std::find(poly.begin(), poly.end(), itr4->first)-poly.begin(),
                                iA = std::find(poly.begin(), poly.end(), eA_.edge[0])-poly.begin(),
                                iB = std::find(poly.begin(), poly.end(), eB_.edge[0])-poly.begin();

                            double dA = static_cast<double>(Mod(iA-i, numPts))+eA_.t,
                                dB = static_cast<double>(Mod(iB-i, numPts))+eB_.t;

                            if (i == iA && a_.t > eA_.t) {
                               dA += static_cast<double>(numPts);
                            }

                            if (i == iB && b_.t > eB_.t) {
                               dB += static_cast<double>(numPts);
                            }

#ifdef DEBUG
                            // std::cout << "dA=" << dA << ", dB=" << dB << std::endl;
#endif

                            return dB < dA;
                        } else {
                            RefsType poly_;

                            if (a.side == Side::START) {
                                poly_.insert(poly_.end(), stripA.begin(), stripA.end());
                            } else {
                                poly_.insert(poly_.end(), stripA.rbegin(), stripA.rend());
                            }

                            if (b.side == Side::START) {
                                poly_.insert(poly_.end(), stripB.rbegin()+1, stripB.rend()-1);
                            } else {
                                poly_.insert(poly_.end(), stripB.begin()+1, stripB.end()-1);
                            }

                            double n[3];
                            ComputeNormal(pts, poly_, n);

                            double ang = vtkMath::Dot(pStrips.n, n);

#ifdef DEBUG
                            // std::cout << "ang=" << ang*180/PI << std::endl;
#endif

                            return ang < .999999;

                        }
                    } else {
                        // gleicher strip

                        StripType &strip = strips[a.strip];

                        if (HasArea(strip)) {
                            RefsType poly_(strip.begin(), strip.end()-1);

                            double n[3];
                            ComputeNormal(pts, poly_, n);

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
            for (auto& p : edge) {
                std::cout << p << std::endl;
            }
#endif

        }

        // baut die strips ein

        std::deque<IdsType> polys;
        polys.push_back(pStrips.poly);

        std::deque<IdsType>::iterator itr9;

        for (itr2 = strips.begin(); itr2 != strips.end(); ++itr2) {
            StripType &strip = *itr2;

            StripPtR &start = strip.front(),
                &end = strip.back();

#ifdef DEBUG
            std::cout << "strip " << start.strip
                << " , refs (" << start.ref << ", " << end.ref << ")"
                << std::endl;
#endif

            std::size_t cycle = 0;

            while (true) {

                if (cycle == polys.size()) {
                    break;
                }

                IdsType next = polys.front();
                polys.pop_front();

                std::vector<IdsType> newPolys(2);

                if (std::find(next.begin(), next.end(), start.ref) != next.end()) {
                    if (start.ref == end.ref) {
                        for (itr5 = next.begin(); itr5 != next.end(); ++itr5) {
                            newPolys[0].push_back(*itr5);

#ifdef DEBUG
                            std::cout << "adding " << *itr5 << " to 0" << std::endl;
#endif

                            if (*itr5 == start.ref) {
                                for (itr3 = strip.begin(); itr3 != strip.end(); ++itr3) {
                                    newPolys[0].push_back(itr3->desc[0]);

#ifdef DEBUG
                                    std::cout << "adding " << itr3->desc[0] << " to 0" << std::endl;
#endif

                                }
                            }
                        }

                        // strip selbst ist ein polygon

                        for (itr6 = strip.rbegin(); itr6 != strip.rend(); ++itr6) {
                            newPolys[1].push_back(itr6->desc[1]);

#ifdef DEBUG
                            std::cout << "adding " << itr6->desc[1] << " to 1" << std::endl;
#endif

                        }

                    } else {
                        std::size_t curr = 0;

                        for (itr5 = next.begin(); itr5 != next.end(); ++itr5) {
                            IdsType &poly = newPolys[curr];

                            poly.push_back(*itr5);

#ifdef DEBUG
                            std::cout << "adding " << *itr5 << " to " << curr << std::endl;
#endif

                            if (*itr5 == start.ref) {
                                for (itr3 = strip.begin(); itr3 != strip.end(); ++itr3) {
                                    poly.push_back(itr3->desc[0]);

#ifdef DEBUG
                                    std::cout << "adding " << itr3->desc[0] << " to " << curr << std::endl;
#endif

                                }

                                curr = curr == 0 ? 1 : 0;

                            } else if (*itr5 == end.ref) {
                                for (itr6 = strip.rbegin(); itr6 != strip.rend(); ++itr6) {
                                    poly.push_back(itr6->desc[1]);

#ifdef DEBUG
                                    std::cout << "adding " << itr6->desc[1] << " to " << curr << std::endl;
#endif

                                }

                                curr = curr == 0 ? 1 : 0;
                            }

                        }
                    }
                }

                if (newPolys[1].size() > 0) {

                    // refs aktualisieren

                    auto idx = std::distance(strips.begin(), itr2);

#ifdef DEBUG
                    std::cout << "idx " << idx << std::endl;
#endif

                    for (itr4 = edges.begin(); itr4 != edges.end(); ++itr4) {
                        RefsType &edge = itr4->second;

                        RefsType::iterator itrA;

                        for (itrA = edge.begin()+1; itrA != edge.end(); ++itrA) {
                            StripPtR &sp = *itrA;

                            if (sp.strip > idx) {
#ifdef DEBUG
                                std::cout << "sp: ind " << sp.ind << ", strip " << sp.strip << std::endl;
#endif

                                RefsType::const_reverse_iterator itrB(itrA);

                                vtkIdType _ind {-1},
                                    _strip {-1};

                                for (; itrB != edge.rend(); ++itrB) {
                                    const StripPtR &p = *itrB;

                                    if (p.strip != sp.strip) {
                                        if (p.strip <= idx) {
#ifdef DEBUG
                                            std::cout << "ref " << sp.ref;
#endif

                                            if (p.side == Side::END) {
                                                sp.ref = p.desc[0];
                                            } else {
                                                sp.ref = p.desc[1];
                                            }

#ifdef DEBUG
                                            std::cout << " -> " << sp.ref << " (from strip " << p.strip << ", ind " << p.ind << ")" << std::endl;
#endif

                                            _ind = p.ind;
                                            _strip = p.strip;

                                            break;

                                        }
                                    } else {
#ifdef DEBUG
                                        std::cout << "~1 ref " << sp.ref << " -> " << p.ref << " (from strip " << p.strip << ", ind " << p.ind << ")" << std::endl;
#endif

                                        sp.ref = p.ref;
                                        break;
                                    }
                                }

                                RefsType::const_iterator itrC(itrA);

                                ++itrC;

                                for (; itrC != edge.end(); ++itrC) {
                                    const StripPtR &p = *itrC;

                                    if (p.ind != sp.ind) {
                                        break;
                                    }

                                    if (p.strip <= idx) {
                                        if (p.ind == _ind && p.strip < _strip) {
                                            break;
                                        }

#ifdef DEBUG
                                        std::cout << "~2 ref " << sp.ref;
#endif

                                        if (p.side == Side::START) {
                                            sp.ref = p.desc[0];
                                        } else {
                                            sp.ref = p.desc[1];
                                        }

#ifdef DEBUG
                                        std::cout << " -> " << sp.ref << " (from strip " << p.strip << ", ind " << p.ind << ")" << std::endl;
#endif

                                        break;
                                    }

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
                                && b.strip == idx
                                && pts[a.ind].capt == Capt::A) { // sollte weg

#ifdef DEBUG
                                std::cout << "~3 ref " << a.ref;
#endif

                                if (b.side == Side::START) {
                                    a.ref = b.desc[0];
                                } else {
                                    a.ref = b.desc[1];
                                }

#ifdef DEBUG
                                std::cout << " -> " << a.ref << " (from strip " << b.strip << ", ind " << b.ind << ")" << std::endl;
#endif

                            }
                        }

                    }

                    // doppelte punkte entfernen

                    double pt[3];

                    for (auto &newPoly : newPolys) {
                        IdsType _newPoly;

                        std::map<vtkIdType, Point3d> _pts;

                        for (vtkIdType id : newPoly) {
                            pd->GetPoint(id, pt);
                            _pts.emplace(id, Point3d{pt[0], pt[1], pt[2]});
                        }

                        IdsType::const_iterator itrA, itrB;

                        for (itrA = newPoly.begin(); itrA != newPoly.end(); ++itrA) {
                            itrB = itrA+1;
                            if (itrB == newPoly.end()) {
                                itrB = newPoly.begin();
                            }

                            auto _a = _pts.find(*itrA);
                            auto _b = _pts.find(*itrB);

                            if (_a->second == _b->second) {
                                // doppelt
#ifdef DEBUG
                                std::cout << "removing " << *itrA << std::endl;
#endif
                            } else {
                                _newPoly.push_back(*itrA);
                            }
                        }

                        newPoly.swap(_newPoly);

                    }

                    // prüft, ob die erstellten polys gültig sind

                    if (newPolys[0].size() > 2) {
                        polys.push_back(newPolys[0]);
                    }

                    if (HasArea(strip) && newPolys[1].size() > 2) {
                        polys.push_back(newPolys[1]);
                    }

                    cycle = 0;

                    break;

                } else {
                    polys.push_back(next);

                    cycle++;
                }

            }

        }

        // erzeugte polys hinzufügen

        IdsType descIds;

        for (itr9 = polys.begin(); itr9 != polys.end(); ++itr9) {
            IdsType &p = *itr9;

            std::size_t num = p.size();

            vtkIdList *cell = vtkIdList::New();
            cell->SetNumberOfIds(num);

            vtkIdType i = 0;

            for (vtkIdType id : p) {
                cell->SetId(i++, id);
            }

            descIds.push_back(pd->InsertNextCell(VTK_POLYGON, cell));

            origCellIds->InsertNextValue(origId);

            cell->Delete();
        }

        pd->DeleteCell(polyInd);

        // holes verarbeiten

        // if (!holes.empty()) {
        //     _Wrapper w(pd, descIds, origId);

        //     for (auto& hole : holes) {
        //         w.Add(hole);
        //     }

        //     w.MergeAll();
        // }

    }

    pd->RemoveDeletedCells();
    pd->BuildCells();

}

void vtkPolyDataBooleanFilter::RestoreOrigPoints (vtkPolyData *pd, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "RestoreOrigPoints()" << std::endl;
#endif

    pd->BuildLinks();

    vtkKdTreePointLocator *loc = vtkKdTreePointLocator::New();
    loc->SetDataSet(pd);
    loc->BuildLocator();

    struct Cmp {
        bool operator() (const StripPt &l, const StripPt &r) const {
            return l.ind < r.ind;
        }
    };

    PolyStripsType::const_iterator itr;
    StripPtsType::const_iterator itr2;

    std::set<std::reference_wrapper<const StripPt>, Cmp> ends;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
        const PStrips &pStrips = itr->second;

        for (itr2 = pStrips.pts.begin(); itr2 != pStrips.pts.end(); ++itr2) {
            const StripPt &sp = itr2->second;

            if (sp.capt != Capt::NOT) {
                ends.emplace(sp);
            }

        }
    }

    vtkIdList *pts = vtkIdList::New();

    vtkIdType i, numPts;

    for (const StripPt &sp : ends) {
        FindPoints(loc, sp.cutPt, pts);
        numPts = pts->GetNumberOfIds();

        for (i = 0; i < numPts; i++) {
            pd->GetPoints()->SetPoint(pts->GetId(i), sp.pt);
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

    vtkKdTreePointLocator *loc = vtkKdTreePointLocator::New();
    loc->SetDataSet(pd);
    loc->BuildLocator();

    struct Cmp {
        bool operator() (const StripPt &l, const StripPt &r) const {
            return l.ind < r.ind;
        }
    };

    PolyStripsType::const_iterator itr;
    StripPtsType::const_iterator itr2;

    std::set<std::reference_wrapper<const StripPt>, Cmp> ends;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
        const PStrips &pStrips = itr->second;

        for (itr2 = pStrips.pts.begin(); itr2 != pStrips.pts.end(); ++itr2) {
            const StripPt &sp = itr2->second;

            if (sp.capt == Capt::A) {
                ends.emplace(sp);
            }
        }
    }

    vtkIdList *pts = vtkIdList::New();
    vtkIdList *cells = vtkIdList::New();

    vtkIdType i, j, numPts, numCells;

    for (const StripPt &sp : ends) {
        FindPoints(loc, sp.pt, pts);
        numPts = pts->GetNumberOfIds();

        for (i = 0; i < numPts; i++) {
            pd->GetPointCells(pts->GetId(i), cells);
            numCells = cells->GetNumberOfIds();

            if (numCells > 1) {
                for (j = 0; j < numCells; j++) {
                    pd->ReplaceCellPoint(cells->GetId(j), pts->GetId(i), pd->GetPoints()->InsertNextPoint(sp.pt));
                }
            }
        }
    }

    cells->Delete();
    pts->Delete();

    loc->FreeSearchStructure();
    loc->Delete();

}

void vtkPolyDataBooleanFilter::ResolveOverlaps (vtkPolyData *pd, vtkIdTypeArray *conts, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "ResolveOverlaps()" << std::endl;
#endif

    pd->BuildCells();
    pd->BuildLinks();

    contLines->BuildLinks();

    vtkKdTreePointLocator *loc = vtkKdTreePointLocator::New();
    loc->SetDataSet(pd);
    loc->BuildLocator();

    struct Cmp {
        bool operator() (const StripPt &l, const StripPt &r) const {
            return l.ind < r.ind;
        }
    };

    PolyStripsType::const_iterator itr;
    StripPtsType::const_iterator itr2;

    std::set<std::reference_wrapper<const StripPt>, Cmp> ends;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
        const PStrips &pStrips = itr->second;

        for (itr2 = pStrips.pts.begin(); itr2 != pStrips.pts.end(); ++itr2) {
            const StripPt &sp = itr2->second;

            if (sp.capt == Capt::EDGE) {
                ends.emplace(sp);
            }
        }
    }

    vtkIdList *links = vtkIdList::New(),
        *ptsA = vtkIdList::New(),
        *ptsB = vtkIdList::New(),
        *cells = vtkIdList::New(),
        *poly = vtkIdList::New();

    double ptA[3], ptB[3];

    vtkIdType i, j, k, numPtsA, numPtsB, numCells, numPts, idA, idB;

    std::map<Pair, IdsType> sharedPts;

    for (const StripPt &sp : ends) {
        contLines->GetPointCells(sp.ind, links);

        if (links->GetNumberOfIds() == 2
            && conts->GetValue(links->GetId(0)) != conts->GetValue(links->GetId(1))) {

            std::vector<Pair> history(sp.history.rbegin(), sp.history.rend());

            for (const Pair &p : history) {
                pd->GetPoint(p.f, ptA);
                pd->GetPoint(p.g, ptB);

                FindPoints(loc, ptA, ptsA);
                FindPoints(loc, ptB, ptsB);

                numPtsA = ptsA->GetNumberOfIds();
                numPtsB = ptsB->GetNumberOfIds();

                std::vector<Pair> cellsA, cellsB;

                for (i = 0; i < numPtsA; i++) {
                    pd->GetPointCells(ptsA->GetId(i), cells);
                    numCells = cells->GetNumberOfIds();

                    for (j = 0; j < numCells; j++) {
                        cellsA.emplace_back(ptsA->GetId(i), cells->GetId(j));
                    }
                }

                for (i = 0; i < numPtsB; i++) {
                    pd->GetPointCells(ptsB->GetId(i), cells);
                    numCells = cells->GetNumberOfIds();

                    for (j = 0; j < numCells; j++) {
                        cellsB.emplace_back(ptsB->GetId(i), cells->GetId(j));
                    }
                }

                for (auto &cellA : cellsA) {
                    for (auto &cellB : cellsB) {
                        if (cellA.g == cellB.g) {
                            // kante/poly existiert noch

                            pd->GetCellPoints(cellA.g, poly);
                            numPts = poly->GetNumberOfIds();

                            for (k = 0; k < numPts; k++) {
                                idA = poly->GetId(k);
                                idB = k+1 == numPts ? poly->GetId(0) : poly->GetId(k+1);

                                if (cellB.f == idA
                                    && cellA.f == idB) {

                                    IdsType &pts = sharedPts[Pair{sp.ind, cellA.g}];

                                    pts.push_back(cellA.f);
                                    pts.push_back(cellB.f);

                                    break;
                                }
                            }

                        }
                    }
                }


            }

        }
    }

    double pt[3];

    decltype(sharedPts)::iterator itr3;

    for (itr3 = sharedPts.begin(); itr3 != sharedPts.end(); ++itr3) {
        const Pair &p = itr3->first;

        IdsType &pts = itr3->second;

        std::sort(pts.begin(), pts.end());

        auto found = std::adjacent_find(pts.begin(), pts.end());

        if (found != pts.end()) {
            contLines->GetPoint(p.f, pt);

            i = pd->GetPoints()->InsertNextPoint(pt);
            pd->ReplaceCellPoint(p.g, *found, i);
        }
    }

    links->Delete();
    ptsA->Delete();
    ptsB->Delete();
    cells->Delete();
    poly->Delete();

    loc->FreeSearchStructure();
    loc->Delete();

}

void vtkPolyDataBooleanFilter::AddAdjacentPoints (vtkPolyData *pd, vtkIdTypeArray *conts, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "AddAdjacentPoints()" << std::endl;
#endif

    pd->BuildLinks();

    vtkIdTypeArray *origCellIds = vtkIdTypeArray::SafeDownCast(pd->GetCellData()->GetScalars("OrigCellIds"));

    struct Cmp {
        bool operator() (const StripPt &l, const StripPt &r) const {
            return l.t < r.t;
        }
    };

    vtkKdTreePointLocator *loc = vtkKdTreePointLocator::New();
    loc->SetDataSet(pd);
    loc->BuildLocator();

    vtkIdList *cells = vtkIdList::New(),
        *ptsA = vtkIdList::New(),
        *ptsB = vtkIdList::New(),
        *poly = vtkIdList::New(),
        *newPoly = vtkIdList::New();

    vtkIdType i, j, numCells, numPtsA, numPtsB, numPts, idA, idB;

    PolyStripsType::const_iterator itr;
    StripPtsType::const_iterator itr2;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
        const PStrips &pStrips = itr->second;

        std::map<Pair, std::set<std::reference_wrapper<const StripPt>, Cmp>> edgePts;

        for (itr2 = pStrips.pts.begin(); itr2 != pStrips.pts.end(); ++itr2) {
            const StripPt &sp = itr2->second;

            if (sp.capt == Capt::EDGE) {
                edgePts[{sp.edge[0], sp.edge[1]}].emplace(sp);
            }
        }

        decltype(edgePts)::iterator itr3;

        for (itr3 = edgePts.begin(); itr3 != edgePts.end(); ++itr3) {
            const Pair &edge = itr3->first;
            auto pts = itr3->second;

            StripPt spA, spB;

            pd->GetPoint(edge.f, spA.pt);
            pd->GetPoint(edge.g, spB.pt);

            spA.t = 0;
            spB.t = 1;

            pts.emplace(spA);
            pts.emplace(spB);

            std::vector<decltype(pts)::value_type> _pts(pts.rbegin(), pts.rend());

            decltype(_pts)::const_iterator itrA, itrB, itrC;

            itrA = _pts.begin();

            while (itrA != _pts.end()-1) {
                itrB = itrA+1;

                while (itrB != _pts.end()-1) {
                    contLines->GetPointCells(itrB->get().ind, cells);
                    numCells = cells->GetNumberOfIds();

                    // std::cout << "[";
                    // for (i = 0; i < numCells; i++) {
                    //     std::cout << conts->GetValue(cells->GetId(i)) << ", ";
                    // }
                    // std::cout << "]" << std::endl;

                    // break;

                    ++itrB;
                }

                if (std::distance(itrA, itrB) > 1) {
                    FindPoints(loc, itrA->get().pt, ptsA);
                    FindPoints(loc, itrB->get().pt, ptsB);

                    numPtsA = ptsA->GetNumberOfIds();
                    numPtsB = ptsB->GetNumberOfIds();

                    std::vector<Pair> polysA, polysB;

                    for (i = 0; i < numPtsA; i++) {
                        pd->GetPointCells(ptsA->GetId(i), cells);
                        numCells = cells->GetNumberOfIds();

                        for (j = 0; j < numCells; j++) {
                            polysA.emplace_back(cells->GetId(j), ptsA->GetId(i));
                        }
                    }

                    for (i = 0; i < numPtsB; i++) {
                        pd->GetPointCells(ptsB->GetId(i), cells);
                        numCells = cells->GetNumberOfIds();

                        for (j = 0; j < numCells; j++) {
                            polysB.emplace_back(cells->GetId(j), ptsB->GetId(i));
                        }
                    }

                    /*for (const Pair &pA : polysA) {
                        std::cout << "pA " << pA << std::endl;
                    }

                    for (const Pair &pB : polysB) {
                        std::cout << "pB " << pB << std::endl;
                    }*/

                    for (const Pair &pA : polysA) {
                        for (const Pair &pB : polysB) {
                            if (pA.f == pB.f && pd->GetCellType(pA.f) != VTK_EMPTY_CELL) {

                                pd->GetCellPoints(pA.f, poly);
                                numPts = poly->GetNumberOfIds();

                                newPoly->Reset();

                                for (j = 0; j < numPts; j++) {
                                    newPoly->InsertNextId(poly->GetId(j));

                                    idA = poly->GetId(j);
                                    idB = j+1 == numPts ? poly->GetId(0) : poly->GetId(j+1);

                                    if (pA.g == idA
                                        && pB.g == idB) {

                                        for (itrC = itrA+1; itrC != itrB; ++itrC) {
                                            // newPoly->InsertNextId(pd->InsertNextLinkedPoint(itrC->get().pt, 1));

                                            pd->InsertNextLinkedPoint(1);

                                            newPoly->InsertNextId(pd->GetPoints()->InsertNextPoint(itrC->get().pt));
                                        }

                                    }

                                    pd->RemoveReferenceToCell(idA, pA.f);
                                }

                                pd->DeleteCell(pA.f);

                                pd->InsertNextLinkedCell(VTK_POLYGON, newPoly->GetNumberOfIds(), newPoly->GetPointer(0));

                                origCellIds->InsertNextValue(origCellIds->GetValue(pA.f));

                                break;
                            }
                        }
                    }
                }

                itrA = itrB;
            }
        }
    }

    cells->Delete();
    ptsA->Delete();
    ptsB->Delete();
    poly->Delete();
    newPoly->Delete();

    loc->FreeSearchStructure();
    loc->Delete();

    pd->RemoveDeletedCells();

}

void vtkPolyDataBooleanFilter::MergePoints (vtkPolyData *pd, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "MergePoints()" << std::endl;
#endif

    pd->BuildCells();
    pd->BuildLinks();

    contLines->BuildLinks();

    vtkKdTreePointLocator *loc = vtkKdTreePointLocator::New();
    loc->SetDataSet(pd);
    loc->BuildLocator();

    PolyStripsType::const_iterator itr;
    StripsType::const_iterator itr2;

    std::map<vtkIdType, std::set<vtkIdType>> neighPts;

    vtkIdList *pts = vtkIdList::New();

    vtkIdType i, j, numPts;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
        const PStrips &pStrips = itr->second;

        for (itr2 = pStrips.strips.begin(); itr2 != pStrips.strips.end(); ++itr2) {
            const StripType &strip = *itr2;

            const StripPtR &spA = strip.front(),
                &spB = strip.back();

            const auto beforeA = pStrips.pts.find((strip.begin()+1)->ind),
                beforeB = pStrips.pts.find((strip.end()-2)->ind);

            FindPoints(loc, beforeA->second.pt, pts);
            numPts = pts->GetNumberOfIds();

            for (i = 0; i < numPts; i++) {
                neighPts[spA.ind].insert(pts->GetId(i));
            }

            FindPoints(loc, beforeB->second.pt, pts);
            numPts = pts->GetNumberOfIds();

            for (i = 0; i < numPts; i++) {
                neighPts[spB.ind].insert(pts->GetId(i));
            }
        }
    }

    double pt[3];

    vtkIdList *polys = vtkIdList::New(),
        *poly = vtkIdList::New();

    vtkIdType ind, polyId, _numPts, before, after;

    decltype(neighPts)::const_iterator itr3;

    for (itr3 = neighPts.begin(); itr3 != neighPts.end(); ++itr3) {
        const auto &inds = itr3->second;

        std::map<Point3d, std::vector<MergePt>> pairs;

        contLines->GetPoint(itr3->first, pt);

        FindPoints(loc, pt, pts);
        numPts = pts->GetNumberOfIds();

        for (i = 0; i < numPts; i++) {
            ind = pts->GetId(i);

            pd->GetPointCells(ind, polys);

            if (polys->GetNumberOfIds() > 0) {
                polyId = polys->GetId(0);

                pd->GetCellPoints(polyId, poly);
                _numPts = poly->GetNumberOfIds();

                for (j = 0; j < _numPts; j++) {
                    if (poly->GetId(j) == ind) {
                        break;
                    }
                }

                // wieder davor und danach ermitteln

                before = poly->GetId(j == 0 ? _numPts-1 : j-1);
                after = poly->GetId(j+1 == _numPts ? 0 : j+1);

                if (std::find(inds.begin(), inds.end(), before) == inds.end()) {
                    pd->GetPoint(before, pt);
                    pairs[{pt[0], pt[1], pt[2]}].emplace_back(polyId, ind);
                }

                if (std::find(inds.begin(), inds.end(), after) == inds.end()) {
                    pd->GetPoint(after, pt);
                    pairs[{pt[0], pt[1], pt[2]}].emplace_back(polyId, ind);
                }

            }
        }

        std::deque<std::deque<std::reference_wrapper<const MergePt>>> Pairs;

        decltype(pairs)::const_iterator itr4;

        for (itr4 = pairs.begin(); itr4 != pairs.end(); ++itr4) {
            const auto &p = itr4->second;

            if (p.size() == 2) {
                auto _pts = {std::ref(p.front()), std::ref(p.back())}; // std::initializer_list
                Pairs.emplace_back(_pts);
            }
        }

        decltype(Pairs)::iterator itr5;

        /*for (itr5 = Pairs.begin(); itr5 != Pairs.end(); ++itr5) {
            for (auto &p : *itr5) {
                std::cout << p.get() << ", ";
            }
            std::cout << std::endl;
        }*/

        decltype(Pairs)::value_type group;

        decltype(group)::const_iterator itr6;

        while (!Pairs.empty()) {
            if (group.empty()) {
                group = Pairs.front();
                Pairs.pop_front();
            }

            itr5 = Pairs.begin();

            while (itr5 != Pairs.end()) {
                const auto &next = *itr5;

                if (next.front().get() == group.front().get()) {
                    group.emplace_front(next.back());
                    Pairs.erase(itr5);
                    itr5 = Pairs.begin();
                } else if (next.front().get() == group.back().get()) {
                    group.emplace_back(next.back());
                    Pairs.erase(itr5);
                    itr5 = Pairs.begin();
                } else if (next.back().get() == group.front().get()) {
                    group.emplace_front(next.front());
                    Pairs.erase(itr5);
                    itr5 = Pairs.begin();
                } else if (next.back().get() == group.back().get()) {
                    group.emplace_back(next.front());
                    Pairs.erase(itr5);
                    itr5 = Pairs.begin();
                } else {
                    ++itr5;
                }
            }

            if (itr5 == Pairs.end()) {
                for (itr6 = group.begin()+1; itr6 != group.end(); ++itr6) {
                    pd->ReplaceCellPoint(itr6->get().polyId, itr6->get().ind, group.front().get().ind);
                }

                group.clear();
            }
        }

    }

    polys->Delete();
    poly->Delete();
    pts->Delete();

    loc->FreeSearchStructure();
    loc->Delete();

}

enum class Congr {
    EQUAL,
    OPPOSITE,
    NOT
};

class PolyAtEdge {
    vtkPolyData *pd;

public:
    PolyAtEdge (vtkPolyData *_pd, vtkIdType _polyId, vtkIdType _ptIdA, vtkIdType _ptIdB) : pd(_pd), polyId(_polyId), ptIdA(_ptIdA), ptIdB(_ptIdB), loc(Loc::NONE) {

        double ptA[3], ptB[3];

        pd->GetPoint(ptIdA, ptA);
        pd->GetPoint(ptIdB, ptB);

        vtkMath::Subtract(ptB, ptA, e);
        vtkMath::Normalize(e);

        const vtkIdType *poly;

        vtkIdType numPts;
        pd->GetCellPoints(polyId, numPts, poly);

        ComputeNormal(pd->GetPoints(), n, numPts, poly);

        vtkMath::Cross(e, n, r);

    }

    vtkIdType polyId, ptIdA, ptIdB;
    double n[3], e[3], r[3];

    Loc loc;

    friend std::ostream& operator<< (std::ostream &out, const PolyAtEdge &p) {
        out << "polyId " << p.polyId << ", ptIdA " << p.ptIdA << ", ptIdB " << p.ptIdB;
        return out;
    }

    static constexpr double eps {.99999999}; // ~.0081deg

    Congr IsCongruent (const PolyAtEdge &p) const {
        double cong = vtkMath::Dot(n, p.n);

        if (cong > eps || cong < -eps) {
            double ang = vtkMath::Dot(r, p.r);

            if (ang > eps) {
                if (cong > eps) {
                    // normalen sind gleich ausgerichtet
                    return Congr::EQUAL;
                } else {
                    return Congr::OPPOSITE;
                }
            }
        }

        return Congr::NOT;
    }

};


class PolyPair {
public:
    PolyPair (PolyAtEdge _pA, PolyAtEdge _pB) : pA(_pA), pB(_pB) {}

    PolyAtEdge pA, pB;

    void GetLoc (PolyAtEdge &pT, vtkIdType mode) {

        Congr cA = pA.IsCongruent(pT),
            cB = pB.IsCongruent(pT);

#ifdef DEBUG
        std::cout << "GetLoc() -> polyId " << pT.polyId
                  << ", cA " << cA
                  << ", cB " << cB
                  << std::endl;

        if (cA != Congr::NOT || cB != Congr::NOT) {
            assert(cA != cB);
        }

#endif

        if (cA == Congr::EQUAL || cA == Congr::OPPOSITE) {
            if (cA == Congr::OPPOSITE) {
                // normalen sind entgegengesetzt gerichtet

                if (mode == OPER_INTERSECTION) {
                    pA.loc = Loc::OUTSIDE;
                    pT.loc = Loc::OUTSIDE;
                } else {
                    pA.loc = Loc::INSIDE;
                    pT.loc = Loc::INSIDE;
                }
            } else if (mode == OPER_UNION || mode == OPER_INTERSECTION) {
                pA.loc = Loc::INSIDE;
                pT.loc = Loc::OUTSIDE;
            }

        } else if (cB == Congr::EQUAL || cB == Congr::OPPOSITE) {
            if (cB == Congr::OPPOSITE) {
                // normalen sind entgegengesetzt gerichtet

                if (mode == OPER_INTERSECTION) {
                    pB.loc = Loc::OUTSIDE;
                    pT.loc = Loc::OUTSIDE;
                } else {
                    pB.loc = Loc::INSIDE;
                    pT.loc = Loc::INSIDE;
                }
            } else if (mode == OPER_UNION || mode == OPER_INTERSECTION) {
                pB.loc = Loc::INSIDE;
                pT.loc = Loc::OUTSIDE;
            }

        } else {
            double alpha = GetAngle(pA.r, pB.r, pA.e),
                beta = GetAngle(pA.r, pT.r, pA.e);

            if (beta > alpha) {
                pT.loc = Loc::INSIDE;
            } else {
                pT.loc = Loc::OUTSIDE;
            }
        }
    }

};


std::shared_ptr<PolyPair> GetEdgePolys (vtkPolyData *pd, vtkIdList *ptsA, vtkIdList *ptsB) {

#ifdef DEBUG
    std::cout << "GetEdgePolys()" << std::endl;
#endif

    std::vector<Pair> p;

    vtkIdType numPtsA = ptsA->GetNumberOfIds(),
        numPtsB = ptsB->GetNumberOfIds();

    vtkIdList *polys = vtkIdList::New();

    vtkIdType i, j, numCells;

    for (i = 0; i < numPtsA; i++) {
        pd->GetPointCells(ptsA->GetId(i), polys);
        numCells = polys->GetNumberOfIds();

        for (j = 0; j < numCells; j++) {
            p.emplace_back(ptsA->GetId(i), polys->GetId(j));
        }
    }

    for (i = 0; i < numPtsB; i++) {
        pd->GetPointCells(ptsB->GetId(i), polys);
        numCells = polys->GetNumberOfIds();

        for (j = 0; j < numCells; j++) {
            p.emplace_back(ptsB->GetId(i), polys->GetId(j));
        }
    }

    polys->Delete();

    std::map<vtkIdType, IdsType> pEdges;

    std::vector<Pair>::const_iterator itr;
    for (itr = p.begin(); itr != p.end(); ++itr) {
        pEdges[itr->g].push_back(itr->f);
    }

    std::vector<PolyAtEdge> opp;

    vtkIdList *poly = vtkIdList::New();

    vtkIdType numPts, a, b;

    std::map<vtkIdType, IdsType>::const_iterator itr2;

    for (itr2 = pEdges.begin(); itr2 != pEdges.end(); ++itr2) {
        const IdsType &pts = itr2->second;

        if (pts.size() > 1) {
            pd->GetCellPoints(itr2->first, poly);
            numPts = poly->GetNumberOfIds();

            for (i = 0; i < numPts; i++) {
                a = poly->GetId(i);
                b = i+1 == numPts ? poly->GetId(0) : poly->GetId(i+1);

                if (std::find(pts.begin(), pts.end(), a) != pts.end()
                    && std::find(pts.begin(), pts.end(), b) != pts.end()) {

                    opp.emplace_back(pd, itr2->first, a, b);
                }
            }

        }
    }

    poly->Delete();

#ifdef DEBUG
    for (auto &op : opp) {
        std::cout << op << std::endl;
    }
#endif

    if (opp.size() != 2) {
        return nullptr;
    }

    return std::make_shared<PolyPair>(opp[0], opp[1]);

}

bool vtkPolyDataBooleanFilter::CombineRegions () {

#ifdef DEBUG
    std::cout << "CombineRegions()" << std::endl;
#endif

    vtkSmartPointer<vtkPolyData> filterdA = vtkSmartPointer<vtkPolyData>::New();
    filterdA->DeepCopy(modPdA);

    vtkSmartPointer<vtkPolyData> filterdB = vtkSmartPointer<vtkPolyData>::New();
    filterdB->DeepCopy(modPdB);

    // ungenutzte punkte löschen
    vtkSmartPointer<vtkCleanPolyData> cleanA = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanA->PointMergingOff();
    cleanA->SetInputData(filterdA);

    vtkSmartPointer<vtkCleanPolyData> cleanB = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanB->PointMergingOff();
    cleanB->SetInputData(filterdB);

    // regionen mit skalaren ausstatten
    vtkSmartPointer<vtkPolyDataConnectivityFilter> cfA = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    cfA->SetExtractionModeToAllRegions();
    cfA->ColorRegionsOn();
    cfA->SetInputConnection(cleanA->GetOutputPort());

    vtkSmartPointer<vtkPolyDataConnectivityFilter> cfB = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    cfB->SetExtractionModeToAllRegions();
    cfB->ColorRegionsOn();
    cfB->SetInputConnection(cleanB->GetOutputPort());

    cfA->Update();
    cfB->Update();

    vtkPolyData *pdA = cfA->GetOutput();
    vtkPolyData *pdB = cfB->GetOutput();

#ifdef DEBUG
    std::cout << "Exporting modPdA_8.vtk" << std::endl;
    WriteVTK("modPdA_8.vtk", pdA);

    std::cout << "Exporting modPdB_8.vtk" << std::endl;
    WriteVTK("modPdB_8.vtk", pdB);
#endif

    resultA->ShallowCopy(contLines);

    if (OperMode == OPER_NONE) {
        resultB->ShallowCopy(pdA);
        resultC->ShallowCopy(pdB);

        return false;
    }

    vtkSmartPointer<vtkKdTreePointLocator> plA = vtkSmartPointer<vtkKdTreePointLocator>::New();
    plA->SetDataSet(pdA);
    plA->BuildLocator();

    vtkSmartPointer<vtkKdTreePointLocator> plB = vtkSmartPointer<vtkKdTreePointLocator>::New();
    plB->SetDataSet(pdB);
    plB->BuildLocator();

    pdA->BuildLinks();
    pdB->BuildLinks();

    vtkIdTypeArray *scalarsA = vtkIdTypeArray::SafeDownCast(pdA->GetPointData()->GetScalars());
    vtkIdTypeArray *scalarsB = vtkIdTypeArray::SafeDownCast(pdB->GetPointData()->GetScalars());

    vtkSmartPointer<vtkIdList> line = vtkSmartPointer<vtkIdList>::New();

    double ptA[3], ptB[3];

    vtkSmartPointer<vtkIdList> fptsA = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> lptsA = vtkSmartPointer<vtkIdList>::New();

    vtkSmartPointer<vtkIdList> fptsB = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> lptsB = vtkSmartPointer<vtkIdList>::New();

    std::map<vtkIdType, Loc> locsA, locsB;

    vtkIdType i, j, numLines = contLines->GetNumberOfCells();

    for (i = 0; i < numLines; i++) {

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

        vtkIdType notLocated = 0;

        for (j = 0; j < fptsA->GetNumberOfIds(); j++) {
            if (locsA.count(scalarsA->GetValue(fptsA->GetId(j))) == 0) {
                notLocated++;
            }
        }

        for (j = 0; j < fptsB->GetNumberOfIds(); j++) {
            if (locsB.count(scalarsB->GetValue(fptsB->GetId(j))) == 0) {
                notLocated++;
            }
        }

        if (notLocated == 0) {
            continue;
        }

#endif

        FindPoints(plA, ptB, lptsA);
        FindPoints(plB, ptB, lptsB);

        auto ppA = GetEdgePolys(pdA, fptsA, lptsA);
        auto ppB = GetEdgePolys(pdB, fptsB, lptsB);

        if (ppA && ppB) {

            ppB->GetLoc(ppA->pA, OperMode);
            ppB->GetLoc(ppA->pB, OperMode);

            ppA->GetLoc(ppB->pA, OperMode);
            ppA->GetLoc(ppB->pB, OperMode);

            vtkIdType fsA = scalarsA->GetValue(ppA->pA.ptIdA);
            vtkIdType lsA = scalarsA->GetValue(ppA->pB.ptIdA);

            vtkIdType fsB = scalarsB->GetValue(ppB->pA.ptIdA);
            vtkIdType lsB = scalarsB->GetValue(ppB->pB.ptIdA);

#ifdef DEBUG
            std::cout << "polyId " << ppA->pA.polyId << ", sA " << fsA << ", loc " << ppA->pA.loc << std::endl;
            std::cout << "polyId " << ppA->pB.polyId << ", sA " << lsA << ", loc " << ppA->pB.loc << std::endl;
            std::cout << "polyId " << ppB->pA.polyId << ", sB " << fsB << ", loc " << ppB->pA.loc << std::endl;
            std::cout << "polyId " << ppB->pB.polyId << ", sB " << lsB << ", loc " << ppB->pB.loc << std::endl;

            if (locsA.count(fsA) > 0 && locsA[fsA] != ppA->pA.loc) {
                std::cout << "sA " << fsA << ": " << locsA[fsA] << " -> " << ppA->pA.loc << std::endl;
            }

            if (locsA.count(lsA) > 0 && locsA[lsA] != ppA->pB.loc) {
                std::cout << "sA " << lsA << ": " << locsA[lsA] << " -> " << ppA->pB.loc << std::endl;
            }

            if (locsB.count(fsB) > 0 && locsB[fsB] != ppB->pA.loc) {
                std::cout << "sB " << fsB << ": " << locsB[fsB] << " -> " << ppB->pA.loc << std::endl;
            }

            if (locsB.count(lsB) > 0 && locsB[lsB] != ppB->pB.loc) {
                std::cout << "sB " << lsB << ": " << locsB[lsB] << " -> " << ppB->pB.loc << std::endl;
            }

#endif

            locsA.emplace(fsA, ppA->pA.loc);
            locsA.emplace(lsA, ppA->pB.loc);

            locsB.emplace(fsB, ppB->pA.loc);
            locsB.emplace(lsB, ppB->pB.loc);

        } else {
            return true;
        }

    }

    // reale kombination der ermittelten regionen

    Loc comb[] = {Loc::OUTSIDE, Loc::OUTSIDE};

    if (OperMode == OPER_INTERSECTION) {
        comb[0] = Loc::INSIDE;
        comb[1] = Loc::INSIDE;
    } else if (OperMode == OPER_DIFFERENCE) {
        comb[1] = Loc::INSIDE;
    } else if (OperMode == OPER_DIFFERENCE2) {
        comb[0] = Loc::INSIDE;
    }

    vtkIdType numA = cfA->GetNumberOfExtractedRegions(),
        numB = cfB->GetNumberOfExtractedRegions();

    cfA->SetExtractionModeToSpecifiedRegions();
    cfB->SetExtractionModeToSpecifiedRegions();

    std::map<vtkIdType, Loc>::const_iterator itr;

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

    // nicht beteiligte regionen hinzufügen

    if (OperMode == OPER_UNION || OperMode == OPER_DIFFERENCE) {
        for (i = 0; i < numA; i++) {
            if (locsA.count(i) == 0) {
                cfA->AddSpecifiedRegion(i);
            }
        }
    }

    if (OperMode == OPER_UNION || OperMode == OPER_DIFFERENCE2) {
        for (i = 0; i < numB; i++) {
            if (locsB.count(i) == 0) {
                cfB->AddSpecifiedRegion(i);
            }
        }
    }

    // nach innen zeigende normalen umkehren

    cfA->Update();
    cfB->Update();

    vtkPolyData *regsA = cfA->GetOutput();
    vtkPolyData *regsB = cfB->GetOutput();

    scalarsA = vtkIdTypeArray::SafeDownCast(regsA->GetPointData()->GetScalars());
    scalarsB = vtkIdTypeArray::SafeDownCast(regsB->GetPointData()->GetScalars());

    if (OperMode != OPER_INTERSECTION) {
        if (comb[0] == Loc::INSIDE) {
            for (i = 0; i < regsA->GetNumberOfCells(); i++) {
                if (locsA.count(scalarsA->GetValue(i)) == 1) {
                    regsA->ReverseCell(i);
                }
            }
        }

        if (comb[1] == Loc::INSIDE) {
            for (i = 0; i < regsB->GetNumberOfCells(); i++) {
                if (locsB.count(scalarsB->GetValue(i)) == 1) {
                    regsB->ReverseCell(i);
                }
            }
        }
    }

    // OrigCellIds und CellData

    vtkIdTypeArray *origCellIdsA = vtkIdTypeArray::SafeDownCast(regsA->GetCellData()->GetScalars("OrigCellIds"));
    vtkIdTypeArray *origCellIdsB = vtkIdTypeArray::SafeDownCast(regsB->GetCellData()->GetScalars("OrigCellIds"));

    vtkIdTypeArray *newOrigCellIdsA = vtkIdTypeArray::New();
    newOrigCellIdsA->SetName("OrigCellIdsA");

    vtkIdTypeArray *newOrigCellIdsB = vtkIdTypeArray::New();
    newOrigCellIdsB->SetName("OrigCellIdsB");

    vtkCellData *newCellDataA = vtkCellData::New();
    newCellDataA->CopyAllocate(cellDataA);

    vtkCellData *newCellDataB = vtkCellData::New();
    newCellDataB->CopyAllocate(cellDataB);

    vtkIdType cellId;

    for (i = 0; i < regsA->GetNumberOfCells(); i++) {
        cellId = cellIdsA->GetValue(origCellIdsA->GetValue(i));

        newOrigCellIdsA->InsertNextValue(cellId);
        newOrigCellIdsB->InsertNextValue(-1);

        newCellDataA->CopyData(cellDataA, cellId, i);
    }

    for (i = 0; i < regsB->GetNumberOfCells(); i++) {
        cellId = cellIdsB->GetValue(origCellIdsB->GetValue(i));

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

    // resultC bleibt hier leer

    resultB->ShallowCopy(cfPd);

    resultB->GetCellData()->AddArray(newOrigCellIdsA);
    resultB->GetCellData()->AddArray(newOrigCellIdsB);

    // aufräumen

    cfApp->Delete();
    cleanApp->Delete();
    app->Delete();

    newOrigCellIdsB->Delete();
    newOrigCellIdsA->Delete();

    plB->FreeSearchStructure();
    plA->FreeSearchStructure();

    return false;

}

// void _Wrapper::MergeAll () {
//     vtkPoints *pdPts = pd->GetPoints();

//     vtkIdTypeArray *origCellIds = vtkIdTypeArray::SafeDownCast(pd->GetCellData()->GetScalars("OrigCellIds"));

//     vtkIdList *cell = vtkIdList::New();

//     // descendants in holes einfügen

//     vtkIdType i, numPts;

//     std::set<vtkIdType> outerIds;

//     for (vtkIdType id : descIds) {
//         pd->GetCellPoints(id, cell);
//         numPts = cell->GetNumberOfIds();

//         IdsType hole;

//         for (i = 0; i < numPts; i++) {
//             hole.push_back(cell->GetId(i));
//         }

//         outerIds.insert(hole.begin(), hole.end());

//         holes.push_back(std::move(hole));

//         // löscht
//         pd->DeleteCell(id);
//     }

//     base = Base(pdPts, cell);

//     double a[3], b[2];

//     Merger m;

//     for (const IdsType &hole : holes) {
//         PolyType casted;

//         for (vtkIdType id : hole) {
//             pd->GetPoint(id, a);

//             Transform(a, b, base);

//             casted.emplace_back(b, id);
//         }

//         if (TestCW(casted)) {
//             std::reverse(casted.begin(), casted.end());
//         }

//         m.AddPoly(casted);
//     }

//     PolysType merged;
//     m.GetMerged(merged);

//     std::set<vtkIdType> usedIds;

//     for (auto &poly : merged) {
//         cell->Reset();

//         // poly ist immer ccw

//         std::map<vtkIdType, vtkIdType> repl;

//         for (auto &p : poly) {
//             if (outerIds.count(p.id) == 1 || usedIds.count(p.id) == 0) {
//                 repl.emplace(p.id, p.id);
//             } else {
//                 pd->GetPoint(p.id, a);

//                 repl.emplace(p.id, pdPts->InsertNextPoint(a));

//             }
//         }

//         for (auto &p : poly) {
//             cell->InsertNextId(repl[p.id]);
//             usedIds.insert(p.id);
//         }

//         pd->InsertNextCell(VTK_POLYGON, cell);
//         origCellIds->InsertNextValue(origId);

//     }

//     cell->Delete();

// }
