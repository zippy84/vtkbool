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

#include <map>
#include <deque>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <functional>
#include <queue>
#include <memory>
#include <tuple>

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
#include <vtkAppendPolyData.h>
#include <vtkKdTreePointLocator.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkSmartPointer.h>
#include <vtkModifiedBSPTree.h>
#include <vtkCellArrayIterator.h>
#include <vtkKdTree.h>
#include <vtkCellIterator.h>

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

int vtkPolyDataBooleanFilter::RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

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

            auto cleanA = vtkSmartPointer<vtkCleanPolyData>::New();
            cleanA->SetOutputPointsPrecision(DOUBLE_PRECISION);
            cleanA->SetTolerance(1e-6);
            cleanA->SetInputData(pdA);
            cleanA->Update();

            auto cleanB = vtkSmartPointer<vtkCleanPolyData>::New();
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

            auto cl = vtkSmartPointer<vtkPolyDataContactFilter>::New();
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
                return 1;
            }

            // in den CellDatas steht drin, welche polygone einander schneiden

            contsA = vtkIdTypeArray::SafeDownCast(contLines->GetCellData()->GetScalars("cA"));
            contsB = vtkIdTypeArray::SafeDownCast(contLines->GetCellData()->GetScalars("cB"));

            vtkIdTypeArray *sourcesA = vtkIdTypeArray::SafeDownCast(contLines->GetCellData()->GetScalars("sourcesA"));
            vtkIdTypeArray *sourcesB = vtkIdTypeArray::SafeDownCast(contLines->GetCellData()->GetScalars("sourcesB"));

            vtkIdType i, numPts = contLines->GetNumberOfPoints();

            auto cells = vtkSmartPointer<vtkIdList>::New();

            for (i = 0; i < numPts; i++) {
                contLines->GetPointCells(i, cells);

                if (cells->GetNumberOfIds() == 1) {
                    vtkErrorMacro("Contact ends suddenly.");
                    return 1;
                }

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

            // löscht bestimmte strips

            if (CleanStrips()) {
                vtkErrorMacro("There is no contact.");
                return 1;
            }

            times.push_back(clock::now()-start);

            // trennt die polygone an den linien

            start = clock::now();

            if (CutCells(modPdA, polyStripsA) ||
                CutCells(modPdB, polyStripsB)) {

                vtkErrorMacro("CutCells failed.");
                return 1;
            }

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

            ResolveOverlaps(modPdA, polyStripsA);
            ResolveOverlaps(modPdB, polyStripsB);

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

    double a[3], b[3], sA[3], sB[3], u[3], v[3], w[3], n, t, d;

    std::map<vtkIdType, vtkIdType> allPts;

    std::map<vtkIdType, vtkIdType> links;

    auto line = vtkSmartPointer<vtkIdList>::New();

    for (auto lineId : lines) {
        contLines->GetCellPoints(lineId, line);

        // std::cout << "? " << contsA->GetValue(lineId)
        //     << ", " << contsB->GetValue(lineId)
        //     << ", [" << line->GetId(0)
        //     << ", " << line->GetId(1)
        //     << "]"
        //     << std::endl;

        allPts.emplace(line->GetId(0), sources->GetTypedComponent(lineId, 0));
        allPts.emplace(line->GetId(1), sources->GetTypedComponent(lineId, 1));

        links[line->GetId(0)]++;
        links[line->GetId(1)]++;
    }

    decltype(allPts)::const_iterator itr;

    for (itr = allPts.begin(); itr != allPts.end(); ++itr) {
        StripPt sp;
        sp.ind = itr->first;

        // die koordinaten
        contLines->GetPoint(sp.ind, sp.pt);

        IdsType::const_iterator itrA, itrB;

        for (itrA = poly.begin(); itrA != poly.end(); ++itrA) {
            itrB = itrA+1;

            if (itrB == poly.end()) {
                itrB = poly.begin();
            }

            if (itr->second != NOTSET && *itrA != itr->second) {
                continue;
            }

            pd->GetPoint(*itrA, a);
            pd->GetPoint(*itrB, b);

            vtkMath::Subtract(a, sp.pt, sA);
            vtkMath::Subtract(b, sp.pt, sB);

            // richtungsvektor und länge der kante

            vtkMath::Subtract(b, a, u);
            n = vtkMath::Norm(u);

            // d und t zur kante

            vtkMath::Subtract(sp.pt, a, v);
            t = vtkMath::Dot(v, u)/(n*n);

            vtkMath::Cross(v, u, w);
            d = vtkMath::Norm(w)/n;

            if (d < 1e-5 && t > -1e-5 && t < 1+1e-5) {
                sp.edge[0] = *itrA;
                sp.edge[1] = *itrB;

                sp.t = std::min(1., std::max(0., t));

                if (vtkMath::Norm(sA) < 1e-5) {
                    Cpy(sp.captPt, a, 3);
                    sp.capt = Capt::A;

                } else if (vtkMath::Norm(sB) < 1e-5) {
                    Cpy(sp.captPt, b, 3);
                    sp.capt = Capt::B;

                } else {
                    // u ist nicht normiert
                    vtkMath::MultiplyScalar(u, t);

                    double x[3];
                    vtkMath::Add(a, u, x);

                    // projektion
                    Cpy(sp.captPt, x, 3);

                    sp.capt = Capt::EDGE;

                }
            }

            // std::cout << "? "
            //     << sp.ind
            //     << ", " << d
            //     << ", " << t
            //     << ", " << sp.capt
            //     << std::endl;
        }

        if (itr->second != NOTSET && sp.edge[0] == NOTSET) {
            sp.catched = false;
        }

        if (sp.capt == Capt::NOT && links[sp.ind] > 2) {
            sp.capt = Capt::BRANCHED;
        }

        pts.emplace(sp.ind, std::move(sp));

    }

    StripPtsType::iterator itr2;

    for (itr2 = pts.begin(); itr2 != pts.end(); ++itr2) {
        StripPt &sp = itr2->second;

        if (sp.capt & Capt::BOUNDARY) {
            if (sp.capt == Capt::B) {
                sp.t = 0;

                sp.edge[0] = sp.edge[1];

                auto itrA = std::find(poly.begin(), poly.end(), sp.edge[0]),
                    itrB = itrA+1;

                if (itrB == poly.end()) {
                    itrB = poly.begin();
                }

                sp.edge[1] = *itrB;

                sp.capt = Capt::A;

            }

            // für den schnitt werden die eingerasteten koordinaten verwendet

            Cpy(sp.cutPt, sp.captPt, 3);
        } else {

            Cpy(sp.cutPt, sp.pt, 3);
        }

    }

#ifdef DEBUG
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

        polyStrips.emplace(std::piecewise_construct,
            std::forward_as_tuple(itr->first),
            std::forward_as_tuple(pd, itr->first));

        PStrips &pStrips = polyStrips.at(itr->first);

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
            const PStrips &pStrips = polyStrips.at(itr->first);

            try {
                const StripPt &corr = pStrips.pts.at(sp.ind);

                if (&corr != &sp) {
                    if (corr.capt == Capt::A) {
                        sp.capt = Capt::A;
                        sp.edge[0] = corr.edge[0];
                        sp.edge[1] = Next(polyStrips.at(sp.polyId).poly, sp.edge[0]);

                        sp.t = 0;

                        Cpy(sp.captPt, corr.captPt, 3);
                        Cpy(sp.cutPt, sp.captPt, 3);

                        sp.catched = true;

                    }
                }
            } catch (...) {}

        }

#ifdef DEBUG
        if (!sp.catched) {
            std::cout << sp << std::endl;
        }
#endif

        assert(sp.catched);

    }

    // sucht nach gleichen captPts

    {
        std::map<Point3d, std::map<vtkIdType, std::vector<std::reference_wrapper<StripPt>>>> collapsed;

        PolyStripsType::iterator itr;
        StripPtsType::iterator itr2;

        for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
            PStrips &pStrips = itr->second;

            StripPtsType &pts = pStrips.pts;

            for (itr2 = pts.begin(); itr2 != pts.end(); ++itr2) {
                StripPt &sp = itr2->second;

                if (sp.capt & Capt::BOUNDARY) {
                    collapsed[{sp.cutPt[0], sp.cutPt[1], sp.cutPt[2]}][sp.ind].push_back(sp);
                }
            }
        }

        for (auto &[pt, map] : collapsed) {
            if (map.size() > 1) {
                // std::cout << pt << std::endl;

                // for (auto &[ind, pts] : map) {
                //     for (auto &p : pts) {
                //         std::cout << p << std::endl;
                //     }
                // }

                // if (map.size() == 2 && std::all_of(map.begin(), map.end(), [](const auto &m) { return m.second.size() == 1 && m.second.back().get().capt == Capt::EDGE; })) {

                //     for (auto &[ind, pts] : map) {
                //         StripPt &sp = pts.back();

                //         sp.capt = Capt::NOT;
                //         sp.edge[0] = sp.edge[1] = NOTSET;
                //         Cpy(sp.cutPt, sp.pt, 3);
                //     }

                // } else {
                //     return true;
                // }

                return true;

            }
        }
    }

    for (itr = polyLines.begin(); itr != polyLines.end(); ++itr) {
        PStrips &pStrips = polyStrips.at(itr->first);

        const IdsType &lines = itr->second;
        const StripPtsType &pts = pStrips.pts;

        StripsType &strips = pStrips.strips;

        // zusammensetzen

        std::deque<Pair> _lines;

        vtkIdList *linePts = vtkIdList::New();

        for (auto &i : lines) {
            contLines->GetCellPoints(i, linePts);
            _lines.emplace_back(linePts->GetId(0), linePts->GetId(1));
        }

        linePts->Delete();

        decltype(_lines)::iterator _itr;

        auto FindRight = [&pts, &_lines, &_itr](StripType &strip, const std::size_t &id) -> bool {
            auto &right = strip.back();

            if (pts.at(right.ind).capt == Capt::NOT) {
                for (_itr = _lines.begin(); _itr != _lines.end(); ++_itr) {
                    if (_itr->f == right.ind) {
                        strip.emplace_back(_itr->g, id);

                        _lines.erase(_itr);
                        return true;
                    } else if (_itr->g == right.ind) {
                        strip.emplace_back(_itr->f, id);

                        _lines.erase(_itr);
                        return true;
                    }
                }
            }

            return false;
        };

        auto FindLeft = [&pts, &_lines, &_itr](StripType &strip, const std::size_t &id) -> bool {
            auto &left = strip.front();

            if (pts.at(left.ind).capt == Capt::NOT) {
                for (_itr = _lines.begin(); _itr != _lines.end(); ++_itr) {
                    if (_itr->f == left.ind) {
                        strip.emplace_front(_itr->g, id);

                        _lines.erase(_itr);
                        return true;
                    } else if (_itr->g == left.ind) {
                        strip.emplace_front(_itr->f, id);

                        _lines.erase(_itr);
                        return true;
                    }
                }
            }

            return false;
        };

        std::size_t stripId {0};

        while (!_lines.empty()) {
            auto &last = _lines.back();

            StripType strip {{last.f, stripId}, {last.g, stripId}};
            _lines.pop_back();

            while (FindRight(strip, stripId)) {}
            while (FindLeft(strip, stripId)) {}

            strips.push_back(std::move(strip));

            stripId++;

        }

        CompleteStrips(pStrips);

    }

    // sucht nach schnitten zw. den strips

    {

        PolyStripsType::const_iterator itr;
        StripType::const_iterator itr2;

        for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
            const PStrips &pStrips = itr->second;

            const StripsType &strips = pStrips.strips;
            const StripPtsType &pts = pStrips.pts;
            const Base &base = pStrips.base;

            auto treePts = vtkSmartPointer<vtkPoints>::New();

            auto treePd = vtkSmartPointer<vtkPolyData>::New();
            treePd->Allocate(1);

            std::map<vtkIdType, vtkIdType> ptIds;
            
            double pt[2];

            for (const auto &p : pts) {
                Transform(p.second.pt, pt, base);

                ptIds.emplace(p.first, treePts->InsertNextPoint(pt[0], pt[1], 0));
            }

            for (const StripType &strip : strips) {
                for (itr2 = strip.begin(); itr2 != strip.end()-1; ++itr2) {

                    vtkIdList *line = vtkIdList::New();
                    line->InsertNextId(ptIds[itr2->ind]);
                    line->InsertNextId(ptIds[(itr2+1)->ind]);

                    treePd->InsertNextCell(VTK_LINE, line);

                    line->Delete();
                }
            }

            treePd->SetPoints(treePts);

            auto tree = vtkSmartPointer<vtkModifiedBSPTree>::New();
            tree->SetDataSet(treePd);
            tree->BuildLocator();

            vtkIdType numA, numB;
            const vtkIdType *lineA, *lineB;

            auto lineItr = vtk::TakeSmartPointer(treePd->GetLines()->NewIterator());

            for (lineItr->GoToFirstCell(); !lineItr->IsDoneWithTraversal(); lineItr->GoToNextCell()) {
                lineItr->GetCurrentCell(numA, lineA);

                double ptA[3], ptB[3];

                treePts->GetPoint(lineA[0], ptA);
                treePts->GetPoint(lineA[1], ptB);

                auto lineIds = vtkSmartPointer<vtkIdList>::New();

                tree->IntersectWithLine(ptA, ptB, 1e-5, nullptr, lineIds);

                for (vtkIdType i = 0; i < lineIds->GetNumberOfIds(); i++) {
                    treePd->GetCellPoints(lineIds->GetId(i), numB, lineB);

                    if (lineB[0] != lineA[0] && lineB[1] != lineA[0] && lineB[0] != lineA[1] && lineB[1] != lineA[1]) {
                        // schnitt gefunden

                        return true;
                    }
                }
            }

        }

    }

    return false;

}

void vtkPolyDataBooleanFilter::RemoveDuplicates (IdsType &lines) {

    typedef std::tuple<vtkIdType, vtkIdType, vtkIdType> LineType;

    std::vector<LineType> _lines;
    _lines.reserve(lines.size());

    auto line = vtkSmartPointer<vtkIdList>::New();

    vtkIdType a, b;

    for (const vtkIdType &id : lines) {
        contLines->GetCellPoints(id, line);

        a = line->GetId(0);
        b = line->GetId(1);

        if (std::find_if(_lines.begin(), _lines.end(), [&](const LineType &_line) {
            return (std::get<1>(_line) == a && std::get<2>(_line) == b)
                || (std::get<1>(_line) == b && std::get<2>(_line) == a);
        }) == _lines.end()) {
            _lines.emplace_back(id, a, b);
        }
    }

    if (_lines.size() != lines.size()) {
        lines.clear();

        for (const auto &_line : _lines) {
            lines.push_back(std::get<0>(_line));
        }

        lines.shrink_to_fit();
    }

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

bool vtkPolyDataBooleanFilter::CleanStrips () {
#ifdef DEBUG
    std::cout << "CleanStrips()" << std::endl;
#endif

    _IdsType inds;

    auto FindHoles = [&](PolyStripsType &polyStrips) {
        PolyStripsType::iterator itr;

        for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
            PStrips &pStrips = itr->second;
            StripsType &strips = pStrips.strips;
            StripPtsType &pts = pStrips.pts;

            strips.erase(std::remove_if(strips.begin(), strips.end(), [&](const StripType &strip) {
                if (pts.at(strip.front().ind).capt == Capt::NOT
                    && pts.at(strip.back().ind).capt == Capt::NOT
                    && !HasArea(strip)) {

                    for (const StripPtR &p : strip) {
                        inds.emplace(p.ind);
                    }

                    return true;
                }

                return false;
            }), strips.end());
        }
    };

    FindHoles(polyStripsA);
    FindHoles(polyStripsB);

#ifdef DEBUG
    std::cout << "inds: [";
    for (auto &ind : inds) {
        std::cout << ind << ", ";
    }
    std::cout << "]" << std::endl;
#endif

    auto CleanOther = [&](PolyStripsType &polyStrips) {
        PolyStripsType::iterator itr;

        for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
            PStrips &pStrips = itr->second;
            StripsType &strips = pStrips.strips;

            strips.erase(std::remove_if(strips.begin(), strips.end(), [&](const StripType &strip) {
                auto found = std::find_if(strip.begin(), strip.end(), [&](const StripPtR &p) {
                    return inds.find(p.ind) != inds.end();
                });

                return found != strip.end();

            }), strips.end());
        }
    };

    CleanOther(polyStripsA);
    CleanOther(polyStripsB);

    auto lines = vtkSmartPointer<vtkIdList>::New();

    vtkIdType i, j, numLines;

    for (vtkIdType ind : inds) {
        contLines->GetPointCells(ind, lines);
        numLines = lines->GetNumberOfIds();

        for (i = 0; i < numLines; i++) {
            contLines->DeleteCell(lines->GetId(i));
        }
    }

    j = 0;

    numLines = contLines->GetNumberOfCells();

    for (i = 0; i < numLines; i++) {
        if (contLines->GetCellType(i) == VTK_EMPTY_CELL) {
            j++;
        }
    }

    if (j == numLines) {
        return true;
    }

    return false;

}

template<typename _RefsType>
void ComputeNormal (const StripPtsType &pts, const _RefsType &poly, double *n) {
    n[0] = 0; n[1] = 0; n[2] = 0;

    typename _RefsType::const_iterator itrA, itrB;

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

void CleanPoly (vtkPolyData *pd, IdsType &poly) {
    IdsType newPoly;
    newPoly.reserve(poly.size());

    double pt[3];

    std::map<vtkIdType, Point3d> _pts;

    for (vtkIdType id : poly) {
        pd->GetPoint(id, pt);

        _pts.emplace(std::piecewise_construct,
            std::forward_as_tuple(id),
            std::forward_as_tuple(pt[0], pt[1], pt[2]));
    }

    IdsType::const_iterator itrA, itrB;

    for (itrA = poly.begin(); itrA != poly.end(); ++itrA) {
        itrB = itrA+1;
        if (itrB == poly.end()) {
            itrB = poly.begin();
        }

        auto _a = _pts.find(*itrA);
        auto _b = _pts.find(*itrB);

        if (_a->second == _b->second) {} else {
            newPoly.push_back(*itrA);
        }
    }

    newPoly.shrink_to_fit();

    poly.swap(newPoly);

}

bool vtkPolyDataBooleanFilter::CutCells (vtkPolyData *pd, PolyStripsType &polyStrips) {
#ifdef DEBUG
    std::cout << "CutCells()" << std::endl;
#endif

    vtkPoints *pdPts = pd->GetPoints();

    vtkIdTypeArray *origCellIds = vtkIdTypeArray::SafeDownCast(pd->GetCellData()->GetScalars("OrigCellIds"));

    PolyStripsType::iterator itrA;

    for (itrA = polyStrips.begin(); itrA != polyStrips.end(); ++itrA) {
        const vtkIdType &polyInd = itrA->first;
        PStrips &pStrips = itrA->second;

        StripsType &strips = pStrips.strips;
        StripPtsType &pts = pStrips.pts;

        IdsType &poly = pStrips.poly;

        vtkIdType origId = origCellIds->GetValue(polyInd);

        if (std::all_of(pts.begin(), pts.end(), [](const auto &p) { return p.second.capt & Capt::A || p.second.capt & Capt::B; })) {
            Poly _poly;

            for (auto &id : poly) {
                auto pt = pd->GetPoint(id);
                _poly.emplace_back(pt[0], pt[1], pt[2]);
            }

            std::set<Point3d> ptsA(_poly.begin(), _poly.end()), ptsB;

            for (const auto& [ind, sp] : pts) {
                ptsB.emplace(sp.cutPt[0], sp.cutPt[1], sp.cutPt[2]);
            }

            if (ptsA == ptsB) {
                vtkIdList *cell = vtkIdList::New();

                for (auto &p : _poly) {
                    cell->InsertNextId(pdPts->InsertNextPoint(p.x, p.y, p.z));
                }

                pd->InsertNextCell(VTK_POLYGON, cell);
                origCellIds->InsertNextValue(origId);

                cell->Delete();

                pd->DeleteCell(polyInd);

                continue;
            }

        }

        double _t = 0;
        std::map<vtkIdType, double> absoluteT;

        for (auto &id : poly) {
            absoluteT.emplace(id, _t++);
        }

#ifdef DEBUG
        std::cout << "polyInd " << polyInd << ", poly [";
        for (auto &id : poly) {
            std::cout << id << ", ";
        }
        std::cout << "]" << std::endl;
#endif

        // alle strips gültig?

        if (std::find_if(strips.begin(), strips.end(), [&](const StripType &s) {
            return pts[s.front().ind].capt == Capt::BRANCHED && pts[s.back().ind].capt == Capt::BRANCHED;
        }) != strips.end()) {
            return true;
        }

        // holes sichern
        StripsType holes;

        auto fct = [&](const StripType &s) {
            return pts[s.front().ind].capt == Capt::NOT && pts[s.back().ind].capt == Capt::NOT;
        };

        std::copy_if(strips.begin(), strips.end(), std::back_inserter(holes), fct);

        strips.erase(std::remove_if(strips.begin(), strips.end(), fct), strips.end());

        // init

        for (auto &strip : strips) {
#ifdef DEBUG
            std::cout << "strip [";
            for (auto &p : strip) {
                std::cout << p.ind << ", ";
            }
            std::cout << "]" << std::endl;
#endif

            // enden auf gleichem edge
            if (pts[strip.front().ind].edge[0] == pts[strip.back().ind].edge[0]
                && strip.front().ind != strip.back().ind
                && pts[strip.front().ind].t > pts[strip.back().ind].t) {

                std::reverse(strip.begin(), strip.end());
            }

            // branched strip
            if (pts[strip.front().ind].capt == Capt::BRANCHED
                && pts[strip.back().ind].capt & Capt::BOUNDARY) {

                std::reverse(strip.begin(), strip.end());
            }

            StripPt &start = pts[strip.front().ind],
                &end = pts[strip.back().ind];

            strip.front().side = Side::START;
            strip.front().ref = start.edge[0];

            if (end.capt & Capt::BOUNDARY) {
                strip.back().side = Side::END;
                strip.back().ref = end.edge[0];
            }

            for (auto &p : strip) {
                StripPt &sp = pts[p.ind];

                p.desc[0] = pdPts->InsertNextPoint(sp.cutPt);
                p.desc[1] = pdPts->InsertNextPoint(sp.cutPt);

#ifdef DEBUG
                std::cout << sp << " => " << p << std::endl;
#endif

            }
        }

        std::deque<IdsType> polys;
        polys.push_back(poly);

        // gruppiert die branched strips

        std::map<vtkIdType, _StripsType> groups;

        for (auto &strip : strips) {
            if (pts[strip.back().ind].capt == Capt::BRANCHED) {
                groups[strip.back().ind].emplace_back(strip);
            }
        }

        std::vector<std::size_t> assembled;

        double n[3], ang, pt[3], proj[2];

        decltype(groups)::iterator itrB;

        for (itrB = groups.begin(); itrB != groups.end(); ++itrB) {
            _StripsType &_strips = itrB->second;

            // sortiert die strips

            std::sort(_strips.begin(), _strips.end(), [&](const StripType &a, const StripType &b) {
                if (a.front().ind == b.front().ind) {
                    ConstRefsType poly_(b.begin(), b.end());
                    poly_.insert(poly_.end(), a.rbegin(), a.rend());

                    ComputeNormal(pts, poly_, n);

                    ang = vtkMath::Dot(pStrips.n, n);

                    return ang > .999999;

                } else {
                    const StripPt &pA = pts[a.front().ind],
                        &pB = pts[b.front().ind];

                    return absoluteT[pA.edge[0]]+pA.t < absoluteT[pB.edge[0]]+pB.t;
                }
            });

            // die klassische sanduhr

            auto next = std::find_if(polys.begin(), polys.end(), [&_strips](const IdsType &p) {
                return std::find(p.begin(), p.end(), _strips.front().get().front().ref) != p.end();
            });

            assert(next != polys.end());

            std::for_each(_strips.begin(), _strips.end(), [&assembled](const StripType &s) {
                assembled.push_back(s.front().strip);
            });

            std::vector<IdsType> newPolys;
            newPolys.reserve(_strips.size()+1);

            _StripsType::const_iterator _itrA, _itrB;

            StripType::const_iterator _itrC;
            StripType::const_reverse_iterator _itrD;

            for (_itrA = _strips.begin(); _itrA != _strips.end(); ++_itrA) {
                _itrB = _itrA+1;

                if (_itrB == _strips.end()) {
                    _itrB = _strips.begin();
                }

                const StripType &stripA = *_itrA,
                    &stripB = *_itrB;

                IdsType newPoly;

                for (_itrC = stripB.begin(); _itrC != stripB.end(); ++_itrC) {
                    newPoly.push_back(_itrC->desc[0]);
                }

                for (_itrD = stripA.rbegin()+1; _itrD != stripA.rend(); ++_itrD) {
                    newPoly.push_back(_itrD->desc[1]);
                }

                // punkte zw. den enden einfügen

                if (stripA.front().ref != stripB.front().ref) {
                    auto posA = std::find(next->begin(), next->end(), stripA.front().ref);
                    auto posB = std::find(next->begin(), next->end(), stripB.front().ref);

                    for (;;) {
                        posA++;

                        if (posA == next->end()) {
                            posA = next->begin();
                        }

                        newPoly.push_back(*posA);

                        if (posA == posB) {
                            break;
                        }
                    }
                }

                CleanPoly(pd, newPoly);

                Poly _poly;

                for (vtkIdType &id : newPoly) {
                    pd->GetPoint(id, pt);
                    Transform(pt, proj, pStrips.base);
                    _poly.emplace_back(proj[0], proj[1], 0);
                }

                // refs aktualisieren

                const StripPt &pA = pts[stripA.front().ind],
                    &pB = pts[stripB.front().ind];

                for (auto &s : strips) {
                    if (std::find(assembled.begin(), assembled.end(), s.front().strip) != assembled.end()) {
                        continue;
                    }

                    // noch nicht eingebaut

                    const StripPt &endA = pts[s.front().ind],
                        &endB = pts[s.back().ind];

                    if (endA.capt & Capt::BOUNDARY
                        && pA.edge[0] == endA.edge[0]
                        && endA.t > pA.t
                        && (pA.edge[0] != pB.edge[0] || endA.t < pB.t)) {
                        s.front().ref = stripA.front().desc[1];

                        if (endB.ind == pA.ind) {
                            s.back().ref = stripA.front().desc[1];
                        } else if (endB.ind == pB.ind) {
                            s.back().ref = stripB.front().desc[0];
                        }
                    }

                    if (endB.capt & Capt::BOUNDARY
                        && pA.edge[0] == endB.edge[0]
                        && endB.t > pA.t
                        && (pA.edge[0] != pB.edge[0] || endB.t < pB.t)) {
                        s.back().ref = stripA.front().desc[1];

                        if (endA.ind == pA.ind) {
                            s.front().ref = stripA.front().desc[1];
                        } else if (endA.ind == pB.ind) {
                            s.front().ref = stripB.front().desc[0];
                        }
                    }

                    if (endA.ind == pA.ind && endB.ind == pB.ind) {
                        s.front().ref = stripA.front().desc[1];
                        s.back().ref = stripB.front().desc[0];
                    } else if (endB.ind == pA.ind && endA.ind == pB.ind) {
                        s.back().ref = stripA.front().desc[1];
                        s.front().ref = stripB.front().desc[0];
                    }

                    if (endB.capt == Capt::BRANCHED) {
                        Transform(endB.pt, proj, pStrips.base);

                        if (PointInPoly(_poly, {proj[0], proj[1], 0})) {
                            if (endA.ind == pA.ind) {
                                s.front().ref = stripA.front().desc[1];
                            } else if (endA.ind == pB.ind) {
                                s.front().ref = stripB.front().desc[0];
                            }
                        }
                    }

                }

                newPolys.push_back(std::move(newPoly));

            }

            polys.erase(next);

            polys.insert(polys.end(), newPolys.begin(), newPolys.end());

        }

        // restliche strips einbauen

        std::vector<IdsType> newPolys;

        for (auto &next : polys) {
            _StripsType _strips;

            for (auto &strip : strips) {
                if (pts[strip.back().ind].capt != Capt::BRANCHED
                    && std::find(next.begin(), next.end(), strip.front().ref) != next.end()) {

                    _strips.emplace_back(strip);
                }
            }

            if (_strips.empty()) {
                newPolys.push_back(next);
                continue;
            }

            std::deque<IdsType> _newPolys {next};

            std::map<vtkIdType, RefsType> edges;

            for (auto &s : _strips) {
                StripType &strip = s;

                const StripPt &a = pts[strip.front().ind],
                    &b = pts[strip.back().ind];

                edges[a.edge[0]].push_back(std::ref(strip.front()));
                edges[b.edge[0]].push_back(std::ref(strip.back()));
            }

            // sortiert die punkte auf den kanten

            double n[3], ang, r, rA, rB;

            decltype(edges)::iterator itrC;

            for (itrC = edges.begin(); itrC != edges.end(); ++itrC) {
                const vtkIdType &id = itrC->first;
                RefsType &edge = itrC->second;

#ifdef DEBUG
                std::cout << "edge (" << id << ", _)" << std::endl;
#endif

                std::sort(edge.begin(), edge.end(), [&](const StripPtR &a, const StripPtR &b) {
                    const StripPt &a_ = pts[a.ind],
                        &b_ = pts[b.ind];

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

                            if (eA_.ind != eB_.ind) {
                                r = absoluteT[id]+a_.t;
                                rA = absoluteT[eA_.edge[0]]+eA_.t;
                                rB = absoluteT[eB_.edge[0]]+eB_.t;

                                rA = rA > r ? rA-r : rA+_t-r,
                                rB = rB > r ? rB-r : rB+_t-r;

                                return rB < rA;

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

                                ComputeNormal(pts, poly_, n);
                                ang = vtkMath::Dot(pStrips.n, n);

                                return ang < .999999;

                            }
                        } else {
                            // gleicher strip

                            StripType &strip = strips[a.strip];

                            if (HasArea(strip)) {
                                RefsType poly_(strip.begin(), strip.end()-1);

                                ComputeNormal(pts, poly_, n);
                                ang = vtkMath::Dot(pStrips.n, n);

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
                for (auto& e : edge) {
                    std::cout << e << std::endl;
                }
#endif

            }

            // strips einbauen

            IdsType::iterator _itrA;
            StripType::reverse_iterator _itrB;

            for (auto &s : _strips) {
                StripType &strip = s;

                const StripPtR &start = strip.front(),
                    &end = strip.back();

#ifdef DEBUG
                std::cout << "strip " << start.strip
                    << " , refs (" << start.ref << ", " << end.ref << ")"
                    << std::endl;
#endif

                std::size_t cycle = 0;

                while (true) {
                    if (cycle == _newPolys.size()) {
                        break;
                    }

                    IdsType _next = _newPolys.front();
                    _newPolys.pop_front();

                    std::vector<IdsType> splitted(2);

                    if (std::find(_next.begin(), _next.end(), start.ref) != _next.end()) {
                        if (start.ref == end.ref) {
                            for (_itrA = _next.begin(); _itrA != _next.end(); ++_itrA) {
                                splitted[0].push_back(*_itrA);

#ifdef DEBUG
                                std::cout << "adding " << *_itrA << " to 0" << std::endl;
#endif

                                if (*_itrA == start.ref) {
                                    for (auto &p : strip) {
                                        splitted[0].push_back(p.desc[0]);

#ifdef DEBUG
                                        std::cout << "adding " << p.desc[0] << " to 0" << std::endl;
#endif

                                    }
                                }
                            }

                            // strip selbst ist ein polygon

                            for (_itrB = strip.rbegin(); _itrB != strip.rend(); ++_itrB) {
                                splitted[1].push_back(_itrB->desc[1]);

#ifdef DEBUG
                                std::cout << "adding " << _itrB->desc[1] << " to 1" << std::endl;
#endif

                            }

                        } else {
                            std::size_t curr = 0;

                            for (_itrA = _next.begin(); _itrA != _next.end(); ++_itrA) {
                                IdsType &poly = splitted[curr];

                                poly.push_back(*_itrA);

#ifdef DEBUG
                                std::cout << "adding " << *_itrA << " to " << curr << std::endl;
#endif

                                if (*_itrA == start.ref) {
                                    for (auto &p : strip) {
                                        poly.push_back(p.desc[0]);

#ifdef DEBUG
                                        std::cout << "adding " << p.desc[0] << " to " << curr << std::endl;
#endif

                                    }

                                    curr = curr == 0 ? 1 : 0;

                                } else if (*_itrA == end.ref) {
                                    for (_itrB = strip.rbegin(); _itrB != strip.rend(); ++_itrB) {
                                        poly.push_back(_itrB->desc[1]);

#ifdef DEBUG
                                        std::cout << "adding " << _itrB->desc[1] << " to " << curr << std::endl;
#endif

                                    }

                                    curr = curr == 0 ? 1 : 0;
                                }

                            }
                        }
                    }

                    if (!splitted[1].empty()) {
                        // refs aktualisieren

                        for (itrC = edges.begin(); itrC != edges.end(); ++itrC) {
                            RefsType &edge = itrC->second;

                            RefsType::iterator _itrA;

                            for (_itrA = edge.begin()+1; _itrA != edge.end(); ++_itrA) {
                                StripPtR &sp = *_itrA;

                                if (sp.strip > start.strip) {
#ifdef DEBUG
                                    std::cout << "ind " << sp.ind << ", strip " << sp.strip << std::endl;
#endif

                                    RefsType::const_reverse_iterator _itrB(_itrA);

                                    std::shared_ptr<StripPtR> _p;

                                    for (; _itrB != edge.rend(); ++_itrB) {
                                        const StripPtR &p = *_itrB;

                                        if (p.strip != sp.strip) {
                                            if (p.strip <= start.strip) {
#ifdef DEBUG
                                                std::cout << "*1 ref " << sp.ref;
#endif

                                                if (p.side == Side::END) {
                                                    sp.ref = p.desc[0];
                                                } else {
                                                    sp.ref = p.desc[1];
                                                }

#ifdef DEBUG
                                                std::cout << " -> " << sp.ref << " (from strip " << p.strip << ", ind " << p.ind << ")" << std::endl;
#endif

                                                _p = std::make_shared<StripPtR>(p);

                                                break;

                                            }
                                        } else {
#ifdef DEBUG
                                            std::cout << "*2 ref " << sp.ref << " -> " << p.ref << " (from strip " << p.strip << ", ind " << p.ind << ")" << std::endl;
#endif

                                            sp.ref = p.ref;
                                            break;
                                        }
                                    }

                                    RefsType::const_iterator _itrC(_itrA);

                                    ++_itrC;

                                    for (; _itrC != edge.end(); ++_itrC) {
                                        const StripPtR &p = *_itrC;

                                        if (p.ind != sp.ind) {
                                            break;
                                        }

                                        if (p.strip <= start.strip) {
                                            if (_p && p.ind == _p->ind && p.strip < _p->strip) {
                                                break;
                                            }

#ifdef DEBUG
                                            std::cout << "*3 ref " << sp.ref;
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

                            // sonderfall
                            if (edge.size() > 1) {
                                StripPtR &a = edge.front(),
                                    &b = *(edge.begin()+1);

                                if (a.ind == b.ind
                                    && b.strip == start.strip
                                    && pts[a.ind].capt == Capt::A) { // sollte weg

#ifdef DEBUG
                                    std::cout << "*4 ref " << a.ref;
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

                        for (auto &p : splitted) {
                            CleanPoly(pd, p);
                        }

                        // prüft, ob die erstellten _newPolys gültig sind

                        if (splitted[0].size() > 2) {
                            _newPolys.push_back(splitted[0]);
                        }

                        if (HasArea(strip) && splitted[1].size() > 2) {
                            _newPolys.push_back(splitted[1]);
                        }

                        cycle = 0;

                        break;

                    } else {
                        _newPolys.push_back(_next);

                        cycle++;
                    }

                }

            }

            newPolys.insert(newPolys.end(), _newPolys.begin(), _newPolys.end());
        }

        // erzeugte polys hinzufügen

        IdsType descIds;
        descIds.reserve(newPolys.size());

        for (auto &p : newPolys) {
            vtkIdList *cell = vtkIdList::New();

            for (vtkIdType &id : p) {
                cell->InsertNextId(id);
            }

            descIds.emplace_back(pd->InsertNextCell(VTK_POLYGON, cell));
            origCellIds->InsertNextValue(origId);

            cell->Delete();
        }

        pd->DeleteCell(polyInd);

        // holes einbauen
        if (!holes.empty()) {
            try {
                Merger(pd, pStrips, holes, descIds, origId).Run();
            } catch (const std::runtime_error &e) {
                return true;
            }
        }

    }

    pd->RemoveDeletedCells();
    pd->BuildCells();

    return false;

}

void vtkPolyDataBooleanFilter::RestoreOrigPoints (vtkPolyData *pd, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "RestoreOrigPoints()" << std::endl;
#endif

    pd->DeleteLinks(); pd->BuildLinks();

    vtkKdTreePointLocator *loc = vtkKdTreePointLocator::New();
    loc->SetDataSet(pd);
    loc->BuildLocator();

    PolyStripsType::const_iterator itr;
    StripPtsType::const_iterator itr2;

    auto pts = vtkSmartPointer<vtkIdList>::New();
    vtkIdType i, numPts;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
        const PStrips &pStrips = itr->second;

        for (itr2 = pStrips.pts.begin(); itr2 != pStrips.pts.end(); ++itr2) {
            const StripPt &sp = itr2->second;

            if (sp.capt & Capt::BOUNDARY) {
                FindPoints(loc, sp.cutPt, pts);
                numPts = pts->GetNumberOfIds();

                for (i = 0; i < numPts; i++) {
                    pd->GetPoints()->SetPoint(pts->GetId(i), sp.pt);
                }

            }

        }
    }

    loc->FreeSearchStructure();
    loc->Delete();

}

void vtkPolyDataBooleanFilter::DisjoinPolys (vtkPolyData *pd, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "DisjoinPolys()" << std::endl;
#endif

    pd->DeleteLinks(); pd->BuildLinks();

    vtkKdTreePointLocator *loc = vtkKdTreePointLocator::New();
    loc->SetDataSet(pd);

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

void vtkPolyDataBooleanFilter::ResolveOverlaps (vtkPolyData *pd, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "ResolveOverlaps()" << std::endl;
#endif

    typedef std::pair<std::reference_wrapper<const StripPtsType>, std::reference_wrapper<const StripPt>> _Type;

    std::map<vtkIdType, std::vector<_Type>> _pts;

    PolyStripsType::const_iterator itr;
    StripPtsType::const_iterator itr2;

    for (itr = polyStrips.begin(); itr != polyStrips.end(); ++itr) {
        const PStrips &pStrips = itr->second;

        const StripPtsType &pts = pStrips.pts;

        for (itr2 = pts.begin(); itr2 != pts.end(); ++itr2) {
            const StripPt &sp = itr2->second;

            if (sp.capt == Capt::EDGE) {
                _pts[sp.ind].emplace_back(pts, sp);
            }
        }
    }

    decltype(_pts)::iterator itr3;

    StripPtsType::const_iterator itr4;

    for (itr3 = _pts.begin(); itr3 != _pts.end(); ++itr3) {
        auto &pairs = itr3->second;

        if (pairs.size() == 2) {
            auto &pairA = pairs[0];
            auto &pairB = pairs[1];

            if (pairA.second.get().edge[1] != pairB.second.get().edge[0]) {
                pairA.swap(pairB);
            }

            auto edgeA = pairA.second.get().edge;
            auto edgeB = pairB.second.get().edge;

            assert(edgeA[1] == edgeB[0]);

            if (edgeA[1] == edgeB[0] && edgeA[0] != edgeB[1]) {

#ifdef DEBUG
                std::cout << itr3->first << ": "
                    << edgeA[0] << ", "
                    << edgeA[1] << ", "
                    << edgeB[0] << ", "
                    << edgeB[1] << std::endl;
#endif

                auto &ptsA = pairA.first.get();
                auto &ptsB = pairB.first.get();

                std::deque<std::reference_wrapper<const StripPt>> _edgeA, _edgeB;

                for (itr4 = ptsA.begin(); itr4 != ptsA.end(); ++itr4) {
                    auto &sp = itr4->second;
                    if (sp.edge[0] == edgeA[0] && sp.edge[1] == edgeA[1]) {
                        _edgeA.emplace_back(sp);
                    }
                }

                for (itr4 = ptsB.begin(); itr4 != ptsB.end(); ++itr4) {
                    auto &sp = itr4->second;
                    if (sp.edge[0] == edgeB[0] && sp.edge[1] == edgeB[1]) {
                        _edgeB.emplace_back(sp);
                    }
                }

                std::sort(_edgeA.begin(), _edgeA.end(), [](const StripPt &l, const StripPt &r) { return l.t < r.t; });
                std::sort(_edgeB.begin(), _edgeB.end(), [](const StripPt &l, const StripPt &r) { return l.t < r.t; });

#ifdef DEBUG
                for (const StripPt &sp : _edgeA) {
                    std::cout << sp << std::endl;
                }

                for (const StripPt &sp : _edgeB) {
                    std::cout << sp << std::endl;
                }
#endif

                assert(_edgeA.back().get().ind == itr3->first);
                assert(_edgeB.front().get().ind == itr3->first);

                _edgeA.pop_back();
                _edgeB.pop_front();

                double pA[3], pB[3];

                if (_edgeA.empty()) {
                    pd->GetPoint(edgeA[0], pA);
                } else {
                    Cpy(pA, _edgeA.back().get().pt, 3);
                }

                if (_edgeB.empty()) {
                    pd->GetPoint(edgeB[1], pB);
                } else {
                    Cpy(pB, _edgeB.front().get().pt, 3);
                }

                Point3d a {pA[0], pA[1], pA[2]};
                Point3d b {pB[0], pB[1], pB[2]};

                auto cells = vtkSmartPointer<vtkIdList>::New(),
                    cell = vtkSmartPointer<vtkIdList>::New();

                pd->GetPointCells(edgeA[1], cells);
                vtkIdType i, j, numCells = cells->GetNumberOfIds();

                double pt[3];

                vtkIdType id;

                for (i = 0; i < numCells; i++) {
                    pd->GetCellPoints(cells->GetId(i), cell);
                    Poly poly;
                    for (j = 0; j < cell->GetNumberOfIds(); j++) {
                        pd->GetPoint(cell->GetId(j), pt);
                        poly.emplace_back(pt[0], pt[1], pt[2]);
                    }
                    if (std::find(poly.begin(), poly.end(), a) != poly.end()
                        && std::find(poly.begin(), poly.end(), b) != poly.end()) {

                        contLines->GetPoint(itr3->first, pt);

                        id = pd->GetPoints()->InsertNextPoint(pt);

#ifdef DEBUG
                        std::cout << cells->GetId(i)
                            << ": " << edgeA[1] << " -> " << id
                            << std::endl;
#endif

                        pd->ReplaceCellPoint(cells->GetId(i), edgeA[1], id);

                        break;

                    }

                }

            }

        }
    }

}

void vtkPolyDataBooleanFilter::AddAdjacentPoints (vtkPolyData *pd, vtkIdTypeArray *conts, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "AddAdjacentPoints()" << std::endl;
#endif

    pd->DeleteLinks(); pd->BuildLinks();

    vtkIdTypeArray *origCellIds = vtkIdTypeArray::SafeDownCast(pd->GetCellData()->GetScalars("OrigCellIds"));

    struct Cmp {
        bool operator() (const StripPt &l, const StripPt &r) const {
            return l.t < r.t;
        }
    };

    auto loc = vtkSmartPointer<vtkKdTreePointLocator>::New();
    loc->SetDataSet(pd);
    loc->BuildLocator();

    auto lines = vtkSmartPointer<vtkIdList>::New();

    vtkIdType i, j, numLines;

    auto ptsA = vtkSmartPointer<vtkIdList>::New();
    auto ptsB = vtkSmartPointer<vtkIdList>::New();

    vtkIdType numPtsA, numPtsB;

    auto cells = vtkSmartPointer<vtkIdList>::New();

    vtkIdType numCells;

    auto poly = vtkSmartPointer<vtkIdList>::New();
    auto newPoly = vtkSmartPointer<vtkIdList>::New();

    vtkIdType numPts;

    vtkIdType idA, idB;

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
                    contLines->GetPointCells(itrB->get().ind, lines);
                    numLines = lines->GetNumberOfIds();

                    _IdsType involved;

                    for (i = 0; i < numLines; i++) {
                        involved.emplace(conts->GetValue(lines->GetId(i)));
                    }

                    if (involved.size() > 1) {
                        break;
                    }

                    ++itrB;
                }

                if (itrA+1 != itrB) {
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

                                for (i = 0; i < numPts; i++) {
                                    newPoly->InsertNextId(poly->GetId(i));

                                    idA = poly->GetId(i);
                                    idB = i+1 == numPts ? poly->GetId(0) : poly->GetId(i+1);

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

    loc->FreeSearchStructure();

    pd->RemoveDeletedCells();

}

void vtkPolyDataBooleanFilter::MergePoints (vtkPolyData *pd, PolyStripsType &polyStrips) {

#ifdef DEBUG
    std::cout << "MergePoints()" << std::endl;
#endif

    pd->BuildCells();
    pd->DeleteLinks(); pd->BuildLinks();

    contLines->DeleteLinks(); contLines->BuildLinks();

    auto loc = vtkSmartPointer<vtkKdTreePointLocator>::New();
    loc->SetDataSet(pd);

    PolyStripsType::const_iterator itr;
    StripsType::const_iterator itr2;

    std::map<vtkIdType, _IdsType> neighPts;

    auto pts = vtkSmartPointer<vtkIdList>::New();

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

    auto polys = vtkSmartPointer<vtkIdList>::New();
    auto poly = vtkSmartPointer<vtkIdList>::New();

    vtkIdType ind, polyId, _numPts, before, after;

    decltype(neighPts)::const_iterator itr3;

    for (itr3 = neighPts.begin(); itr3 != neighPts.end(); ++itr3) {
        const auto &inds = itr3->second;

        std::map<Point3d, std::vector<Pair>> pairs;

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

        std::deque<std::deque<std::reference_wrapper<const Pair>>> Pairs;

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
                    pd->ReplaceCellPoint(itr6->get().f, itr6->get().g, group.front().get().g);
                }

                group.clear();
            }
        }

    }

    loc->FreeSearchStructure();

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

    auto filterdA = vtkSmartPointer<vtkPolyData>::New();
    filterdA->DeepCopy(modPdA);

    auto filterdB = vtkSmartPointer<vtkPolyData>::New();
    filterdB->DeepCopy(modPdB);

    // ungenutzte punkte löschen
    auto cleanA = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanA->PointMergingOff();
    cleanA->SetInputData(filterdA);

    auto cleanB = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanB->PointMergingOff();
    cleanB->SetInputData(filterdB);

    // regionen mit skalaren ausstatten
    auto cfA = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    cfA->SetExtractionModeToAllRegions();
    cfA->ColorRegionsOn();
    cfA->SetInputConnection(cleanA->GetOutputPort());

    auto cfB = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
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

    if (OperMode == OPER_NONE) {
        resultA->ShallowCopy(pdA);
        resultB->ShallowCopy(pdB);

        contLines->RemoveDeletedCells();
        resultC->ShallowCopy(contLines);

        return false;
    }

    auto plA = vtkSmartPointer<vtkKdTreePointLocator>::New();
    plA->SetDataSet(pdA);
    plA->BuildLocator();

    auto plB = vtkSmartPointer<vtkKdTreePointLocator>::New();
    plB->SetDataSet(pdB);
    plB->BuildLocator();

    pdA->DeleteLinks(); pdA->BuildLinks();
    pdB->DeleteLinks(); pdB->BuildLinks();

    vtkIdTypeArray *scalarsA = vtkIdTypeArray::SafeDownCast(pdA->GetPointData()->GetScalars());
    vtkIdTypeArray *scalarsB = vtkIdTypeArray::SafeDownCast(pdB->GetPointData()->GetScalars());

    auto line = vtkSmartPointer<vtkIdList>::New();

    double ptA[3], ptB[3];

    auto fptsA = vtkSmartPointer<vtkIdList>::New();
    auto lptsA = vtkSmartPointer<vtkIdList>::New();

    auto fptsB = vtkSmartPointer<vtkIdList>::New();
    auto lptsB = vtkSmartPointer<vtkIdList>::New();

    std::map<vtkIdType, Loc> locsA, locsB;

    vtkIdType i, j, numLines = contLines->GetNumberOfCells();

    IdsType _failed;

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
            _failed.push_back(i);

            // return true;
        }

    }

    if (_failed.size() > 0) {
#ifdef DEBUG
        for (auto i : _failed) {
            std::cout << "failed at " << i
                << std::endl;
        }
#endif

        return true;
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
        vtkCellIterator *cellItr;

        vtkIdType cellId;
        vtkIdList *ptIds;

        if (comb[0] == Loc::INSIDE) {
            cellItr = regsA->NewCellIterator();

            for (cellItr->InitTraversal(); !cellItr->IsDoneWithTraversal(); cellItr->GoToNextCell()) {
                cellId = cellItr->GetCellId();
                ptIds = cellItr->GetPointIds();

                if (locsA.count(scalarsA->GetValue(ptIds->GetId(0))) == 1) {
                    regsA->ReverseCell(cellId);
                }
            }

            cellItr->Delete();
        }

        if (comb[1] == Loc::INSIDE) {
            cellItr = regsB->NewCellIterator();

            for (cellItr->InitTraversal(); !cellItr->IsDoneWithTraversal(); cellItr->GoToNextCell()) {
                cellId = cellItr->GetCellId();
                ptIds = cellItr->GetPointIds();

                if (locsB.count(scalarsB->GetValue(ptIds->GetId(0))) == 1) {
                    regsB->ReverseCell(cellId);
                }
            }

            cellItr->Delete();
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

    // resultB bleibt hier leer

    resultA->ShallowCopy(cfPd);

    resultA->GetCellData()->AddArray(newOrigCellIdsA);
    resultA->GetCellData()->AddArray(newOrigCellIdsB);

    contLines->RemoveDeletedCells();
    resultC->ShallowCopy(contLines);

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

Merger::Merger (vtkPolyData *pd, const PStrips &pStrips, const StripsType &strips, const IdsType &descIds, vtkIdType origId) : pd(pd), pStrips(pStrips), origId(origId) {

    const StripPtsType &pts = pStrips.pts;
    const Base &base = pStrips.base;

    vtkIdType i, num;
    const vtkIdType *cell;

    double pt[3];

    for (auto id : descIds) {
        pd->GetCellPoints(id, num, cell);

        Poly p;

        for (i = 0; i < num; i++) {
            pd->GetPoint(cell[i], pt);

            double proj[2];
            Transform(pt, proj, base);

            p.emplace_back(proj[0], proj[1], 0, cell[i]);

        }

        polys.push_back(p);

        pd->DeleteCell(id);
    }

    for (auto &strip : strips) {
        Poly p;

        for (auto &sp : strip) {
            const double *pt = pts.at(sp.ind).pt;

            double proj[2];
            Transform(pt, proj, base);

            p.emplace_back(proj[0], proj[1], 0);
        }

        p.pop_back();

        double n[3];
        ComputeNormal(p, n);

        if (n[2] < 0) {
            Poly q(p.rbegin(), p.rend());
            p.swap(q);
        }

        innerIds.push_back(polys.size());

        polys.push_back(p);
    }

}

void Merger::Run () {
    // mergen mit hilfe von vtkKdTree und vtkModifiedBSPTree

    vtkPoints *pdPts = pd->GetPoints();
    vtkIdTypeArray *origCellIds = vtkIdTypeArray::SafeDownCast(pd->GetCellData()->GetScalars("OrigCellIds"));

    assert(origCellIds != nullptr);

    const Base &base = pStrips.base;

    std::vector<GroupType> groups(polys.size());

    PolysType::const_iterator itrA, itrB;

    std::size_t i {0};

    for (itrA = polys.begin(); itrA != polys.end(); ++itrA) {
        if (std::find(innerIds.begin(), innerIds.end(), i) != innerIds.end()) {
            std::size_t j {0};
            for (itrB = polys.begin(); itrB != polys.end(); ++itrB) {
                if (itrA != itrB && PointInPoly(*itrB, *itrA->begin())) {
                    groups[j].push_back(i);
                }
                j++;
            }
        }
        i++;
    }

    std::size_t parent = 0;

    for (auto &group : groups) {
        GroupType parents;

        for (auto &index : group) {
            const GroupType &_group = groups[index];
            parents.insert(parents.end(), _group.begin(), _group.end());
        }

        std::sort(group.begin(), group.end());
        std::sort(parents.begin(), parents.end());

        GroupType _group {parent++};
        std::set_difference(group.begin(), group.end(), parents.begin(), parents.end(), std::back_inserter(_group));

#ifdef DEBUG
        std::cout << "[";
        for (auto &index : _group) {
            std::cout << index << ", ";
        }
        std::cout << "]" << std::endl;
#endif

        PolysType merged;

        MergeGroup(_group, merged);

        std::map<Point3d, vtkIdType> newIds;

        for (auto &poly : merged) {
            auto newCell = vtkSmartPointer<vtkIdList>::New();

            for (auto &p : poly) {
                double in[] = {p.x, p.y},
                    out[3];

                BackTransform(in, out, base);

                vtkIdType id = p.id;

                if (id == NOTSET) {
                    auto itr = newIds.find(p);

                    if (itr == newIds.end()) {
                        id = pdPts->InsertNextPoint(out);
                        newIds.emplace(p, id);
                    } else {
                        id = itr->second;
                    }
                }

                newCell->InsertNextId(id);

            }

            pd->InsertNextCell(VTK_POLYGON, newCell);
            origCellIds->InsertNextValue(origId);
        }

    }

}

void Merger::MergeGroup (const GroupType &group, PolysType &merged) {
    if (group.size() == 1) {
        merged.push_back(polys.at(group.back()));

        return;
    }

    auto pts = vtkSmartPointer<vtkPoints>::New();

    IndexedPolysType indexedPolys;

    ReferencedPointsType refPts;

    SourcesType sources;
    std::size_t src = 0;

    for (auto &index : group) {
        const Poly &poly = polys.at(index);

        IndexedPoly ids;

        for (auto &p : poly) {
            vtkIdType id = pts->InsertNextPoint(p.x, p.y, p.z);

            ids.push_back(id);
            sources.emplace(id, src);

            refPts.emplace(id, p);
        }

        indexedPolys.push_back(std::move(ids));
        src++;
    }

    auto kdTree = vtkSmartPointer<vtkKdTree>::New();
    kdTree->OmitZPartitioning();
    kdTree->BuildLocatorFromPoints(pts);

    auto linesA = vtkSmartPointer<vtkPolyData>::New();
    linesA->SetPoints(pts);
    linesA->Allocate(1);

    IndexedPoly::const_iterator itrA, itrB;

    for (const auto &ids : indexedPolys) {
        for (itrA = ids.begin(); itrA != ids.end(); ++itrA) {
            itrB = itrA+1;
            if (itrB == ids.end()) {
                itrB = ids.begin();
            }

            vtkIdList *line = vtkIdList::New();
            line->InsertNextId(*itrA);
            line->InsertNextId(*itrB);

            linesA->InsertNextCell(VTK_LINE, line);

            line->Delete();
        }
    }

#ifdef DEBUG
    WriteVTK("linesA.vtk", linesA);
#endif

    auto bspTreeA = vtkSmartPointer<vtkModifiedBSPTree>::New();
    bspTreeA->SetDataSet(linesA);

    int n = 0;

    PolyConnsType polyConns;

    FindConns(linesA, kdTree, bspTreeA, polyConns, indexedPolys, sources, n);

    PolyConnsType connected {{0, {}}};
    _IdsType restricted; // keine der conns darf im gleichen punkt beginnen

    auto linesB = vtkSmartPointer<vtkPolyData>::New();
    linesB->SetPoints(pts);
    linesB->Allocate(1);

    auto bspTreeB = vtkSmartPointer<vtkModifiedBSPTree>::New();
    bspTreeB->SetDataSet(linesB);

    ConnsType firstConns;

    std::size_t i, numPolys = indexedPolys.size();

    double ptA[3], ptB[3];

    while (connected.size() < numPolys) {

        bool foundOne = false;

        for (i = 1; i < numPolys; i++) {
            if (connected.count(i) == 0) {
                const ConnsType &conns = polyConns[i];

                for (auto &conn : conns) {
                    if (connected.count(sources.at(conn.j)) == 1
                        && restricted.count(conn.j) == 0) {

                        pts->GetPoint(conn.i, ptA);
                        pts->GetPoint(conn.j, ptB);

                        if (bspTreeB->IntersectWithLine(ptA, ptB, 1e-5, nullptr, nullptr) == 0) {
                            connected[sources.at(conn.i)].push_back(conn);

                            // das andere poly auch aktualisieren
                            connected[sources.at(conn.j)].emplace_back(conn.d, conn.j, conn.i);

                            restricted.insert(conn.i);
                            restricted.insert(conn.j);

                            vtkIdList *line = vtkIdList::New();
                            line->InsertNextId(conn.i);
                            line->InsertNextId(conn.j);

                            linesB->InsertNextCell(VTK_LINE, line);

                            line->Delete();

                            bspTreeB->Modified();

                            foundOne = true;

                            firstConns.push_back(conn);

                            break;
                        }

                    }
                }
            }
        }

        if (!foundOne) {
            if (!FindConns(linesA, kdTree, bspTreeA, polyConns, indexedPolys, sources, n)) {
                throw std::runtime_error("Merging failed.");
            }
        }
    }

    std::map<std::size_t, std::vector<std::size_t>> chains;

    PolyConnsType::const_iterator itrC;

    for (itrC = connected.begin(); itrC != connected.end(); ++itrC) {
        auto &chain = chains[itrC->first];
        chain.push_back(itrC->first);

        while (chain.back() != 0) {
            chain.push_back(sources.at(connected.at(chain.back()).front().j));
        }
    }

#ifdef DEBUG
    std::cout << connected;

    decltype(chains)::const_iterator itrD;

    for (itrD = chains.begin(); itrD != chains.end(); ++itrD) {
        std::cout << itrD->first << ": [";
        for (auto &id : itrD->second) {
            std::cout << id << ", ";
        }
        std::cout << "]" << std::endl;
    }
#endif

    std::set<std::size_t> solved {0};

    std::deque<std::size_t> searchInds;

    for (i = 1; i < numPolys; i++) {
        if (connected.at(i).size() == 1) {
            searchInds.push_back(i);
        }
    }

    while (!searchInds.empty()) {
        PriosType prios;

        for (auto ind : searchInds) {
            PolyPriosType polyPrios;

#ifdef DEBUG
            std::cout << "ind " << ind << std::endl;
#endif

            const Conn &first = connected.at(ind).back();

            for (auto &conn : polyConns.at(ind)) {
                auto &src = sources.at(conn.j);

                if (polyPrios.count(src) == 1) {
                    continue;
                }

                if (conn.i != first.i
                    && conn.j != first.j
                    && restricted.count(conn.i) == 0
                    && restricted.count(conn.j) == 0) {

                    pts->GetPoint(conn.i, ptA);
                    pts->GetPoint(conn.j, ptB);

                    if (bspTreeB->IntersectWithLine(ptA, ptB, 1e-5, nullptr, nullptr) == 0) {
                        auto &chainA = chains.at(ind),
                            &chainB = chains.at(src);

                        std::set<std::size_t> _chainA(chainA.begin(), chainA.end()),
                            _chainB(chainB.begin(), chainB.end());

                        // gemeinsame eltern
                        std::set<std::size_t> shared;

                        std::set_intersection(_chainA.begin(), _chainA.end(), _chainB.begin(), _chainB.end(), std::inserter(shared, shared.end()));

                        // gemeinsame eltern müssen sich alle in solved befinden
                        if (std::includes(solved.begin(), solved.end(), shared.begin(), shared.end())) {
                            std::set<std::size_t> solvable;

                            std::set_difference(_chainA.begin(), _chainA.end(), solved.begin(), solved.end(), std::inserter(solvable, solvable.end()));
                            std::set_difference(_chainB.begin(), _chainB.end(), solved.begin(), solved.end(), std::inserter(solvable, solvable.end()));

                            polyPrios.emplace(std::piecewise_construct,
                                std::forward_as_tuple(src),
                                std::forward_as_tuple(conn, solvable, -conn.d));
                        }
                    }
                }
            }

            PolyPriosType::const_iterator itr;
            for (itr = polyPrios.begin(); itr != polyPrios.end(); ++itr) {
                prios.insert(itr->second);
            }
        }

        if (!prios.empty()) {
            auto &prio = *prios.rbegin();

#ifdef DEBUG
            std::cout << "found " << prio << std::endl;
#endif

            auto &conns = connected.at(sources.at(prio.conn.i));

            conns.push_back(prio.conn);

            connected.at(sources.at(prio.conn.j)).emplace_back(prio.conn.d, prio.conn.j, prio.conn.i);

            restricted.insert(prio.conn.i);
            restricted.insert(prio.conn.j);

            vtkIdList *line = vtkIdList::New();
            line->InsertNextId(prio.conn.i);
            line->InsertNextId(prio.conn.j);

            linesB->InsertNextCell(VTK_LINE, line);

            line->Delete();

            bspTreeB->Modified();

            solved.insert(prio.solvable.begin(), prio.solvable.end());

            searchInds.erase(std::find(searchInds.begin(), searchInds.end(), sources.at(prio.conn.i)));

            auto itr = std::find(searchInds.begin(), searchInds.end(), sources.at(prio.conn.j));

            if (itr != searchInds.end()) {
                searchInds.erase(itr);
            }
        } else {
            if (!FindConns(linesA, kdTree, bspTreeA, polyConns, indexedPolys, sources, n)) {
                break;
            }
        }
    }

#ifdef DEBUG
    std::cout << connected;
#endif

    // fallback

    double pt[3];

    if (!searchInds.empty()) {
        for (auto ind : searchInds) {
            std::vector<std::size_t> newChain;

            for (auto c : chains.at(ind)) {
                if (solved.find(c) != solved.end()) {
                    break;
                }

                newChain.push_back(c);
            }

            ConnsType &conns = connected.at(ind);

            decltype(newChain)::const_reverse_iterator itr;

            for (itr = newChain.rbegin(); itr != newChain.rend(); itr++) {
                // gesucht ist hier die kürzeste verbindung

#ifdef DEBUG
                std::cout << "itr " << *itr << std::endl;
#endif

                // polyConns.at(*itr) ist nach d sortiert

                std::shared_ptr<Conn> found;

                for (auto &conn : polyConns.at(*itr)) {
                    auto &src = sources.at(conn.j);

                    if (solved.find(src) != solved.end()) {
                        if (restricted.count(conn.i) == 0
                            // && restricted.count(conn.j) == 0
                            && std::find_if(conns.begin(), conns.end(), [&conn](const Conn &other) { return conn.i == other.i || conn.j == other.j; }) == conns.end()) {

                            pts->GetPoint(conn.i, ptA);
                            pts->GetPoint(conn.j, ptB);

                            auto intersPts = vtkSmartPointer<vtkPoints>::New();

                            auto c = bspTreeB->IntersectWithLine(ptA, ptB, 1e-5, intersPts, nullptr);

                            if (c == 0) {
                                found = std::make_shared<Conn>(conn);

                                break;
                            }

                            // wenn schnittpunkte existieren, dann müssen alle mit ptB übereinstimmen

                            vtkIdType i, numPts = intersPts->GetNumberOfPoints();

                            std::set<Point3d> foundPts {{ptB[0], ptB[1], ptB[2]}};

                            for (i = 0; i < numPts; i++) {
                                intersPts->GetPoint(i, pt);
                                foundPts.emplace(pt[0], pt[1], pt[2]);
                            }

                            if (foundPts.size() == 1) {
                                found = std::make_shared<Conn>(conn);

                                break;
                            }

                        }
                    }
                }

                if (found) {
#ifdef DEBUG
                    std::cout << "found " << *found << std::endl;
#endif

                    conns.push_back(*found);

                    connected.at(sources.at(found->j)).emplace_back(found->d, found->j, found->i);

                    restricted.insert(found->i);
                    restricted.insert(found->j);

                    vtkIdList *line = vtkIdList::New();
                    line->InsertNextId(found->i);
                    line->InsertNextId(found->j);

                    linesB->InsertNextCell(VTK_LINE, line);

                    line->Delete();

                    bspTreeB->Modified();

                    solved.insert(*itr);

                } else {
                    throw std::runtime_error("Merging failed.");
                }

            }
        }
    }

#ifdef DEBUG
    WriteVTK("linesB.vtk", linesB);
#endif

    ConnsType2 usedConns(firstConns.begin(), firstConns.end());

    IndexedPoly polyA {indexedPolys.front()};

    MergeStage1(indexedPolys, refPts, sources, firstConns, polyA);

    IndexedPolysType splitted {polyA};

    ConnsType2 leftConns;

    for (itrC = connected.begin(); itrC != connected.end(); ++itrC) {
        if (itrC->first == 0) {
            continue;
        }

        auto &conns = itrC->second;

        ConnsType::const_iterator itr;

        for (itr = conns.begin()+1; itr != conns.end(); ++itr) {
            Conn conn(0, itr->j, itr->i);

            if (usedConns.find(conn) == usedConns.end()) {
                if (itr->i < itr->j) {
                    leftConns.emplace(0, itr->i, itr->j);
                } else {
                    leftConns.insert(std::move(conn));
                }
            }
        }

    }

#ifdef DEBUG
    std::cout << "leftConns: [";
    for (auto &conn : leftConns) {
        std::cout << conn << ", ";
    }
    std::cout << "]" << std::endl;
#endif

    MergeStage2(leftConns, refPts, usedConns, splitted);

    PolysType newPolys;
    GetPolys(refPts, splitted, newPolys);

#ifdef DEBUG
    WritePolys("merged_stage2.vtk", newPolys);
#endif

    std::move(newPolys.begin(), newPolys.end(), std::back_inserter(merged));

}

bool Merger::FindConns (vtkPolyData *lines, vtkSmartPointer<vtkKdTree> kdTree, vtkSmartPointer<vtkModifiedBSPTree> bspTree, PolyConnsType &polyConns, const IndexedPolysType &indexedPolys, const SourcesType &sources, int &n) {

    vtkPoints *pts = lines->GetPoints();

    if (n > pts->GetNumberOfPoints()) {
        return false;
    }

    n += 10;

    auto foundPts = vtkSmartPointer<vtkIdList>::New();

    vtkIdType i, numPts;

    vtkIdType idB;

    auto lineIds = vtkSmartPointer<vtkIdList>::New();

    double ptA[3], ptB[3];

    bool good;

    vtkIdType j;
    vtkIdType _idA, _idB;

    std::map<std::size_t, std::set<Conn, ConnCmp>> _polyConns;

    auto line = vtkSmartPointer<vtkIdList>::New();

    for (const auto &ids : indexedPolys) {
        for (vtkIdType idA : ids) {
            pts->GetPoint(idA, ptA);

            kdTree->FindClosestNPoints(n, ptA, foundPts);

            numPts = foundPts->GetNumberOfIds();

            for (i = 0; i < numPts; i++) {
                idB = foundPts->GetId(i);

                auto srcA = sources.at(idA),
                    srcB = sources.at(idB);

                if (srcA == srcB) {
                    continue;
                }

                pts->GetPoint(idB, ptB);

                good = true;

                if (bspTree->IntersectWithLine(ptA, ptB, 1e-5, nullptr, lineIds) == 1) {
                    for (j = 0; j < lineIds->GetNumberOfIds(); j++) {
                        lines->GetCellPoints(lineIds->GetId(j), line);

                        _idA = line->GetId(0);
                        _idB = line->GetId(1);

                        if (_idA != idA && _idA != idB
                            && _idB != idA && _idB != idB) {

                            good = false;
                            break;
                        }
                    }
                }

                if (good) {
                    double d = vtkMath::Distance2BetweenPoints(ptA, ptB);

                    _polyConns[srcA].emplace(d, idA, idB);
                    _polyConns[srcB].emplace(d, idB, idA);
                }
            }
        }
    }

    decltype(_polyConns)::const_iterator itr;

    for (itr = _polyConns.begin(); itr != _polyConns.end(); ++itr) {
        auto &_conns = itr->second;

        ConnsType conns(_conns.begin(), _conns.end());
        std::sort(conns.begin(), conns.end());

        polyConns[itr->first].swap(conns);

    }

    return true;
}

void Merger::MergeStage1 (const IndexedPolysType &indexedPolys, [[maybe_unused]] const ReferencedPointsType &refPts, const SourcesType &sources, const ConnsType &conns, IndexedPoly &polyA) {

    for (const auto &conn : conns) {
        auto itrA = std::find(polyA.begin(), polyA.end(), conn.j);

        assert(itrA != polyA.end());

        IndexedPoly polyB(indexedPolys.at(sources.at(conn.i)));

        auto itrB = std::find(polyB.begin(), polyB.end(), conn.i);

        assert(itrB != polyB.end());

        std::rotate(polyA.begin(), itrA, polyA.end());
        std::rotate(polyB.begin(), itrB, polyB.end());

        IndexedPoly newPoly {polyA};
        newPoly.push_back(polyA.front());
        newPoly.push_back(polyB.front());

        newPoly.insert(newPoly.end(), polyB.rbegin(), polyB.rend());

        polyA.swap(newPoly);

    }

#ifdef DEBUG
    PolysType newPolys;
    GetPolys(refPts, {polyA}, newPolys);

    WritePolys("merged_stage1.vtk", newPolys);
#endif

}

void Merger::MergeStage2 (const ConnsType2 &conns, const ReferencedPointsType &refPts, const ConnsType2 &usedConns, IndexedPolysType &splitted) {
    std::set<Point3d> endPts;

    for (const Conn &conn : usedConns) {
        endPts.emplace(refPts.at(conn.i));
        endPts.emplace(refPts.at(conn.j));
    }

    IndexedPolysType::iterator itr;

    double vA[3], vB[3], w[3], ang, phi;

    const double n[] = {0, 0, 1};

    IndexedPoly::iterator itrA, itrB;

    IndexedPoly::iterator prev, next;

    for (auto &conn : conns) {
        for (itr = splitted.begin(); itr != splitted.end(); ++itr) {
            IndexedPoly poly(itr->begin(), itr->end());

            if (endPts.count(refPts.at(conn.i)) == 0)  {
                itrA = std::find(poly.begin(), poly.end(), conn.i);
            } else {
                Point3d::GetVec(refPts.at(conn.i), refPts.at(conn.j), w);

                itrA = poly.begin();
                while ((itrA = std::find(itrA, poly.end(), conn.i)) != poly.end()) {
                    next = itrA+1;
                    if (next == poly.end()) {
                        next = poly.begin();
                    }

                    if (itrA == poly.begin()) {
                        prev = poly.end()-1;
                    } else {
                        prev = itrA-1;
                    }

                    Point3d::GetVec(refPts.at(conn.i), refPts.at(*next), vA);
                    Point3d::GetVec(refPts.at(conn.i), refPts.at(*prev), vB);

                    ang = GetAngle(vA, vB, n);
                    phi = GetAngle(vA, w, n);

                    if (phi < ang) {
                        break;
                    }

                    ++itrA;
                }
            }

            if (itrA == poly.end()) {
                continue;
            }

            std::rotate(poly.begin(), itrA, poly.end());

            if (endPts.count(refPts.at(conn.j)) == 0)  {
                itrB = std::find(poly.begin(), poly.end(), conn.j);
            } else {
                Point3d::GetVec(refPts.at(conn.j), refPts.at(conn.i), w);

                itrB = poly.begin();
                while ((itrB = std::find(itrB, poly.end(), conn.j)) != poly.end()) {
                    next = itrB+1;
                    if (next == poly.end()) {
                        next = poly.begin();
                    }

                    if (itrB == poly.begin()) {
                        prev = poly.end()-1;
                    } else {
                        prev = itrB-1;
                    }

                    Point3d::GetVec(refPts.at(conn.j), refPts.at(*next), vA);
                    Point3d::GetVec(refPts.at(conn.j), refPts.at(*prev), vB);

                    ang = GetAngle(vA, vB, n);
                    phi = GetAngle(vA, w, n);

                    if (phi < ang) {
                        break;
                    }

                    ++itrB;
                }
            }

            if (itrB == poly.end()) {
                continue;
            }

            IndexedPoly newPolyA(poly.begin(), itrB+1);
            IndexedPoly newPolyB(itrB, poly.end());

            newPolyB.push_back(poly.front());

            splitted.erase(itr);

            splitted.push_back(std::move(newPolyA));
            splitted.push_back(std::move(newPolyB));

            endPts.emplace(refPts.at(conn.i));
            endPts.emplace(refPts.at(conn.j));

            break;
        }
    }

}
