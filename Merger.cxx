/*
Copyright 2012-2025 Ronald Römer

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

#include "Merger.h"

#include <vtkCellData.h>

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
