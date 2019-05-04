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

#include "Merger.h"

#include <set>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <memory>
#include <deque>
#include <string>
#include <sstream>

#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkKdTree.h>
#include <vtkMath.h>

#include "Tools.h"
#include "AABB.h"

void Merger::AddPoly (PolyType &poly) {
    polys.push_back(poly);
}

void Merger::GetMerged (PolysType &res) {

    int numPolys = polys.size();

    std::vector<IdsType> all(numPolys);

    for (int i = 0; i < numPolys; i++) {
        IdsType outer;

        for (int j = 0; j < numPolys; j++) {
            if (i != j) {
                if (TestPIP(polys[j], polys[i][0])) {
                    outer.push_back(j);
                }
            }
        }

        for (int j : outer) {
            all[j].push_back(i);
        }

    }

    //std::vector<IdsType> all = {{2}, {}, {1, 3}, {}, {0}};

    int i = 0;

    for (auto& ids : all) {
        std::set<int> excludes;

        for (int id : ids) {
            excludes.insert(all[id].begin(), all[id].end());
        }

        ids.erase(std::remove_if(ids.begin(), ids.end(), [&excludes](int id) {
            return excludes.count(id) > 0;
        }), ids.end());

        // std::cout << ">> [";
        // for (int id : ids) {
        //     std::cout << id << ", ";
        // }
        // std::cout << "]" << std::endl;

        // erstmal nur die, die holes haben

        if (!ids.empty()) {
            PolysType group = {polys[i]};
            for (int id : ids) {
                group.push_back(polys[id]);
            }

            PolyType merged;
            Merge(group, merged);

            res.push_back(merged);

        } else {
            res.push_back(polys[i]);
        }

        i++;

    }

}

void Merger::Merge (PolysType &group, PolyType &merged) {
    int numPolys = group.size();

    int a = 0,
        b = 0;

    IdsType src;

    std::map<int, int> oldPtIds;

    vtkPoints *pts = vtkPoints::New();

    AABB treeA;

    PolyType::iterator itr, itr2;

    for (auto& poly : group) {
        for (itr = poly.begin(); itr != poly.end(); ++itr) {
            // alte id sichern
            if (itr->id != NO_USE) {
                oldPtIds[a] = itr->id;
            }

            itr->id = a++;

            src.push_back(b);

            pts->InsertNextPoint(itr->x, itr->y, 0);
        }

        for (itr = poly.begin(); itr != poly.end(); ++itr) {
            itr2 = itr+1;

            if (itr2 == poly.end()) {
                itr2 = poly.begin();
            }

            std::shared_ptr<Line> line(new Line({*itr}, {*itr2}));

            treeA.InsertObj(line);
        }

        b++;
    }

    vtkKdTree *treeB = vtkKdTree::New();
    treeB->OmitZPartitioning();
    treeB->BuildLocatorFromPoints(pts);

    int numPts = pts->GetNumberOfPoints();

    const Bnds bnds(-E, E, -E, E);

    IdsType ids(numPts);
    std::iota(ids.begin(), ids.end(), 0);

    int curr = 1;

    ResType res;

    std::vector<Pair> cons;

    while (!ids.empty()) {
        for (int idA : ids) {
            double ptA[3],
                ptB[3];

            pts->GetPoint(idA, ptA);

            vtkIdList *cls = vtkIdList::New();
            treeB->FindClosestNPoints(curr*2, ptA, cls);

            int numCls = cls->GetNumberOfIds();

            // sucht nach gültigen verbindungen

            std::vector<int> _ids;

            int srcA = src[idA];

            for (int i = (curr-1)*2; i < numCls; i++) {
                int idB = cls->GetId(i),
                    srcB = src[idB];

                if (srcA != srcB) {
                    _ids.push_back(idB);
                }
            }

            _ids.erase(std::remove_if(_ids.begin(), _ids.end(), [&](int idB) {
                pts->GetPoint(idB, ptB);

                std::shared_ptr<Line> lineA(new Line({ptA}, {ptB}));

                auto found = treeA.Search(lineA);

                for (auto &f : found) {
                    Line &lineB = dynamic_cast<Line&>(*f);

                    if (lineB.pA.id != idA && lineB.pB.id != idA
                        && lineB.pA.id != idB && lineB.pB.id != idB
                        && Intersect2(ptA, ptB, lineB.pA.pt, lineB.pB.pt, bnds)) {

                        return true;
                    }
                }

                return false;

            }), _ids.end());

            for (int idB : _ids) {
                pts->GetPoint(idB, ptB);

                double v[3];
                vtkMath::Subtract(ptA, ptB, v);

                double d = vtkMath::Norm(v);

                int pA = src[idA],
                    pB = src[idB];

                res[pA].emplace(Pair(idA, idB), d);
                res[pB].emplace(Pair(idB, idA), d);
            }

            cls->Delete();
        }

        ids.clear();

        // erstmal die keys betrachten

        // for (const auto& r : res) {
        //     std::cout << r.first << std::endl;
        // }

        cons.clear();

        // sucht an bereits besuchten polygonen nach verbindungen zu anderen polygonen
        // es wird das polygon in viewed aufgenommen, dass die kürzeste distanz hat

        // stellt sicher, dass alle polygone zusammenhängen

        std::set<int> viewed = {0};

        while (viewed.size() < numPolys) {

            std::shared_ptr<G> g;

            for (int v : viewed) {
                for (auto& r : res[v]) {
                    if (viewed.count(src[r.first.g]) == 0) {
                        // das poly wurde noch nicht besucht

                        if (!g || r.second < g->d) {
                            g = std::make_shared<G>(r.second, r.first);
                        }
                    }
                }

            }

            if (g) {
                viewed.insert(src[g->con.g]);
                cons.push_back(g->con);

            } else {

                // std::cout << "not viewed: [";

                int i = 0;

                for (auto& poly : group) {
                    if (viewed.count(i) == 0) {
                        // std::cout << i << ", ";

                        for (Point& p : poly) {
                            ids.push_back(p.id);
                        }
                    }

                    i++;
                }

                // std::cout << "]" << std::endl;

                curr++;

                break;
            }
        }

        // std::cout << "ids: [";
        // for (int id : ids) {
        //     std::cout << id << ", ";
        // }
        // std::cout << "]" << std::endl;

    }

    // dies sind die verbindungen
    // for (auto& con : cons) {
    //     std::cout << con << std::endl;
    // }

    typedef std::deque<Point> PolyTypeD;

    auto FindId = [&](PolyTypeD &poly, int id, int end, int s) -> PolyTypeD::iterator {
        PolyTypeD::iterator itr;

        double pt[3];
        pts->GetPoint(end, pt);

        int num = poly.size();

        for (itr = poly.begin(); itr != poly.end(); ++itr) {
            if (itr->id == id) {
                int i = itr-poly.begin();

                double v[] = {pt[0]-itr->x, pt[1]-itr->y};
                Normalize(v);

                int iA = (i+1)%num,
                    iB = (i+num-1)%num;

                Point &pA = poly[iA],
                    &pB = poly[iB];

                double wA[] = {pA.x-itr->x, pA.y-itr->y},
                    wB[] = {pB.x-itr->x, pB.y-itr->y};

                Normalize(wA);
                Normalize(wB);

                double ref = GetAngle(wA, wB),
                    phi = GetAngle(wA, v);

                if (phi > ref || (s == 0 && phi < ref)) {
                    break;
                }
            }
        }

        assert(itr != poly.end());

        return itr;
    };


    // fügt die polygone zusammen

    int num = numPolys;

    std::map<int, int> repls;

    for (int i = 0; i < numPolys; i++) {
        repls[i] = i;
    }

    for (auto& con : cons) {
        int a = con.f,
            b = con.g;

        int pA = repls.at(src[a]),
            pB = repls.at(src[b]);

        // std::cout << a << "->" << pA << ", " << b << "->" << pB << std::endl;

        PolyType &polyA = group.at(pA),
            &polyB = group.at(pB);

        // std::cout << "polyA " << GetAbsolutePath(polyA) << std::endl;
        // std::cout << "polyB " << GetAbsolutePath(polyB) << std::endl;

        PolyTypeD deqA(polyA.begin(), polyA.end()),
            deqB(polyB.begin(), polyB.end());

        /*
        auto itrA = std::find_if(deqA.begin(), deqA.end(), [&a](const Point& pt) { return pt.id == a; });
        auto itrB = std::find_if(deqB.begin(), deqB.end(), [&b](const Point& pt) { return pt.id == b; });
        */

        auto itrA = FindId(deqA, a, b, pA);
        auto itrB = FindId(deqB, b, a, pB);

        std::rotate(deqA.begin(), itrA, deqA.end());
        std::rotate(deqB.begin(), itrB, deqB.end());

        // std::cout << deqA[0].id << " -> " << a << std::endl;
        // std::cout << deqB[0].id << " -> " << b << std::endl;

        PolyType newPoly;

        if (pA == 0) {
            newPoly.push_back(deqA.front());
            newPoly.insert(newPoly.end(), deqA.rbegin(), deqA.rend());
        } else {
            newPoly.insert(newPoly.end(), deqA.begin(), deqA.end());
            newPoly.push_back(deqA.front());
        }

        if (pB == 0) {
            newPoly.push_back(deqB.front());
            newPoly.insert(newPoly.end(), deqB.rbegin(), deqB.rend());
        } else {
            newPoly.insert(newPoly.end(), deqB.begin(), deqB.end());
            newPoly.push_back(deqB.front());
        }

        // std::cout << GetAbsolutePath(newPoly) << std::endl;

        group.push_back(newPoly);

        // repls[src[a]] = num;
        // repls[src[b]] = num;

        for (auto &v : newPoly) {
            repls[src[v.id]] = num;
        }

        num++;

    }

    for (Point &p : group[num-1]) {
        if (oldPtIds.count(p.id) == 0) {
            merged.push_back({p.pt});
        } else {
            merged.push_back({p.pt, oldPtIds[p.id]});
        }

    }

    std::reverse(merged.begin(), merged.end());

    assert(!TestCW(merged));

    treeB->Delete();
    pts->Delete();

}
