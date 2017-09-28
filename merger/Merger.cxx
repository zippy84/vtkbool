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

        std::cout << ">> [";
        for (int id : ids) {
            std::cout << id << ", ";
        }
        std::cout << "]" << std::endl;

        // erstmal nur die, die holes haben

        if (!ids.empty()) {
            PolysType group = {polys[i]};
            for (int id : ids) {
                group.push_back(polys[id]);
            }

            PolyType merged;
            Merge(group, merged);

            std::cout << GetAbsolutePath(merged) << std::endl;

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

    std::vector<Pair> edges;

    for (auto& poly : group) {
        for (Point& p : poly) {
            // alte id sichern

            if (p.id != NO_USE) {
                oldPtIds[a] = p.id;
            }

            p.id = a++;

            src.push_back(b);
        }

        int num = poly.size();

        for (int i = 0; i < num; i++) {
            int j = (i+1)%num;

            edges.push_back({poly[i].id, poly[j].id});

            pts->InsertNextPoint(poly[i].x, poly[i].y, 0);
        }

        b++;
    }

    vtkKdTree *tree = vtkKdTree::New();
    tree->OmitZPartitioning();
    tree->BuildLocatorFromPoints(pts);

    for (Pair& edge : edges) {
        std::cout << edge << std::endl;
    }

    int numPts = pts->GetNumberOfPoints();

    IdsType _ids(numPts);
    std::iota(_ids.begin(), _ids.end(), 0);

    int curr = 1;

    ResType res;

    std::vector<Pair> cons;

    while (!_ids.empty()) {

        for (int _id  : _ids) {

            double _pt[3];
            pts->GetPoint(_id, _pt);

            vtkIdList *found = vtkIdList::New();
            tree->FindClosestNPoints(curr*2, _pt, found);

            // sucht nach gültigen verbindungen

            IdsType valids;

            for (int i = (curr-1)*2; i < found->GetNumberOfIds(); i++) {
                valids.push_back(found->GetId(i));
            }

            std::cout << "_id: " << _id << ", curr: " << curr << ", valids: [";
            for (int v : valids) {
                std::cout << v << ", ";
            }
            std::cout << "]" << std::endl;

            valids.erase(std::remove_if(valids.begin(), valids.end(), [&](int id) {
                if (src[id] == src[_id]) {
                    return true;
                } else {
                    for (Pair& edge : edges) {
                        // edge mit id oder _id auslassen
                        if (edge.f != id && edge.g != id
                            && edge.f != _id && edge.g != _id) {

                            double eA[3], eB[3];

                            pts->GetPoint(edge.f, eA);
                            pts->GetPoint(edge.g, eB);

                            double pt[3];
                            pts->GetPoint(id, pt);

                            std::shared_ptr<D> d = Intersect2(_pt, pt, eA, eB);
                            if (d) {
                                // es existiert ein schnitt zw. der verbindung und einem der edges
                                return true;
                            }
                        }
                    }

                    return false;
                }
            }), valids.end());

            std::cout << "=> valids: [";
            for (int v : valids) {
                std::cout << v << ", ";
            }
            std::cout << "]" << std::endl;

            for (int id : valids) {
                double pt[3];
                pts->GetPoint(id, pt);

                double v[3];
                vtkMath::Subtract(_pt, pt, v);

                double d = vtkMath::Norm(v);

                int pA = src[_id],
                    pB = src[id];

                res[pA].emplace(Pair(_id, id), d);
                res[pB].emplace(Pair(id, _id), d);

            }

            found->Delete();
        }

        _ids.clear();

        // erstmal die keys betrachten
        for (const auto& r : res) {
            std::cout << r.first << std::endl;
        }

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

                std::cout << "not viewed: [";
                int i = 0;

                for (auto& poly : group) {
                    if (viewed.count(i) == 0) {
                        std::cout << i << ", ";

                        for (Point& p : poly) {
                            _ids.push_back(p.id);
                        }
                    }

                    i++;
                }
                std::cout << "]" << std::endl;

                curr++;

                break;
            }
        }

        std::cout << "_ids: [";
        for (int id : _ids) {
            std::cout << id << ", ";
        }
        std::cout << "]" << std::endl;

    }

    // dies sind die verbindungen
    for (auto& con : cons) {
        std::cout << con << std::endl;
    }

    // fügt die polygone zusammen

    int num = numPolys;

    std::map<int, int> repls;

    for (int i = 0; i < numPolys; i++) {
        repls[i] = i;
    }

    for (auto& con : cons) {
        int a = con.f,
            b = con.g;

        int pA = repls[src[a]],
            pB = repls[src[b]];

        std::cout << pA << ", " << pB << std::endl;

        PolyType &polyA = group[pA],
            &polyB = group[pB];

        std::deque<Point> deqA(polyA.begin(), polyA.end()),
            deqB(polyB.begin(), polyB.end());

        auto itrA = std::find_if(deqA.begin(), deqA.end(), [&a](const Point& pt) { return pt.id == a; });
        auto itrB = std::find_if(deqB.begin(), deqB.end(), [&b](const Point& pt) { return pt.id == b; });

        std::rotate(deqA.begin(), itrA, deqA.end());
        std::rotate(deqB.begin(), itrB, deqB.end());

        std::cout << deqA[0].id << " -> " << a << std::endl;
        std::cout << deqB[0].id << " -> " << b << std::endl;

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

        group.push_back(std::move(newPoly));

        repls[src[a]] = num;
        repls[src[b]] = num;

        num++;

    }

    for (Point &p : group[num-1]) {
        if (oldPtIds.count(p.id) == 0) {
            merged.push_back({p.pt});
        } else {
            merged.push_back({p.pt, oldPtIds[p.id]});
        }

    }

    tree->Delete();
    pts->Delete();

}
