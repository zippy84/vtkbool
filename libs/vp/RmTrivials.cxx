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

#include <memory>
#include <utility>
#include <string>
#include <sstream>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include <cassert>
#include <numeric>
#include <map>
#include <set>

#include "RmTrivials.h"

void TrivialRm::RemovePockets (VertsType3 &good, double *rot, double d, Src src) {

    std::vector<Pair2> pairs;

    VertsType3::iterator itr;
    for (itr = good.begin(); itr != good.end()-1; ++itr) {
        Vert3 &pA = *itr,
            &pB = *(itr+1);

        pairs.push_back(Pair2(pA.i, pB.i,
            pB.t > pA.t ? Dir::FORWARD : Dir::BACKWARD));

        if (std::abs(pB.t-pA.t) < E) {
            pairs.back().dir = Dir::UNDEFINED;
        }
    }

    /*for (auto& p : pairs) {
        std::cout << "?_ pair " << p << std::endl;
    }*/

    std::vector<Pair2>::iterator itr2, itr3;

    auto ResolveUD = [this, &rot, &d, &src](Pair2 &p) -> void {
        IdsType pocket;
        GetPocket(p, pocket);

        IdsType::iterator itr;

        std::vector<int> _s;

        for (itr = pocket.begin()+1; itr != pocket.end()-1; ++itr) {
            double e = (verts[*itr].x*rot[0]+verts[*itr].y*rot[1])-d;

            if (std::abs(e) > E) {
                _s.push_back(int(e/std::abs(e)));
            }
        }

        // _s kann 0en und 1en enthalten
        // es kommt aber nur auf den anfang oder das ende an

        if (HasArea(pocket)) {

            PolyType _poly;

            for (auto& id : pocket) {
                _poly.push_back(verts[id].pt);
            }

            if (src == Src::B) {
                std::reverse(_poly.begin(), _poly.end());
            }

            if (*(_s.begin()) == 1 ^ TestCW(_poly)) {
                p.dir = Dir::BACKWARD;
            } else {
                p.dir = Dir::FORWARD;
            }

        } else {
            std::set<int> ss(_s.begin(), _s.end());

            if (ss.size() == 1) {
                if (*(ss.begin()) == 1) {
                    p.dir = Dir::BACKWARD;
                } else {
                    p.dir = Dir::FORWARD;
                    p.side = Side::IN;
                }
            } else {
                assert(false);
            }
        }
    };

    if (pairs.size() > 0) {
        for (itr2 = pairs.begin(); itr2 != pairs.end(); ++itr2) {
            IdsType pocket;
            GetPocket(*itr2, pocket);

            PolyType _poly;

            for (auto& id : pocket) {
                _poly.push_back(verts[id].pt);
            }

            if (TestPIP(_poly, x)) {

                if ((src == Src::A && !TestCW(_poly))
                    || (src == Src::B && TestCW(_poly))) {
                    // siehe bspw. (complex:3, ind:0)

                    break;
                }

                // std::cout << "?_ erasing from " << (itr2-pairs.begin()) << std::endl;

                pairs.erase(itr2, pairs.end());
                break;

            }

        }

        for (itr2 = pairs.begin(); itr2 != pairs.end(); ++itr2) {
            if (itr2->dir == Dir::UNDEFINED) {
                ResolveUD(*itr2);
            }

        }

    }

    auto AddPair = [&](int i, int j) {
        pairs.push_back(Pair2(i, j));
        return pairs.size()-1;
    };

    if (pairs.size() > 0) {
        std::vector<Pair2> newPairs;

        std::vector<Grp> grps;
        int i = 0;

        itr2 = pairs.begin();

        while (itr2 != pairs.end()) {

            IdsType ids;

            for (itr3 = itr2; itr3 != pairs.end() && itr3->dir == itr2->dir; ++itr3) {
                ids.push_back(i++);
            }

            grps.push_back(Grp(itr2->dir, ids));

            itr2 = itr3;
        }

        std::vector<Grp>::iterator itr4, itr5;

        /*{
            for (itr4 = grps.begin(); itr4 != grps.end(); ++itr4) {
                std::cout << "?_ (" << itr4->dir << ", [";
                for (auto id : itr4->ids) {
                    std::cout << id << ", ";
                }
                std::cout << "])" << std::endl;
            }
        }*/

        if ((grps.begin())->dir == Dir::BACKWARD) {
            if (grps.size() > 1) {
                pairs[(grps.begin()+1)->ids.front()].i = pairs[(grps.begin())->ids.front()].i;

                grps.erase(grps.begin());
            }
        }

        for (itr2 = pairs.begin(); itr2 != pairs.end(); ++itr2) {
            if (itr2->dir == Dir::FORWARD) {
                AssignSide(*itr2, src);
            }
        }

        if (grps.size() > 0) {
            for (;;) {
                std::vector<Grp> grps2;

                itr4 = grps.begin();

                while (itr4 != grps.end()) {
                    IdsType ids2;

                    for (itr5 = itr4; itr5 != grps.end() && itr5->dir == itr4->dir; ++itr5) {
                        IdsType &ids = itr5->ids;
                        ids2.insert(ids2.end(), ids.begin(), ids.end());
                    }

                    grps2.push_back(Grp(itr4->dir, ids2));

                    itr4 = itr5;
                }

                grps2.swap(grps);

                /*{
                    for (itr4 = grps.begin(); itr4 != grps.end(); ++itr4) {
                        std::cout << "?_ > (" << itr4->dir << ", [";
                        for (auto id : itr4->ids) {
                            std::cout << id << ", ";
                        }
                        std::cout << "])" << std::endl;
                    }
                }*/

                auto next = std::find_if(grps.begin(), grps.end(), [](const Grp &g) {
                    return g.dir == Dir::BACKWARD;
                });

                if (next == grps.end()) {
                    break;
                }

                bool _c = true;

                IdsType &ids = next->ids;

                Pair2 &last = pairs[ids.back()];

                vtkbool_throw(grps.size() > 1, "TrivialRm::RemovePockets", "too few grps");

                IdsType &_ids = (next-1)->ids;

                IdsType::reverse_iterator itr6;
                for (itr6 = _ids.rbegin(); itr6 != _ids.rend(); ++itr6) {
                    Pair2 &p = pairs[*itr6];

                    int j = itr6-_ids.rbegin();

                    double tA = verts[p.i].t,
                        tB = verts[p.j].t,
                        tC = verts[last.j].t;

                    bool q = tC-tA > E && tB-tC > E && p.side == Side::OUT;

                    // std::cout << "?_ " << j <<  ", (" << tA << ", " << tB << ", " << tC << ") " << q << std::endl;

                    if (q) {
                        // std::cout << "?_ (1) new pair (" << p.i << ", " << last.j << ")" << std::endl;

                        _ids.erase(_ids.end()-j-1, _ids.end());
                        _ids.push_back(AddPair(p.i, last.j));

                        AssignSide(pairs[_ids.back()], src);

                        _c = false;

                        break;
                    }

                }

                if (_c) {
                    if (next+1 != grps.end()) {
                        IdsType &_ids2 = (next+1)->ids;

                        int jA = pairs[_ids.back()].j,
                            jB = pairs[_ids2.front()].j;

                        /*std::cout << "?_ (2) new pair (" << jA
                            << ", " << jB << ")" << std::endl;*/

                        double tA = verts[jA].t,
                            tB = verts[jB].t;

                        // std::cout << "?_ [" << tA << ", " << tB << "]" << std::endl;

                        // (jA, jB) kann in abh. von t alle drei zustände von Dir annehmen

                        bool fr = tA < tB;

                        if (std::abs(tB-tA) < E) {
                            Pair2 p(jA, jB);
                            ResolveUD(p);

                            // std::cout << "?_ dir -> " << p.dir << std::endl;

                            fr = p.dir == Dir::FORWARD;
                        }

                        // assert(std::abs(tB-tA) > E);

                        if (fr) {
                            _ids.push_back(AddPair(jA, jB));
                            AssignSide(pairs[_ids.back()], src);

                            _ids2.erase(_ids2.begin());
                        } else {
                            ids.clear();
                            ids.push_back(AddPair(jA, jB));

                            pairs[ids.back()].dir = Dir::BACKWARD;

                            _ids2.erase(_ids2.begin());

                            if (_ids2.empty()) {
                                grps.erase(next+1);
                            }

                            continue;
                        }

                    } else {
                        // siehe bspw. (complex:3, ind:0)

                        for (itr6 = _ids.rbegin(); itr6 != _ids.rend(); ++itr6) {
                            Pair2 &p = pairs[*itr6];

                            double tA = verts[p.j].t,
                                tB = verts[last.j].t;

                            if (std::abs(tA-tB) < E) {
                                /*std::cout << "?_ (3) new pair (" << p.i
                                    << ", " << last.j << ")" << std::endl;*/

                                _ids.erase(itr6.base()-1, _ids.end());

                                _ids.push_back(AddPair(p.i, last.j));
                                AssignSide(pairs[_ids.back()], src);

                                break;
                            }

                        }

                    }
                }

                grps.erase(next);

            }

            assert(grps.size() == 1);

            IdsType &result = grps.front().ids;

            std::transform(result.begin(), result.end(), std::back_inserter(newPairs), [&pairs](const int id) {
                return pairs[id];
            });

        }

        for (itr2 = newPairs.begin(); itr2 != newPairs.end(); ++itr2) {
            // std::cout << *itr2 << std::endl;

            if (itr2->side == Side::OUT) {
                IdsType &pocket = itr2->pocket;

                IdsType::iterator itr3;
                for (itr3 = pocket.begin()+1; itr3 != pocket.end()-1; ++itr3) {
                    verts[*itr3].rm = true;
                }
            }
        }

    }

}

void TrivialRm::GetPocket (Pair2 &pair, IdsType &pocket) {
    int num = verts.size();
    int k = pair.i;

    pocket.push_back(k);

    while (k != pair.j) {
        k = (k+1)%num;
        pocket.push_back(k);
    }
}

void TrivialRm::AssignSide (Pair2 &pair, Src src) {
    IdsType &pocket = pair.pocket;

    GetPocket(pair, pocket);

    if (pair.side != Side::NONE) {
        return;
    }

    if (pocket.size() > 2) {
        if (src == Src::B) {
            std::reverse(pocket.begin(), pocket.end());
        }

        PolyType _poly;

        for (auto& id : pocket) {
            _poly.push_back(verts[id].pt);
        }

        pair.side = TestCW(_poly) ? Side::OUT : Side::IN;

    }

}

bool TrivialRm::HasArea (const IdsType &pocket) {
    bool area = true;

    int num = pocket.size();

    if (num%2 == 1) {
        for (int i = 0; i < (num-1)/2; i++) {
            area = !IsNear(verts[pocket[i]].pt, verts[pocket[num-i-1]].pt);
        }
    }

    return area;
}

void TrivialRm::GetSimplified (PolyType &res) {
    auto itr = std::find_if(poly.begin(), poly.end(), [&](const Point &p) { return p.id == ind; });

    int ind2 = itr-poly.begin();

    int numPts = poly.size();

    int iA = (ind2+1)%numPts,
        iB = (ind2+numPts-1)%numPts;

    double rA[] = {poly[iA].x-poly[ind2].x, poly[iA].y-poly[ind2].y},
        rB[] = {poly[iB].x-poly[ind2].x, poly[iB].y-poly[ind2].y};

    int num = numPts;

    std::map<int, std::vector<Src>> srcs;

    for (int i = 0; i < numPts; i++) {
        verts.push_back(poly[i]);
        Vert3 &last = verts.back();

        if (i == ind2) {
            continue;
        }

        int j = (i+1)%numPts;

        std::multimap<double, Vert3> newVerts;

        double t;

        if (Ld(poly[ind2].pt, poly[iA].pt, poly[i].pt)) {
            t = GetT(poly[ind2].pt, poly[iA].pt, poly[i].pt);

            if (t > 1-E) {
                last.src = Src::A;
                last.t = t;

                srcs[last.id].push_back(Src::A);
            }

        } else {
            std::shared_ptr<D> d(Intersect(poly[ind2].pt, rA, poly[i].pt, poly[j].pt));

            if (d && d->t1 > E) {
                // std::cout << *d << std::endl;

                Vert3 _v(d->s, Src::A, d->t1);
                newVerts.emplace(d->t2, _v);
                num++;

                tr.Track(poly[j], poly[i], _v, 1-d->t2);
            }
        }

        if (Ld(poly[ind2].pt, poly[iB].pt, poly[i].pt)) {
            t = GetT(poly[ind2].pt, poly[iB].pt, poly[i].pt);

            if (t > 1-E) {
                last.src = Src::B;
                last.t = t;

                srcs[last.id].push_back(Src::B);
            }

        } else {
            std::shared_ptr<D> d(Intersect(poly[ind2].pt, rB, poly[i].pt, poly[j].pt));

            if (d && d->t1 > E) {
                // std::cout << *d << std::endl;

                Vert3 _v(d->s, Src::B, d->t1);
                newVerts.emplace(d->t2, _v);
                num++;

                tr.Track(poly[j], poly[i], _v, 1-d->t2);
            }
        }

        for (auto& itr : newVerts) {
            verts.push_back(itr.second);
        }

    }

    // es folgt ein trivialer test

    VertsType3::iterator itr2, itr3;

    for (itr2 = verts.begin(); itr2 != verts.end(); ++itr2) {
        if (itr2->src != Src::NONE
            && tr.locs.at(itr2->tag).t < E) {
            itr3 = itr2 != verts.begin() ? itr2-1 : verts.end()-1;

            if (itr3->id == NO_USE
                && itr2->tag == tr.locs.at(itr3->tag).edA) {

                itr3->rm = true;
            }
        }
    }

    auto nxt = std::find_if(verts.begin(), verts.end(), [&](const Vert3 &p) {
        return p.id == ind;
    });

    ind2 = nxt-verts.begin();

    assert(nxt != verts.end());

    nxt = nxt+1;

    if (nxt == verts.end()) {
        nxt = verts.begin();
    }

    if (srcs.count(nxt->id) == 1
        && srcs[nxt->id].size() == 2) {

        double _rA[] = {-rA[0], -rA[1]};

        int bnd = NO_USE;

        for (int i = 0; i < num; i++) {
            int j = (i+1)%num;

            std::shared_ptr<D> d(Intersect(verts[ind2].pt, _rA, verts[i].pt, verts[j].pt));

            if (d && d->t1 > E) {
                bnd = i;
                break;
            }

        }

        // bnd ist nur dann nicht gesetzt, wenn verts[ind2].pt mit einem anderen punkt in verts bzw. poly übereinstimmt
        // das polygon ist malformed!
        // dieses assert nimmt assert(ss.size() == 1); in zeile 104 vorweg
        // siehe special:3, ind:0
        //assert(bnd != NO_USE);

        vtkbool_throw(bnd != NO_USE, "TrivialRm::GetSimplified", "bnd not set");

        for (int i = 0; i < num; i++) {
            int j = (ind2+i)%num;

            if (verts[j].src != Src::NONE) {
                verts[j].src = Src::A;
            }

            if (j == bnd) {
                break;
            }
        }

        for (int i = 0; i < num; i++) {
            int j = (bnd+i)%num;

            if (verts[j].src != Src::NONE
                && verts[j].src == Src::A) {
                verts[j].rm = true;
            }

            if (j == ind2) {
                break;
            }
        }

        verts.erase(std::remove_if(verts.begin(), verts.end(), [](const Vert3 &p) {
            return p.rm;
        }), verts.end());

    }

    /*for (Vert3& v : verts) {
        std::cout << v << std::endl;
    }*/

    {
        // std::cout << "?_ A()" << std::endl;

        for (int i = 0; i < num; i++) {
            verts[i].i = i;
        }

        PolyType _verts;
        std::copy(verts.begin(), verts.end(), std::back_inserter(_verts));

        // std::cout << "?B " << GetAbsolutePath(_verts) << std::endl;

        VertsType3 good;

        std::copy_if(verts.begin(), verts.end(), std::back_inserter(good), [](const Vert3 &p) {
            return p.src == Src::A;
        });

        auto first = std::find_if(good.begin(), good.end(), [&](const Vert3 &p) {
            return p.id == poly[iA].id;
        });

        std::rotate(good.begin(), first, good.end());

        /*std::cout << "[";
        for (const auto& g : good) {
            std::cout << "(" << g.i << ", " << g.t << "), ";
        }
        std::cout << "]" << std::endl;*/

        double rot[] = {-rA[1], rA[0]};
        Normalize(rot);

        double d = rot[0]*x.pt[0]+rot[1]*x.pt[1];

        RemovePockets(good, rot, d, Src::A);

        RemoveRedundants(good);

    }

    {
        // std::cout << "?_ B()" << std::endl;

        std::reverse(verts.begin(), verts.end());

        for (int i = 0; i < num; i++) {
            verts[i].i = i;
        }

        PolyType _verts;
        std::copy(verts.begin(), verts.end(), std::back_inserter(_verts));

        // std::cout << "?C " << GetAbsolutePath(_verts) << std::endl;

        VertsType3 good;

        std::copy_if(verts.begin(), verts.end(), std::back_inserter(good), [](const Vert3 &p) {
            return p.src == Src::B;
        });

        auto first = std::find_if(good.begin(), good.end(), [&](const Vert3 &p) {
            return p.id == poly[iB].id;
        });

        std::rotate(good.begin(), first, good.end());

        /*std::cout << "[";
        for (const auto& g : good) {
            std::cout << "(" << g.i << ", " << g.t << "), ";
        }
        std::cout << "]" << std::endl;*/

        double rot[] = {rB[1], -rB[0]};
        Normalize(rot);

        double d = rot[0]*x.pt[0]+rot[1]*x.pt[1];

        RemovePockets(good, rot, d, Src::B);

        RemoveRedundants(good);

    }

    // abschluss

    verts.erase(std::remove_if(verts.begin(), verts.end(), [](const Vert3 &p) {
        return p.rm;
    }), verts.end());

    // ind ist dann der erste im neuen verts

    std::rotate(verts.begin(), std::find_if(verts.begin(), verts.end(), [&](const Vert3 &p) {
        return ind == p.id;
    }), verts.end());

    std::copy(verts.begin(), verts.end(), std::back_inserter(res));

    assert(res.front().id == ind);

}

void TrivialRm::RemoveRedundants (const VertsType3 &good) {

    std::set<int> tags{x.tag};

    for (const Vert3 &v : good) {
        tags.insert(v.tag);
    }

    int num = verts.size();

    VertsType3::const_iterator itr;

    for (itr = good.begin(); itr != good.end()-1; ++itr) {
        const Vert3 &a = verts[(itr->i+1)%num],
            &b = verts[(itr->i+num-1)%num];

        if ((tags.count(a.tag) == 1 || a.rm)
            && (tags.count(b.tag) == 1 || b.rm)) {

            verts[itr->i].rm = true;

        }
    }

}
