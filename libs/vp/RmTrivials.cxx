/*
Copyright 2012-2018 Ronald Römer

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

void MarkInternals (VertsType4 &verts, int skip) {
    int num = verts.size();

    std::set<int> ids;

    for (int i = 0; i < num; i++) {
        int j = (i+1)%num,
            k = (i+2)%num;

        Vert4 &pA = verts[i],
            &pB = verts[j],
            &pC = verts[k];

        if (j != skip
            && !IsNear(pA.pt, pC.pt)
            && Ld(pA.pt, pB.pt, pC.pt) < 1e-3) {

            pB.marked = std::make_shared<Marked>(i, k);

            ids.insert(j);
        }
    }

    std::set<int> invalids;

    for (int id : ids) {
        for (int i = 0; i < num; i++) {
            if (ids.count(i) == 0
                && IsNear(verts[id].pt, verts[i].pt)) {
                invalids.insert(id);
            }
        }
    }

    for (int id : invalids) {
        ids.erase(id);
        verts[id].marked.reset();
    }

    std::map<Pair, IdsType> col;

    for (int id : ids) {
        auto &m = verts[id].marked;

        while (ids.count(m->a) == 1) {
            m->a = (verts[m->a].marked)->a;
        }

        while (ids.count(m->b) == 1) {
            m->b = (verts[m->b].marked)->b;
        }

        Vert4 &pA = verts[m->a],
            &pB = verts[m->b];

        double v[] = {pB.x-pA.x, pB.y-pA.y},
            w[] = {verts[id].x-pA.x, verts[id].y-pA.y};

        double v_ = Normalize(v),
            w_ = Normalize(w);

        m->t = w_/v_;

        col[{m->a, m->b}].push_back(id);

    }

    if (col.size() > 1) {
        std::map<Pair, IdsType>::iterator itr;

        for (itr = col.begin(); itr != col.end(); ++itr) {
            const Pair &edge = itr->first;
            IdsType nodes(itr->second);

            if (nodes.size() > 1) {

                std::sort(nodes.begin(), nodes.end(), [&verts](const int &a, const int &b) {
                    return (verts[a].marked)->t < (verts[b].marked)->t;
                });

                nodes.insert(nodes.begin(), edge.f);
                nodes.push_back(edge.g);

                IdsType keep = {0, static_cast<int>(nodes.size())-1};

                int e;

                do {
                    e = 0;

                    int l = keep.size();
                    for (int i = 0; i < l-1; i++) {
                        int j = i+1;

                        int a = keep[i],
                            b = keep[j];

                        if (b-a > 1) {
                            int c = a+(b-a)/2;

                            Vert4 &pA = verts[nodes[a]],
                                &pB = verts[nodes[b]],
                                &pC = verts[nodes[c]];

                            if (!(Ld(pA.pt, pB.pt, pC.pt) < 1e-3)) {
                                e++;

                                keep.insert(keep.begin()+j, c);

                                break;
                            }

                        }
                    }

                } while (e > 0);

                if (keep.size() > 2) {
                    IdsType::iterator itr2;

                    for (itr2 = keep.begin()+1; itr2 != keep.end()-1; ++itr2) {
                        std::cout << "Rm " << *itr2 << std::endl;
                        verts[nodes[*itr2]].marked.reset();
                    }

                }

            }
        }
    }

    for (auto& vert : verts) {
        if (vert.marked) {
            std::cout << *(vert.marked) << std::endl;
        }
    }
}

void RemoveInternals (VertsType4 &verts) {
    verts.erase(std::remove_if(verts.begin(), verts.end(), [](const Vert4 &vert) {
        return vert.marked;
    }), verts.end());

}

void AlignPts (VertsType4 &verts, int ind) {
    int num = verts.size();

    double x[] = {verts[ind].x, verts[ind].y};

    VertsType2 verts2;

    for (int i = 0; i < num; i++) {
        if (i != ind
            && !verts[i].marked
            && !IsNear(verts[i].pt, x)) {

            double r[] = {verts[i].pt[0]-x[0], verts[i].pt[1]-x[1]},
                d = Normalize(r);

            verts2.push_back(Vert2(i, r, d));
        }
    }

    std::sort(verts2.begin(), verts2.end());

    /*
    for (auto& v : verts2) {
        std::cout << v << std::endl;
    }
    */

    for (auto& v : verts2) {
        double _v[] = {-v.r[1], v.r[0]},
            _d = _v[0]*x[0]+_v[1]*x[1];

        for (auto& w : verts2) {
            if (&v == &w) {
                continue;
            }

            double d = -(verts[w.i].x*_v[0]+verts[w.i].y*_v[1]-_d);

            double _l = Ld(x, verts[v.i].pt, verts[w.i].pt);

            //std::cout << verts[w.i].id << " -> " << _l << std::endl;

            if (_l < 1e-3) {
                verts[w.i].pt[0] += d*_v[0];
                verts[w.i].pt[1] += d*_v[1];
            }
        }
    }

    // korrektur der internen

    for (auto& vert : verts) {
        auto &m = vert.marked;

        if (m) {
            const Vert4 &pA = verts[m->a],
                &pB = verts[m->b];

            double v[] = {pB.x-pA.x, pB.y-pA.y};

            vert.pt[0] = pA.x+m->t*v[0];
            vert.pt[1] = pA.y+m->t*v[1];
        }
    }

}

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

    for (auto& p : pairs) {
        std::cout << "pair " << p << std::endl;
    }

    std::vector<Pair2>::iterator itr2, itr3;

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

                std::cout << "erasing from " << (itr2-pairs.begin()) << std::endl;

                pairs.erase(itr2, pairs.end());
                break;

            }

        }

        for (itr2 = pairs.begin(); itr2 != pairs.end(); ++itr2) {
            if (itr2->dir == Dir::UNDEFINED) {
                IdsType pocket;
                GetPocket(*itr2, pocket);

                std::set<int> ss;

                IdsType::iterator itr4;

                for (itr4 = pocket.begin()+1; itr4 != pocket.end()-1; ++itr4) {
                    double e = (verts[*itr4].x*rot[0]+verts[*itr4].y*rot[1])-d;

                    if (std::abs(e) > E) {
                        ss.insert(int(e/std::abs(e)));
                    }
                }

                assert(ss.size() == 1);

                if (HasArea(pocket)) {

                    PolyType _poly;

                    for (auto& id : pocket) {
                        _poly.push_back(verts[id].pt);
                    }

                    if (src == Src::B) {
                        std::reverse(_poly.begin(), _poly.end());
                    }

                    if (*(ss.begin()) == 1 ^ TestCW(_poly)) {
                        itr2->dir = Dir::BACKWARD;
                    } else {
                        itr2->dir = Dir::FORWARD;
                    }

                } else {
                    if (*(ss.begin()) == 1) {
                        itr2->dir = Dir::BACKWARD;
                    } else {
                        itr2->dir = Dir::FORWARD;
                        itr2->side = Side::IN;
                    }
                }

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

            for (itr3 = itr2; itr3->dir == itr2->dir && itr3 != pairs.end(); ++itr3) {
                ids.push_back(i++);
            }

            grps.push_back(Grp(itr2->dir, ids));

            itr2 = itr3;
        }

        std::vector<Grp>::iterator itr4, itr5;

        {
            // Ausgabe

            for (itr4 = grps.begin(); itr4 != grps.end(); ++itr4) {
                std::cout << "(" << static_cast<int>(itr4->dir) << ", [";
                for (auto id : itr4->ids) {
                    std::cout << id << ", ";
                }
                std::cout << "]), ";
            }

            std::cout << std::endl;
        }

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

                    for (itr5 = itr4; itr5->dir == itr4->dir && itr5 != grps.end(); ++itr5) {
                        IdsType &ids = itr5->ids;
                        ids2.insert(ids2.end(), ids.begin(), ids.end());
                    }

                    grps2.push_back(Grp(itr4->dir, ids2));

                    itr4 = itr5;
                }

                grps2.swap(grps);

                {
                    // Ausgabe

                    std::cout << ">> ";

                    for (itr4 = grps.begin(); itr4 != grps.end(); ++itr4) {
                        std::cout << "(" << static_cast<int>(itr4->dir) << ", [";
                        for (auto id : itr4->ids) {
                            std::cout << id << ", ";
                        }
                        std::cout << "]), ";
                    }

                    std::cout << std::endl;
                }

                auto next = std::find_if(grps.begin(), grps.end(), [](const Grp &g) {
                    return g.dir == Dir::BACKWARD;
                });

                if (next == grps.end()) {
                    break;
                }

                bool _c = true;

                IdsType &ids = next->ids;

                Pair2 &last = pairs[ids.back()];

                IdsType &_ids = (next-1)->ids;

                IdsType::reverse_iterator itr6;
                for (itr6 = _ids.rbegin(); itr6 != _ids.rend(); ++itr6) {
                    Pair2 &p = pairs[*itr6];

                    int j = itr6-_ids.rbegin();

                    double tA = verts[p.i].t,
                        tB = verts[p.j].t,
                        tC = verts[last.j].t;

                    bool q = tC-tA > E && tB-tC > E && p.side == Side::OUT;

                    std::cout << j <<  ", (" << tA << ", " << tB << ", " << tC << ") " << q << std::endl;

                    if (q) {
                        std::cout << "(1) new pair (" << p.i << ", " << last.j << ")" << std::endl;

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

                        std::cout << "(2) new pair (" << pairs[_ids.back()].j
                                  << ", " << pairs[_ids2.front()].j << ")" << std::endl;

                        _ids.push_back(AddPair(pairs[_ids.back()].j, pairs[_ids2.front()].j));
                        AssignSide(pairs[_ids.back()], src);

                        _ids2.erase(_ids2.begin());

                    } else {
                        // siehe bspw. (complex:3, ind:0)

                        for (itr6 = _ids.rbegin(); itr6 != _ids.rend(); ++itr6) {
                            Pair2 &p = pairs[*itr6];

                            double tA = verts[p.j].t,
                                tB = verts[last.j].t;

                            if (std::abs(tA-tB) < E) {
                                std::cout << "(3) new pair (" << p.i
                                    << ", " << last.j << ")" << std::endl;

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
            std::cout << *itr2 << std::endl;

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
    VertsType4 poly_(poly.begin(), poly.end());

    MarkInternals(poly_, ind);
    AlignPts(poly_, ind);

    poly.clear();
    poly.insert(poly.end(), poly_.begin(), poly_.end());

    RemoveInternals(poly_);

    int ind_ = (std::find_if(poly_.begin(), poly_.end(), [&](const Point &p) { return p.id == ind; }))-poly_.begin();

    int numPts = poly_.size();

    int iA = (ind_+1)%numPts,
        iB = (ind_+numPts-1)%numPts;

    double rA[] = {poly_[iA].x-poly_[ind_].x, poly_[iA].y-poly_[ind_].y},
        rB[] = {poly_[iB].x-poly_[ind_].x, poly_[iB].y-poly_[ind_].y};

    int num = numPts;

    std::map<int, std::vector<Src>> srcs;

    for (int i = 0; i < numPts; i++) {
        verts.push_back(poly_[i]);
        Vert3 &last = verts.back();

        if (i == ind_) {
            continue;
        }

        int j = (i+1)%numPts;

        std::multimap<double, Vert3> newVerts;

        double t;

        if (Ld(poly_[ind_].pt, poly_[iA].pt, poly_[i].pt) < 1e-3) {
            t = GetT(poly_[ind_].pt, poly_[iA].pt, poly_[i].pt);

            if (t > 1-E) {
                last.src = Src::A;
                last.t = t;

                srcs[last.id].push_back(Src::A);
            }

        } else {
            std::shared_ptr<D> d(Intersect(poly_[ind_].pt, rA, poly_[i].pt, poly_[j].pt));

            if (d && d->t1 > E) {
                newVerts.emplace(d->t2, Vert3(d->s, Src::A, d->t1));
                num++;
            }
        }

        if (Ld(poly_[ind_].pt, poly_[iB].pt, poly_[i].pt) < 1e-3) {
            t = GetT(poly_[ind_].pt, poly_[iB].pt, poly_[i].pt);

            if (t > 1-E) {
                last.src = Src::B;
                last.t = t;

                srcs[last.id].push_back(Src::B);
            }

        } else {
            std::shared_ptr<D> d(Intersect(poly_[ind_].pt, rB, poly_[i].pt, poly_[j].pt));

            if (d && d->t1 > E) {
                newVerts.emplace(d->t2, Vert3(d->s, Src::B, d->t1));
                num++;
            }
        }

        for (auto& itr : newVerts) {
            verts.push_back(itr.second);
        }

    }

    int i = 0;
    for (auto& v : verts) {
        std::cout << i++ << " " << v << std::endl;
    }

    auto nxt = std::find_if(verts.begin(), verts.end(), [&](const Vert3 &p) {
        return p.id == ind;
    });

    ind_ = nxt-verts.begin();

    assert(nxt != verts.end());

    nxt = nxt+1;

    if (nxt == verts.end()) {
        nxt = verts.begin();
    }

    if (srcs.count(nxt->id) == 1
        && srcs[nxt->id].size() == 2) {

        double _rA[] = {-rA[0], -rA[1]};

        int bnd;

        for (int i = 0; i < num; i++) {
            int j = (i+1)%num;

            std::shared_ptr<D> d(Intersect(verts[ind_].pt, _rA, verts[i].pt, verts[j].pt));

            if (d && d->t1 > E) {
                bnd = i;
                break;
            }

        }

        for (int i = 0; i < num; i++) {
            int j = (ind_+i)%num;

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

            if (j == ind_) {
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
        std::cout << "A()" << std::endl;

        for (int i = 0; i < num; i++) {
            verts[i].i = i;
        }

        VertsType3 good;

        std::copy_if(verts.begin(), verts.end(), std::back_inserter(good), [](const Vert3 &p) {
            return p.src == Src::A;
        });

        auto first = std::find_if(good.begin(), good.end(), [&](const Vert3 &p) {
            return p.id == poly_[iA].id;
        });

        std::rotate(good.begin(), first, good.end());

        std::cout << "[";
        for (const auto& g : good) {
            std::cout << "(" << g.i << ", " << g.t << "), ";
        }
        std::cout << "]" << std::endl;

        double rot[] = {-rA[1], rA[0]};
        Normalize(rot);

        double d = rot[0]*x.pt[0]+rot[1]*x.pt[1];

        RemovePockets(good, rot, d, Src::A);

    }

    {
        std::cout << "B()" << std::endl;

        std::reverse(verts.begin(), verts.end());

        for (int i = 0; i < num; i++) {
            verts[i].i = i;
        }

        VertsType3 good;

        std::copy_if(verts.begin(), verts.end(), std::back_inserter(good), [](const Vert3 &p) {
            return p.src == Src::B;
        });

        auto first = std::find_if(good.begin(), good.end(), [&](const Vert3 &p) {
            return p.id == poly_[iB].id;
        });

        std::rotate(good.begin(), first, good.end());

        std::cout << "[";
        for (const auto& g : good) {
            std::cout << "(" << g.i << ", " << g.t << "), ";
        }
        std::cout << "]" << std::endl;

        double rot[] = {rB[1], -rB[0]};
        Normalize(rot);

        double d = rot[0]*x.pt[0]+rot[1]*x.pt[1];

        RemovePockets(good, rot, d, Src::B);

    }

    // abschluss

    //std::reverse(verts.begin(), verts.end());

    verts.erase(std::remove_if(verts.begin(), verts.end(), [](const Vert3 &p) {
        return p.rm;
    }), verts.end());

    // ind ist dann der erste im neuen verts

    std::rotate(verts.begin(), std::find_if(verts.begin(), verts.end(), [&](const Vert3 &p) {
        return ind == p.id;
    }), verts.end());

    VertsType4 verts_(verts.begin(), verts.end());

    MarkInternals(verts_, 0);
    RemoveInternals(verts_);

    std::transform(verts_.begin(), verts_.end(), std::back_inserter(res), [](Vert4 &p) {
        return p;
    });

    assert(res.front().id == ind);

}

void AddInternals (PolyType &origin, PolyType &poly) {
    PolyType res;

    int numA = poly.size(),
        numB = origin.size();

    for (int i = 0; i < numA; i++) {
        int j = (i+1)%numA;

        Point &pA = poly[i],
            &pB = poly[j];

        double u[] = {pB.x-pA.x, pB.y-pA.y},
            r = Normalize(u);

        res.push_back(pA);

        if (IsNear(pA.pt, pB.pt)) {
            // spezialfall (special:0, ind:13)

            assert(pA.id != NO_USE && pB.id != NO_USE);

            int k = pB.id;
            while (k != pA.id) {
                k = (k+1)%numB;

                res.push_back(origin[k]);
            }

            res.pop_back();

            continue;
        }

        VertsType5 verts;

        for (auto& p : origin) {
            if (IsOnSeg(pA.pt, pB.pt, p.pt)) {
                double v[] = {p.x-pA.x, p.y-pA.y},
                    t = Normalize(v)/r;
                verts.push_back({p, t});
            }
        }

        std::sort(verts.begin(), verts.end());

        /*
        for (auto& v : verts) {
            std::cout << v << std::endl;
        }
        */

        VertsType5::iterator itr;

        if (verts.size() > 1) {
            std::map<int, double> all;

            all[pA.id] = 0;
            all[pB.id] = 1;

            for (auto& v : verts) {
                all[v.id] = v.t;
            }

            for (itr = verts.begin(); itr != verts.end()-1; ++itr) {
                if (((itr+1)->t-itr->t) < E) {
                    int idA = itr->id,
                        idB = (itr+1)->id;

                    int cA, cB;

                    {
                        int a = (idA+1)%numB,
                            b = (idA+numB-1)%numB;

                        //cA = all.count(a)+all.count(b);
                        cA = (all.count(a) == 1 && all[a] < all[idA])
                            +(all.count(b) == 1 && all[b] > all[idA]);
                    }

                    {
                        int a = (idB+1)%numB,
                            b = (idB+numB-1)%numB;

                        //cB = all.count(a)+all.count(b);
                        cB = (all.count(a) == 1 && all[a] < all[idB])
                            +(all.count(b) == 1 && all[b] > all[idB]);
                    }

                    std::cout << "(" << idA << "->" << cA << ", "
                        << idB << "->" << cB << ") => " << poly.front().id << std::endl;

                    assert(cA > 0 || cB > 0);

                    if (cA > cB) {
                        (itr+1)->valid = false;
                    } else {
                        itr->valid = false;
                    }

                }
            }
        }

        for (auto& v : verts) {
            if (v.valid) {
                res.push_back(v);
            }
        }

    }

    poly.swap(res);

}
