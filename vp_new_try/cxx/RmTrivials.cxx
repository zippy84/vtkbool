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

#include "RmTrivials.h"

void AlignPts (PolyType &poly, int ind) {
    // verändert poly

    int num = poly.size();

    double x[] = {poly[ind].x, poly[ind].y};

    double ei[] = {1, 0};

    VertsType2 verts;

    for (int i = 0; i < num; i++) {
        if (i != ind
            && !IsNear(poly[i].pt, x)) {

            double r[] = {poly[i].pt[0]-x[0], poly[i].pt[1]-x[1]};

            double d = Normalize(r),
                phi = GetAngle(ei, r);

            verts.push_back(Vert2(i, r, d, phi));
        }
    }

    std::map<std::string, IdsType> phis;

    for (int i = 0; i < verts.size(); i++) {
        std::ostringstream out;
        out << std::fixed << std::setprecision(4) << verts[i].phi;

        phis[out.str()].push_back(i);
    }

    for (auto& itr : phis) {
        IdsType &ids = itr.second;
        if (ids.size() > 1) {
            std::sort(ids.begin(), ids.end(), [&verts](int &a, int &b) {
                return verts[b].d > verts[a].d;
            });

            for (int id : ids) {
                Vert2 &v = verts[id];

                double phi_ = std::stod(itr.first)-v.phi;

                // dreht den ortsvektor um den kleinen winkel

                double w[] = {v.d*v.r[0], v.d*v.r[1]},
                    r[] = {w[0]-w[1]*phi_, w[0]*phi_+w[1]};

                poly[v.i].pt[0] = x[0]+r[0];
                poly[v.i].pt[1] = x[1]+r[1];

            }

        }
    }

}

void RemoveInterals (PolyType &poly, int skip) {
    // modifiert poly

    int num = poly.size();

    IdsType internal;

    for (int i = 0; i < num; i++) {
        int j = (i+1)%num,
            k = (i+2)%num;

        Point &pA = poly[i],
            &pB = poly[j],
            &pC = poly[k];

        if (Ld(pA.pt, pB.pt, pC.pt) < 1e-2) {
            internal.push_back(j);
        }
    }

    internal.erase(std::remove(internal.begin(), internal.end(), skip), internal.end());

    std::cout << "internal=" << internal.size() << std::endl;

    /*
    IdsType ids(num);
    std::iota(ids.begin(), ids.end(), 0);

    IdsType::const_iterator itr;
    for (itr = internal.begin(); itr != internal.end(); ++itr) {
        ids.erase(std::find(ids.begin(), ids.end(), *itr));
    }
    */

    // entfernen
    assert(std::is_sorted(internal.begin(), internal.end()));

    IdsType::reverse_iterator itr2;
    for (itr2 = internal.rbegin(); itr2 != internal.rend(); ++itr2) {
        std::cout << "del " << *itr2 << std::endl;

        poly.erase(poly.begin()+(*itr2));
    }

}

void AddPocket (Pair2 &pair, VertsType3 &verts, double *rot, double d) {
    int num = verts.size();

    IdsType &pocket = pair.pocket;

    int k = pair.i;

    pocket.push_back(k);

    for (;;) {
        k = (k+1)%num;
        pocket.push_back(k);
        if (k == pair.j) {
            break;
        }
    }

    if (pocket.size() > 2) {
        IdsType::iterator itr;
        for (itr = pocket.begin()+1; itr != pocket.end()-1; ++itr) {
            double *pt = verts[*itr].pt,
                e = (pt[0]*rot[0]+pt[1]*rot[1])-d;

            if (std::abs(e) > E) {
                pair.side = static_cast<Side>(e > 0);
                break;
            }
        }

    }

}

void RemovePockets (VertsType3 &verts, VertsType3 &good, double *rot, double d) {

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

    auto AddPair = [&](int i, int j) {
        pairs.push_back(Pair2(i, j));
        return pairs.size()-1;
    };

    if (pairs.size() > 0) {
        std::vector<Pair2> newPairs;

        std::vector<Pair2>::iterator itr2, itr3;
        for (itr2 = pairs.begin(); itr2 != pairs.end(); ++itr2) {
            if (itr2->dir == Dir::FORWARD) {
                AddPocket(*itr2, verts, rot, d);
            }
        }

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

        for (itr4 = grps.begin(); itr4 != grps.end(); ++itr4) {
            if (itr4->dir == Dir::UNDEFINED) {
                if (itr4 != grps.begin() && (itr4-1)->dir == Dir::BACKWARD) {
                    itr4->dir = Dir::BACKWARD;
                } else {
                    AddPocket(pairs[itr4->ids[0]], verts, rot, d);
                    itr4->dir = Dir::FORWARD;
                }
            }
        }

        if (grps.size() > 0 && grps[0].dir == Dir::FORWARD) {
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

                        AddPocket(pairs[_ids.back()], verts, rot, d);

                        _c = false;

                        break;
                    }

                }

                if (_c && next+1 != grps.end()) {
                    IdsType &_ids2 = (next+1)->ids;

                    std::cout << "(2) new pair (" << pairs[_ids.back()].j
                              << ", " << pairs[_ids2.front()].j << ")" << std::endl;

                    _ids.push_back(AddPair(pairs[_ids.back()].j, pairs[_ids2.front()].j));
                    AddPocket(pairs[_ids.back()], verts, rot, d);

                    _ids2.erase(_ids2.begin());
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

void RemoveTrivials (PolyType &poly, PolyType &res, int ind) {
    // poly ist eine kopie

    AlignPts(poly, ind);

    int numPts = poly.size();

    int iA = (ind+1)%numPts,
        iB = (ind+numPts-1)%numPts;

    double rA[] = {poly[iA].x-poly[ind].x, poly[iA].y-poly[ind].y},
        rB[] = {poly[iB].x-poly[ind].x, poly[iB].y-poly[ind].y};

    VertsType3 verts;

    int num = numPts;

    for (int i = 0; i < numPts; i++) {
        poly[i].id = i;

        verts.push_back(poly[i]);
        Vert3 &last = verts.back();

        if (i == ind) {
            continue;
        }

        int j = (i+1)%numPts;

        std::map<double, Vert3> newVerts;

        double t;

        if (Ld(poly[ind].pt, poly[iA].pt, poly[i].pt) < 1e-2) {
            t = GetT(poly[ind].pt, poly[iA].pt, poly[i].pt);

            if (t > E) {
                last.src = Src::A;
                last.t = t;
            }

        } else {
            std::shared_ptr<D> d(Intersect(poly[ind].pt, rA, poly[i].pt, poly[j].pt));

            if (d && d->t1 > E) {
                newVerts.emplace(d->t2, Vert3(d->s, Src::A, d->t1));
                num++;
            }
        }

        if (Ld(poly[ind].pt, poly[iB].pt, poly[i].pt) < 1e-2) {
            t = GetT(poly[ind].pt, poly[iB].pt, poly[i].pt);

            if (t > E) {
                last.src = Src::B;
                last.t = t;
            }

        } else {
            std::shared_ptr<D> d(Intersect(poly[ind].pt, rB, poly[i].pt, poly[j].pt));

            if (d && d->t1 > E) {
                newVerts.emplace(d->t2, Vert3(d->s, Src::B, d->t1));
                num++;
            }
        }

        for (auto& itr : newVerts) {
            verts.push_back(itr.second);
        }

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
            return p.id == iA;
        });

        std::rotate(good.begin(), first, good.end());

        std::cout << "[";
        for (const auto& g : good) {
            std::cout << "(" << g.i << ", " << g.t << "), ";
        }
        std::cout << "]" << std::endl;

        double rot[] = {-rA[1], rA[0]};
        Normalize(rot);

        double d = rot[0]*poly[ind].pt[0]+rot[1]*poly[ind].pt[1];

        RemovePockets(verts, good, rot, d);

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
            return p.id == iB;
        });

        std::rotate(good.begin(), first, good.end());

        std::cout << "[";
        for (const auto& g : good) {
            std::cout << "(" << g.i << ", " << g.t << "), ";
        }
        std::cout << "]" << std::endl;

        double rot[] = {rB[1], -rB[0]};
        Normalize(rot);

        double d = rot[0]*poly[ind].pt[0]+rot[1]*poly[ind].pt[1];

        RemovePockets(verts, good, rot, d);

    }

    // abschluss

    //std::reverse(verts.begin(), verts.end());

    verts.erase(std::remove_if(verts.begin(), verts.end(), [](const Vert3 &p) {
        return p.rm;
    }), verts.end());

    // ind ist dann der erste im neuen verts

    std::rotate(verts.begin(), std::find_if(verts.begin(), verts.end(), [&ind](const Vert3 &p) {
        return ind == p.id;
    }), verts.end());

    std::transform(verts.begin(), verts.end(), std::back_inserter(res), [](Vert3 &p) {
        return p;
    });

    RemoveInterals(res, 0);

    assert(res.front().id == ind);

}
