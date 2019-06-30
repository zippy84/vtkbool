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

#include <vector>
#include <iostream>
#include <memory>
#include <set>
#include <map>
#include <cassert>
#include <deque>
#include <iterator>

#include "Tools.h"
#include "VisPoly.h"
#include "RmTrivials.h"
#include "AABB.h"

void GetVisPoly (PolyType &poly, Tracker &tr, PolyType &res, SavedPtsType &savedPts, int ind) {
    if (poly.size() < 3) {
        return;
    }

    // std::cout << "?_ --" << std::endl;

    double x[] = {poly[ind].x, poly[ind].y};

    VertsType verts;

    int num = poly.size();

    for (int i = 0; i < num-1; i++) {
        int _i = (ind+i+1)%num;
        verts.push_back(Vert(x, poly[_i]));
    }

    double ref[] = {verts[0].pt[0]-x[0], verts[0].pt[1]-x[1]};
    Normalize(ref);

    for (Vert& v : verts) {
        v.phi = GetAngle(ref, v.r);

        // std::cout << v.phi*180/PI << std::endl;
    }

    for (int i = 0; i < num-2; i++) {
        verts[i].nxt = i+1;
    }

    /*for (Vert& v : verts) {
        std::cout << v.tag << " -> " << verts[v.nxt].tag << std::endl;
    }*/

    IdsType vp = {0, 1};

    int t = 0, u, v;

    std::vector<Bag> leftBags;

    std::map<int, int> repl;

    try {
        for (int _i = 0;; _i++) {
            u = verts.at(t).nxt;
            v = verts.at(u).nxt;

            if (v == NO_USE) {
                break;
            }

            // std::cout << "?_ > " << (u+1) << ", " << (v+1) << std::endl;

            double ptU[2], ptV[2];

            Cpy(ptU, verts.at(u).pt);
            Cpy(ptV, verts.at(v).pt);

            if (Ld(x, ptU, ptV)) {
                // std::cout << "?_ skipping" << std::endl;

                double lU = GetSqDis(x, verts.at(u)),
                    lV = GetSqDis(x, verts.at(v));

                if (lU > lV) {
                    verts.at(u).id = NO_USE;
                } else {
                    verts.at(v).id = NO_USE;

                    repl.emplace(v, u);
                }

                t = u;
                continue;
            }

            double cA = Cross(x, ptU, ptV),
                cB = Cross(verts.at(t).pt, ptU, ptV);

            /*std::cout << "?_ cA " << cA << std::endl;
            std::cout << "?_ cB " << cB << std::endl;*/

            if (cA < 0) {
                // std::cout << "?_ vis" << std::endl;

                if (vp.back() != u) {
                    vp.push_back(u);
                }

                vp.push_back(v);

                t = u;
            } else {
                if (cB > 0 || IsNear(verts.at(t).pt, ptV)) {

                    int w = v;
                    for (;;) {
                        Vert &a = verts.at(w),
                            &b = verts.at(a.nxt);

                        double *ptA = a.pt,
                            *ptB = b.pt;

                        // std::cout << "?_ 1> " << (w+1) << ", " << (a.nxt+1) << std::endl;

                        std::shared_ptr<D> d(Intersect(x, verts.at(u).r, ptA, ptB));

                        if (d
                            && d->t1 > E
                            && IsFrontfaced(verts.at(u).r, ptA, ptB)) {

                            // steht im zusammenhang mit dem skipping

                            int _u = u;

                            for (;;) {
                                try {
                                    _u = repl.at(_u);
                                } catch (...) {
                                    break;
                                }
                            }

                            if (d->t2 < E) {
                                verts.at(_u).nxt = w;
                                vp.push_back(w);
                                t = _u;
                                leftBags.push_back(Bag(_u, w, verts.at(_u).phi));

                                verts.at(w).id = NO_USE;

                            } else {

                                Point _p(d->s);

                                Vert _v(x, _p, ref, a.nxt);
                                verts.push_back(_v);

                                int k = verts.size()-1;

                                verts.at(_u).nxt = k;

                                vp.push_back(k);

                                t = _u;

                                leftBags.push_back(Bag(_u, k, verts.at(_u).phi));

                                tr.Track(a, b, _v, d->t2);

                            }

                            break;

                        } else {
                            w = a.nxt;
                        }

                        // tritt meistens dann ein, wenn det() zum nullptr führt oder d->s zu nahe an x ist
                        vtkbool_throw(w != NO_USE, "GetVisPoly", "w not set");

                    }

                } else if (cB < 0) {
                    // schnitt mit leftBags?

                    std::shared_ptr<Bag> bag;
                    std::shared_ptr<D> d;

                    while (leftBags.size() > 0 && !d) {
                        bag = std::make_shared<Bag>(leftBags.back());

                        if (bag->phi > verts.at(v).phi || std::abs(bag->phi-verts.at(v).phi) < E) {
                            d = Intersect2(verts.at(bag->f).pt, verts.at(bag->g).pt, ptU, ptV);

                            leftBags.pop_back();
                        } else {
                            break;
                        }
                    }

                    if (d) {
                        // std::cout << "?_ bag " << *bag << std::endl;

                        assert(std::find(vp.begin(), vp.end(), bag->f) != vp.end());

                        while (vp.size() > 0 && vp.back() != bag->f) {
                            // std::cout << "?_ popping_1 " << (vp.back()+1) << std::endl;
                            vp.pop_back();
                        }

                        vtkbool_throw(vp.size() > 1, "GetVisPoly", "too many pop's");

                        int _x = v;

                        int i = 0;

                        for (;;) {
                            Vert &a = verts.at(_x),
                                &b = verts.at(a.nxt);

                            double *ptA = a.pt,
                                *ptB = b.pt;

                            // std::cout << "?_ 2>" << (_x+1) << ", " << (a.nxt+1) << std::endl;

                            std::shared_ptr<D> _d;

                            if (IsNear(verts.at(bag->f).pt, ptV)) {
                                _d = std::make_shared<D>(verts.at(bag->f).pt);
                            } else {
                                _d = Intersect2(verts.at(bag->f).pt, d->s, ptB, ptA);
                            }

                            if (_d
                                && IsFrontfaced(verts.at(bag->f).r, ptA, ptB)
                                && (i > 0 || Cross(ptA, ptU, ptB) < 0)) {

                                if (IsNear(verts.at(bag->f).pt, _d->s)) {
                                    verts.at(bag->f).nxt = a.nxt;

                                    vp.push_back(a.nxt);

                                } else {
                                    if (_d->t2 > 1-E) {
                                        verts.at(bag->f).nxt = _x;

                                        vp.push_back(_x);

                                        leftBags.push_back(Bag(bag->f, _x, bag->phi));

                                        verts.at(_x).id = NO_USE;

                                    } else {

                                        Point _p(_d->s);

                                        Vert _v(x, _p, ref, a.nxt);
                                        verts.push_back(_v);

                                        int k = verts.size()-1;

                                        verts.at(bag->f).nxt = k;

                                        vp.push_back(k);

                                        leftBags.push_back(Bag(bag->f, k, bag->phi));

                                        tr.Track(b, a, _v, _d->t2);

                                    }

                                }

                                t = bag->f;

                                break;

                            } else {
                                _x = a.nxt;
                            }

                            vtkbool_throw(_x != NO_USE, "GetVisPoly", "_x not set");

                            i++;


                        }

                    } else {
                        while (vp.size() > 0) {
                            int a = vp.end()[-2],
                                b = vp.back();

                            // std::cout << "?_ popping_2 " << (vp.back()+1) << std::endl;

                            vp.pop_back();

                            vtkbool_throw(vp.size() > 0, "GetVisPoly", "too many pop's");

                            std::shared_ptr<D> d(Intersect(x, verts.at(v).r, verts.at(a).pt, verts.at(b).pt));

                            if (d) {
                                if (d->t2 < E) {
                                    int c = vp.end()[-2];

                                    if (Ld(x, verts.at(a).pt, verts.at(c).pt) || IsNear(verts.at(a).pt, ptV)) {
                                        vp.pop_back();
                                        t = vp.back();
                                    } else {
                                        t = a;
                                    }

                                } else {
                                    Point _p(d->s);

                                    Vert _v(x, _p, ref);
                                    verts.push_back(_v);

                                    int k = verts.size()-1;

                                    verts.at(a).nxt = k;

                                    vp.push_back(k);

                                    t = k;

                                    tr.Track(verts.at(a), verts.at(b), _v, d->t2);

                                }

                                break;

                            }

                        }

                        int p = v;

                        int w = verts.at(v).nxt;

                        if (w == NO_USE) {
                            break;
                        }

                        if (Ld(x, ptV, verts.at(w).pt)) {
                            p = w;
                            w = verts.at(w).nxt;
                        }

                        // std::cout << "?_ " << (v+1) << " -> " << (p+1) << std::endl;

                        double *ptW = verts.at(w).pt;

                        double cC = Cross(x, ptV, ptW),
                            cD = Cross(ptV, ptU, ptW);

                        /*std::cout << "?_ cC " << cC << std::endl;
                        std::cout << "?_ cD " << cD << std::endl;*/

                        if (cC < 0) {

                            if (cD < 0 || IsNear(ptU, ptW)) {
                                verts.at(vp.back()).nxt = p;

                                vp.push_back(p);

                                verts.at(t).id = NO_USE;

                            } else {
                                int _x = w;

                                for (;;) {
                                    Vert &a = verts.at(_x),
                                        &b = verts.at(a.nxt);

                                    double *ptA = a.pt,
                                        *ptB = b.pt;

                                    // std::cout << "?_ 3> " << (_x+1) << ", " << (a.nxt+1) << std::endl;

                                    std::shared_ptr<D> d(Intersect(x, verts.at(v).r, ptA, ptB));

                                    if (d
                                        && (!IsFrontfaced(verts.at(v).r, ptA, ptB)
                                            || IsNear(ptA, ptV))) { // spezialfall (special:1, ind:1)

                                        if (d->t2 < E) {
                                            verts.at(vp.back()).nxt = _x;
                                            vp.push_back(_x);

                                        } else {
                                            Point _p(d->s);

                                            Vert _v(x, _p, ref, a.nxt);
                                            verts.push_back(_v);
                                            int k = verts.size()-1;

                                            verts.at(vp.back()).nxt = k;

                                            vp.push_back(k);

                                            tr.Track(a, b, _v, d->t2);

                                        }

                                        break;
                                    } else {
                                        _x = a.nxt;
                                    }

                                    vtkbool_throw(_x != NO_USE, "GetVisPoly", "_x not set");

                                }
                            }

                        } else {
                            verts.at(vp.back()).nxt = p;

                            vp.push_back(p);
                        }

                    }
                }
            }

            // -

            PolyType r;
            r.push_back(poly[ind]);

            for (int _v : vp) {
                r.push_back(verts[_v]);
            }

            // std::cout << "?F " << _i << " " << GetAbsolutePath(r) << std::endl;

        }

    } catch (const std::out_of_range&) {
        vtkbool_throw(false, "GetVisPoly", "cannot access element in verts");
    }

    res.push_back(poly[ind]);

    for (int _v : vp) {
        res.push_back(verts[_v]);
    }

}

/*
list([ list(map(float, p.split(','))) for p in 'm 26.402829,29.895027 -2.132521,24.374833 -3.073759,35.133226 22.541594,1.972134 76.814397,6.720388 1.86346,-21.299507 3.34282,-38.208551 -31.800976,-2.782225 -0.800142,-0.07 -7.246298,-0.633968 -2.155836,24.641314 -6.148254,-0.537902 -8.586643,-0.751234 -4.925747,-0.430947 1.112312,-12.713787 0.198,-2.263145 0.192176,-2.196583 0.326117,-3.727536 0.327231,-3.740264 z'[2:-2].split(' ') ])
*/

void Simplify (const PolyType &poly, SavedPtsPtr &savedPts, SpecTagsPtr &specTags, PolyType &res, int skip, bool rev) {

    // der dritte anlauf um es in den griff zu bekommen

    // std::cout << "POLY " << GetAbsolutePath(poly) << std::endl;

    const double tol = 1e-3,
        sol = tol*tol;

    PolyType poly2{{poly.front()}};

    PolyType::const_iterator itr, itr2;

    for (itr = poly.begin()+1; itr != poly.end(); ++itr) {
        double sq = GetSqDis(*itr, *poly2.rbegin());
        if (sq > sol) {
            poly2.push_back(*itr);
        } else if (itr->tag == skip) {
            poly2.pop_back();
            poly2.push_back(*itr);
        }
    }

    if (GetSqDis(poly2.front(), poly2.back()) < sol) {
        poly2.pop_back();
    }

    // std::cout << "POLY2 " << GetAbsolutePath(poly2) << std::endl;

    auto GetMaxDis = [&poly2, &tol](int i, int k, int &id) {
        double last = 0, _t;

        for (int j = i+1; j < k; j++) {
            double t,
                d = GetDis(poly2[i], poly2[k], poly2[j], t);

            double delt = std::abs(last-d);

            // if (i == 77 && k == 82) {
            //     std::cout << "=> " << j << ", " << d << ", " << t;
            // }

            if (id > 0 && delt < tol) {
               if (t < _t) {
                    last = d;
                    id = j;
                    _t = t;

                    // if (i == 77 && k == 82) {
                    //     std::cout << " -> I";
                    // }
                }
            } else if (d > last) {
                last = d;
                id = j;
                _t = t;

                // if (i == 77 && k == 82) {
                //     std::cout << " -> II";
                // }
            }

            // if (i == 77 && k == 82) {
            //     std::cout << std::endl;
            // }

        }
        return last;
    };

    std::vector<PolyType> sects;

    std::deque<Pair> test{{0, static_cast<int>(poly2.size()-1)}}, pairs;

    while (!test.empty()) {
        Pair _p = test.front();
        test.pop_front();

        if (_p.g-_p.f == 1) {
            pairs.push_back(_p); // es könnten sich interne punkt darin befinden
            continue;
        }

        int id = 0;

        double d = GetMaxDis(_p.f, _p.g, id);

        // std::cout << "_" << _p << " -> " << id << ", " << d << std::endl;

        if (d > tol) {
            test.push_back({_p.f, id});
            test.push_back({id, _p.g});

        } else {
            // std::cout << _p << ", d=" << d << std::endl;

            pairs.push_back(_p);
        }
    }

    // hier bräuchte man den abschnitt aus poly zw. a.tag und b.tag
    // um die am anfang entfernten mitzunehmen

    std::set<int> tags{skip};

    for (auto &w : pairs) {
        tags.insert(poly2[w.f].tag);
        tags.insert(poly2[w.g].tag);
    }

    SpecTagsPtr _specTags(new SpecTagsType);

    std::set<Point> pts;

    for (auto &p : poly2) {
        if (tags.count(p.tag) == 1) {
            pts.insert(p);
        }
    }

    for (auto &p : poly2) {
        if (tags.count(p.tag) == 0 && pts.count(p) == 1) {
            _specTags->insert(p.tag);
        }
    }

    PolyType _res;

    std::copy_if(poly2.begin(), poly2.end(), std::back_inserter(_res), [&tags, &_specTags](const Point &p) {
        return tags.count(p.tag) == 1 || _specTags->count(p.tag) == 1;
    });

    vtkbool_throw(_res.size() > 2, "Simplify", "_res has too few points");

    _res.insert(_res.end(), _res.begin(), _res.begin()+2);

    double s;

    for (itr = _res.begin(); itr != _res.end()-2; ++itr) {
        if (_specTags->count((itr+1)->tag) == 0
            && (itr+1)->tag != skip
            && pts.find(*itr) != pts.find(*(itr+2))) {

            double d = GetDis(*itr, *(itr+2), *(itr+1), s);

            if (d < tol) {
                tags.erase((itr+1)->tag);

                if (_specTags->count(itr->tag) == 1) {
                    _specTags->erase(itr->tag);
                    tags.insert(itr->tag);

                } else if (_specTags->count((itr+2)->tag) == 1) {
                    _specTags->erase((itr+2)->tag);
                    tags.insert((itr+2)->tag);

                }

            }
        }
    }

    if (specTags) {
        specTags.swap(_specTags);
    }

    std::copy_if(poly2.begin(), poly2.end(), std::back_inserter(res), [&tags, &specTags](const Point &p) {
        return tags.count(p.tag) == 1 || (specTags && specTags->count(p.tag) == 1);
    });

    // dieser test ist zu einfach

    double areaA = GetArea(poly),
        areaB = GetArea(res);

    vtkbool_throw(areaB/areaA > .95, "Simplify", "failed");

    if (savedPts) {

        for (itr = res.begin(); itr != res.end(); ++itr) {
            itr2 = itr+1;

            if (itr2 == res.end()) {
                itr2 = res.begin();
            }

            PolyType sect(poly); // kopiert
            GetSect(itr->tag, itr2->tag, sect);

            sects.push_back(std::move(sect));
        }

        for (auto &sect : sects) {
            if (sect.size() > 2) {

                const Point &a = *sect.begin(),
                    &b = *sect.rbegin();

                double n[] = {b.x-a.x, b.y-a.y},
                    l = Normalize(n);

                double d = a.x*n[0]+a.y*n[1];

                VertsType4 verts;

                for (itr = sect.begin()+1; itr != sect.end()-1; ++itr) {
                    const Point &_p = *itr;

                    double t = _p.x*n[0]+_p.y*n[1]-d;

                    //assert(t/l > 0 && t/l < 1);

                    Vert4 v(_p, t/l);

                    v.pt[0] = a.x+t*n[0];
                    v.pt[1] = a.y+t*n[1];

                    if (rev) {
                        v.t = 1-v.t;
                    }

                    verts.push_back(std::move(v));
                }

                if (rev) {
                    std::reverse(verts.begin(), verts.end());
                    (*savedPts)[{b.tag, a.tag}] = verts;
                } else {
                    (*savedPts)[{a.tag, b.tag}] = verts;
                }
            }
        }

        // std::cout << "POLY3 " << GetAbsolutePath(res) << std::endl;

    }

}

void Align (PolyType &poly, const Point &p) {
    VertsType2 verts;

    PolyType::iterator itr;

    for (itr = poly.begin(); itr != poly.end(); ++itr) {
        if (itr->id != p.id) {
            double v[] = {itr->x-p.x, itr->y-p.y};
            verts.push_back({static_cast<int>(itr-poly.begin()), Normalize(v)});
        }
    }

    std::sort(verts.begin(), verts.end(), [](const Vert2 &a, const Vert2 &b) {
        return b.l < a.l;
    });

    VertsType2::iterator itr2, itr3;

    for (itr2 = verts.begin(); itr2 != verts.end()-1; ++itr2) {
        for (itr3 = itr2+1; itr3 != verts.end(); ++itr3) {
            double t,
                pro[2],
                d = GetDis(p, poly[itr2->i], poly[itr3->i], t, pro);

            if (d < 1e-3) {
                Cpy(poly[itr3->i].pt, pro);
            }
        }

    }
}

void Restore (const PolyType &poly, const Tracker &tr, const SavedPtsType &savedPts, PolyType &res) {
    std::map<int, Point> rr;

    PolyType::const_iterator itr, itr2;
    for (itr = poly.begin(); itr != poly.end(); ++itr) {
        res.push_back(*itr);

        itr2 = itr+1;

        if (itr2 == poly.end()) {
            itr2 = poly.begin();
        }

        if (itr == poly.begin() || itr2 == poly.begin()) {
            continue;
        }

        const Pos &posA = tr.locs.at(itr->tag),
            &posB = tr.locs.at(itr2->tag);

        if (savedPts.find({ posA.edA, posA.edB }) != savedPts.end()
            && (posA == posB
                || posA.edB == itr2->tag)) {

            double t = posA == posB ? posB.t : 1;

            auto &pts = savedPts.at({ posA.edA, posA.edB });

            for (auto &p : pts) {
                if (p.t > posA.t+E && p.t < t-E) {
                    res.push_back(p);
                } else if (IsNear(p.pt, itr->pt)) {
                    res.back() = p;
                } else if (IsNear(p.pt, itr2->pt)) {
                    rr.emplace(itr2->tag, p);
                }
            }

        }

    }

    for (auto &p : res) {
        /*if (rr.find(p.tag) != rr.end()) {
            p = rr.at(p.tag);
        }*/

        try {
            p = rr.at(p.tag);
        } catch (...) {}
    }

}

void Restore2 (const PolyType &poly, PolyType &res) {
    PolyType::const_iterator itr, itr2;

    PolyType::const_reverse_iterator itr3, itr4;

    std::map<int, Point> found;

    for (itr = res.begin(); itr != res.end(); ++itr) {
        itr2 = itr+1;

        if (itr2 == res.end()) {
            itr2 = res.begin();
        }

        if (IsNear(itr->pt, itr2->pt)) {
            itr3 = std::find_if(poly.rbegin(), poly.rend(), [&itr](const Point &p) {
                return p.tag == itr->tag;
            });

            vtkbool_throw(itr3 != poly.rend(), "Restore2", "itr3 is on end");

            itr4 = std::find_if(poly.rbegin(), poly.rend(), [&itr2](const Point &p) {
                return p.tag == itr2->tag;
            });

            vtkbool_throw(itr4 != poly.rend(), "Restore2", "itr4 is on end");

            std::map<Point, int> counts;

            do {
                counts[*itr3]++;

                if (++itr3 == poly.rend()) {
                    itr3 = poly.rbegin();
                }
            } while (itr3->tag != itr2->tag);

            counts[*itr3]++;

            PolyType single;

            for (auto &c : counts) {
                if (c.second == 1) {
                    single.push_back(c.first);
                }
            }

            vtkbool_throw(single.size() == 1, "Restore2", "too many points at the same place");

            found.emplace(itr2-res.begin(), single.front());

        }
    }

    std::map<int, Point>::const_reverse_iterator itr5;

    for (itr5 = found.rbegin(); itr5 != found.rend(); ++itr5) {
        res.insert(res.begin()+itr5->first, itr5->second);
    }

}

void SimpleRestore (const PolyType &poly, const SavedPtsType &savedPts, PolyType &res) {
    PolyType::const_iterator itr, itr2;
    for (itr = poly.begin(); itr != poly.end(); ++itr) {
        res.push_back(*itr);

        itr2 = itr+1;

        if (itr2 == poly.end()) {
            itr2 = poly.begin();
        }

        try {
            const VertsType4 &pts = savedPts.at({ itr->tag, itr2->tag });
            std::copy(pts.begin(), pts.end(), std::back_inserter(res));
        } catch (...) {}

    }
}

void GetVisPoly_wrapper (PolyType &poly, PolyType &res, int ind) {
    int i = 0;
    for (auto& p : poly) {
        p.id = i++;
    }

    vtkbool_throw(TestCW(poly), "GetVisPoly_wrapper", "poly not clockwise");

    /*std::cout << "?" << std::endl
        << "?X " << ind << std::endl
        << "?A " << GetAbsolutePath(poly) << std::endl;*/

    PolyType poly2, poly3, poly4;

    Point x(poly[ind]);

    SavedPtsPtr savedPts(new SavedPtsType);

    SpecTagsPtr specTags(new SpecTagsType);

    Simplify(poly, savedPts, specTags, poly2, x.tag, true);

    // Align(poly2, x);

    Tracker tr(poly2);

    TrivialRm(poly2, tr, ind, x).GetSimplified(poly3);

    // std::cout << "?D " << GetAbsolutePath(poly3) << std::endl;

    try {
        GetVisPoly(poly3, tr, poly4, *savedPts);

        for (auto &l : tr.locs) {
            vtkbool_throw(l.second.t < 1, "GetVisPoly_wrapper", "t not less than 1");
        }

        /*i = 0;
        for (auto &p : poly4) {
            std::cout << i++ << ". " << p << " => " << tr.locs[p.tag] << std::endl;
        }*/

        if (specTags->count(poly4[1].tag) == 1) {
            poly4.erase(poly4.begin()+1);

            if (specTags->count(poly4[1].tag) == 1) {
                poly4.erase(poly4.begin()+1);
            }
        }

        // std::copy(poly4.begin(), poly4.end(), std::back_inserter(res));

        Restore(poly4, tr, *savedPts, res);

        Restore2(poly, res);

        {
            PolyType::const_iterator itr = std::find_if(poly2.begin(), poly2.end(), [&x](const Point &p) {
                return p.tag == x.tag;
            }), itr2 = itr;

            do {
                ++itr2;

                if (itr2 == poly2.end()) {
                    itr2 = poly2.begin();
                }

                if (specTags->count(itr2->tag) == 0) {
                    if ((res.end()-1)->tag != itr2->tag) {
                        (res.end()-1)->id = NO_USE;
                        res.push_back(*itr2);
                    }

                    break;
                }

            } while (itr2 != itr);
        }

        {
            PolyType::const_reverse_iterator itr = std::find_if(poly2.rbegin(), poly2.rend(), [&x](const Point &p) {
                return p.tag == x.tag;
            }), itr2 = itr;

            do {
                ++itr2;

                if (itr2 == poly2.rend()) {
                    itr2 = poly2.rbegin();
                }

                if (specTags->count(itr2->tag) == 0) {
                    if ((res.begin()+1)->tag != itr2->tag) {
                        (res.begin()+1)->id = NO_USE;
                        res.insert(res.begin()+1, *itr2);
                    }

                    break;
                }

            } while (itr2 != itr);
        }

        /*i = 0;
        for (auto &p : res) {
            std::cout << i++ << ". " << p << std::endl;
        }*/

        /*std::cout << "{";
        for (auto &p : res) {
            std::cout << p.id << ",";
        }
        std::cout << "}" << std::endl;*/

        // trivialer test ob sich die linien innerhalb von res überschneiden

        AABB tree;

        std::vector<std::shared_ptr<Line>> lines;

        PolyType::const_iterator itr, itr2;

        std::set<Point> pts;

        for (itr = res.begin(); itr != res.end(); ++itr) {
            itr2 = itr+1;

            if (itr2 == res.end()) {
                itr2 = res.begin();
            }

            lines.emplace_back(new Line({*itr}, {*itr2}));

            pts.insert(*itr);

        }

        for (auto &line : lines) {
            tree.InsertObj(line);
        }

        Bnds bnds(-E, E, -E, E);

        for (auto &lineA : lines) {
            auto found = tree.Search(lineA);
            const Line &lA = *lineA;

            for (auto &lineB : found) {
                const Line &lB = dynamic_cast<Line&>(*lineB);

                // die linien dürfen sich nicht an den enden berühren

                if (pts.find(lA.pA) != pts.find(lB.pA)
                    && pts.find(lA.pA) != pts.find(lB.pB)
                    && pts.find(lA.pB) != pts.find(lB.pA)
                    && pts.find(lA.pB) != pts.find(lB.pB)

                    && Intersect2(lA.pA.pt, lA.pB.pt, lB.pA.pt, lB.pB.pt, bnds)) {

                    vtkbool_throw(false, "GetVisPoly_wrapper", "res is invalid - intersecting edges");
                }
            }
        }

        // alles gut

        // std::cout << "?E " << GetAbsolutePath(res) << std::endl;

    } catch (...) {
        throw;
    }

    /*
    PolyType res2;
    std::copy_if(res.begin(), res.end(), std::back_inserter(res2), [](const Point &p) {
        return p.id != NO_USE;
    });

    res.swap(res2);
    */

    // std::cout << GetAbsolutePath(res) << std::endl;

}
