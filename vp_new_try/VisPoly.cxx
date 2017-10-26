#include <vector>
#include <iostream>
#include <memory>

#include "Tools.h"
#include "VisPoly.h"
#include "RmTrivials.h"

void GetVisPoly (PolyType &poly, PolyType &res, int ind) {
    double x[] = {poly[ind].x, poly[ind].y};

    VertsType verts;

    int num = poly.size();

    for (int i = 0; i < num-1; i++) {
        int _i = (ind+i+1)%num;
        double pt[] = {poly[_i].x, poly[_i].y};
        verts.push_back(Vert(x, pt, _i));
    }

    double ref[] = {verts[0].pt[0]-x[0], verts[0].pt[1]-x[1]};
    Normalize(ref);

    for (Vert& v : verts) {
        v.phi = GetAngle(ref, v.r);

        std::cout << v.phi*180/PI << std::endl;
    }

    for (int i = 0; i < num-2; i++) {
        verts[i].nxt = i+1;
    }

    for (Vert& v : verts) {
        std::cout << v.id << ", " << v.nxt << std::endl;
    }

    IdsType vp = {0, 1};

    int t = 0, u, v;

    std::vector<Bag> leftBags;

    for (;;) {
        u = verts[t].nxt;
        v = verts[u].nxt;

        if (v == NO_USE) {
            break;
        }

        std::cout << "> " << u << ", " << v << std::endl;

        double *ptU = verts[u].pt,
            *ptV = verts[v].pt;

        std::cout << "orig " << verts[u].id << ", " << verts[v].id << std::endl;

        if (Ld(x, ptU, ptV) < 1e-2) {
            std::cout << "skipping" << std::endl;

            t = u;
            continue;
        }

        double cA = Cross(x, ptU, ptV),
            cB = Cross(verts[t].pt, ptU, ptV);

        std::cout << "cA " << cA << std::endl;
        std::cout << "cB " << cB << std::endl;

        if (cA < 0) {
            std::cout << "vis" << std::endl;

            if (vp.back() != u) {
                vp.push_back(u);
            }

            vp.push_back(v);

            t = u;
        } else {
            if (cB > 0) {

                int w = v;
                for (;;) {
                    Vert &a = verts[w],
                        &b = verts[a.nxt];

                    double *ptA = a.pt,
                        *ptB = b.pt;

                    std::cout << ">> " << a.id << ", " << b.id << std::endl;

                    std::shared_ptr<D> d(Intersect(x, verts[u].r, ptA, ptB));

                    if (d
                        && IsFrontfaced(verts[u].r, ptA, ptB)) {

                        std::cout << "edge (" << a.id << ", " << b.id << ")" << std::endl;

                        if (d->t2 < E) {
                            verts[u].nxt = w;
                            vp.push_back(w);
                            t = u;
                            leftBags.push_back(Bag(u, w, verts[u].phi));

                        } else {

                            Vert _v(x, d->s, ref, a.nxt);
                            verts.push_back(_v);

                            int k = verts.size()-1;

                            verts[u].nxt = k;

                            vp.push_back(k);

                            t = u;

                            leftBags.push_back(Bag(u, k, verts[u].phi));

                        }

                        break;

                    } else {
                        w = a.nxt;
                    }

                }

            } else if (cB < 0) {
                // schnitt mit leftBags?
                std::shared_ptr<Bag> bag;
                while (leftBags.size() > 0 && (leftBags.back().phi-verts[v].phi) > -E) {
                    bag = std::make_shared<Bag>(leftBags.back());
                    leftBags.pop_back();
                }

                std::shared_ptr<D> d;

                if (bag) {
                    std::cout << "bag " << *bag << std::endl;

                    d = Intersect2(verts[bag->f].pt, verts[bag->g].pt, ptU, ptV);
                }

                if (d) {

                    while (vp.size() > 0 && vp.back() != bag->f) {
                        std::cout << "popping_1 " << vp.back() << std::endl;
                        vp.pop_back();
                    }

                    int _x = v;

                    int i = 0;

                    for (;;) {
                        Vert &a = verts[_x],
                            &b = verts[a.nxt];

                        double *ptA = a.pt,
                            *ptB = b.pt;

                        std::cout << ">>" << a.id << ", " << b.id << std::endl;

                        std::shared_ptr<D> _d;

                        if (IsNear(verts[bag->f].pt, ptV)) {
                            _d = std::make_shared<D>(verts[bag->f].pt);
                        } else {
                            _d = Intersect2(verts[bag->f].pt, d->s, ptB, ptA);
                        }

                        if (_d
                            && IsFrontfaced(verts[bag->f].r, ptA, ptB)
                            && (i > 0 || Cross(ptA, ptU, ptB) < 0)) {

                            if (IsNear(verts[bag->f].pt, _d->s)) {
                                verts[bag->f].nxt = a.nxt;

                                vp.push_back(a.nxt);

                            } else {
                                if (_d->t2 > 1-E) {
                                    verts[bag->f].nxt = _x;

                                    vp.push_back(_x);

                                    leftBags.push_back(Bag(bag->f, _x, bag->phi));

                                } else {

                                    Vert _v(x, _d->s, ref, a.nxt);
                                    verts.push_back(_v);

                                    int k = verts.size()-1;

                                    verts[bag->f].nxt = k;

                                    vp.push_back(k);

                                    leftBags.push_back(Bag(bag->f, k, bag->phi));

                                }

                            }

                            t = bag->f;

                            break;

                        } else {
                            _x = a.nxt;
                        }

                        i++;


                    }

                } else {
                    while (vp.size() > 0) {
                        int a = vp.end()[-2],
                            b = vp.back();

                        std::cout << "popping_2 " << vp.back() << std::endl;

                        vp.pop_back();

                        std::shared_ptr<D> d(Intersect(x, verts[v].r, verts[a].pt, verts[b].pt));

                        if (d) {
                            if (d->t2 < E) {
                                int c = vp.end()[-2];

                                if (Ld(x, verts[a].pt, verts[c].pt) < 1e-3) {
                                    vp.pop_back();
                                    t = vp.back();
                                } else {
                                    t = a;
                                }

                            } else {
                                Vert _v(x, d->s, ref);
                                verts.push_back(_v);

                                int k = verts.size()-1;

                                verts[a].nxt = k;

                                vp.push_back(k);

                                t = k;

                            }

                            break;

                        }

                    }

                    int p = v;

                    int w = verts[v].nxt;

                    if (Ld(x, ptV, verts[w].pt) < 1e-3) {
                        p = w;
                        w = verts[w].nxt;
                    }

                    std::cout << v << " -> " << p << std::endl;

                    double *ptW = verts[w].pt,
                        *ptP = verts[p].pt;

                    double cC = Cross(x, ptV, ptW),
                        cD = Cross(ptV, ptU, ptW);

                    std::cout << "cC " << cC << std::endl;
                    std::cout << "cD " << cD << std::endl;

                    if (cC < 0) {

                        if (cD < 0) {
                            verts[vp.back()].nxt = p;

                            vp.push_back(p);
                        } else {
                            int _x = w;

                            for (;;) {
                                Vert &a = verts[_x],
                                    &b = verts[a.nxt];

                                double *ptA = a.pt,
                                    *ptB = b.pt;

                                std::cout << ">> " << a.id << ", " << b.id << std::endl;

                                std::shared_ptr<D> d(Intersect(x, verts[v].r, ptA, ptB));

                                if (d && !IsFrontfaced(verts[v].r, ptA, ptB)) {
                                    std::cout << "x" << std::endl;

                                    if (d->t2 < E) {
                                        verts[vp.back()].nxt = _x;
                                        vp.push_back(_x);

                                    } else {
                                        Vert _v(x, d->s, ref, a.nxt);
                                        verts.push_back(_v);
                                        int k = verts.size()-1;

                                        verts[vp.back()].nxt = k;

                                        vp.push_back(k);

                                    }

                                    break;
                                } else {
                                    _x = a.nxt;
                                }

                            }
                        }

                    } else {
                        verts[vp.back()].nxt = p;

                        vp.push_back(p);
                    }

                }
            }
        }
    }

    res.push_back(poly[ind]);

    for (int _v : vp) {
        int id = verts[_v].id;
        res.push_back({verts[_v].pt, id != NO_USE ? poly[id].id : id});
    }

}

void GetVisPoly_wrapper (PolyType &poly, PolyType &res, int ind) {
    PolyType poly2(poly);

    PolyType poly3;
    TrivialRm(poly2, ind).GetSimplified(res);

    //GetVisPoly(poly3, res);
}
