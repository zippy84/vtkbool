/*
   Copyright 2012-2016 Ronald Römer

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

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cfloat>
#include <stack>
#include <sstream>
#include <set>
#include <deque>
#include <map>
#include <cassert>

#include "Utilities.h"
using Utilities::GetD;
using Utilities::Normalize;
using Utilities::GetNormal;
using Utilities::IdsType;
using Utilities::Pair;

#include "Decomposer.h"

namespace Decomposer {

    Base::Base (double pts[][3], const int num) {
        GetNormal(pts, n, num);

        ei[0] = pts[1][0]-pts[0][0];
        ei[1] = pts[1][1]-pts[0][1];
        ei[2] = pts[1][2]-pts[0][2];

        Normalize(ei);

        ej[0] = n[1]*ei[2]-n[2]*ei[1];
        ej[1] = -n[0]*ei[2]+n[2]*ei[0];
        ej[2] = n[0]*ei[1]-n[1]*ei[0];

        Normalize(ej);

        d = n[0]*pts[0][0]+n[1]*pts[0][1]+n[2]*pts[0][2];
    }

    void Transform (double *in, double *out, Base &base) {
        double x = base.ei[0]*in[0]+base.ei[1]*in[1]+base.ei[2]*in[2],
            y = base.ej[0]*in[0]+base.ej[1]*in[1]+base.ej[2]*in[2],
            z = base.n[0]*in[0]+base.n[1]*in[1]+base.n[2]*in[2];

        out[0] = x;
        out[1] = y;
    }

    void Transform (double in[][3], double out[][2], const int num, Base &base) {
        for (int i = 0; i < num; i++) {
            Transform(in[i], out[i], base);
        }
    }

    bool Intersect (double *o, double *r, double *pA, double *pB, double *s, double *t, double *u) {
        // schnittpunkt einer geraden mit einer linie

        double *oA = o,
            oB[] = {o[0]+r[0], o[1]+r[1]};

        double m11 = oB[0]-oA[0],
            m12 = pA[0]-pB[0],
            m21 = oB[1]-oA[1],
            m22 = pA[1]-pB[1],
            v1 = pA[0]-oA[0],
            v2 = pA[1]-oA[1];

        double det = m11*m22-m12*m21;

        double t1 = (v1*m22-m12*v2)/det,
            t2 = (m11*v2-v1*m21)/det;

        if (t2 > -1e-6 && t2 < 0.999999) {
            s[0] = pA[0]+t2*(pB[0]-pA[0]);
            s[1] = pA[1]+t2*(pB[1]-pA[1]);

            if (t != NULL) {
                *t = t1;
            }

            if (u != NULL) {
                *u = t2;
            }

            return true;
        } else {
            return false;
        }

    }

    bool Intersect2 (double *oA, double *oB, double *pA, double *pB, double *s) {
        // schnittpunkt zweier linien

        double m11 = oB[0]-oA[0],
            m12 = pA[0]-pB[0],
            m21 = oB[1]-oA[1],
            m22 = pA[1]-pB[1],
            v1 = pA[0]-oA[0],
            v2 = pA[1]-oA[1];

        double det = m11*m22-m12*m21;

        double t1 = (v1*m22-m12*v2)/det,
            t2 = (m11*v2-v1*m21)/det;

        if (t1 > 0 && t1 < 1 && t2 > 0 && t2 < 1) {
            s[0] = pA[0]+t2*(pB[0]-pA[0]);
            s[1] = pA[1]+t2*(pB[1]-pA[1]);

            return true;
        } else {
            return false;
        }

    }

    double GetAng (double *vA, double *vB) {
        // http://math.stackexchange.com/questions/878785/how-to-find-an-angle-in-range0-360-between-2-vectors
        double ang = std::atan2(vA[0]*vB[1]-vB[0]*vA[1], vA[0]*vB[0]+vA[1]*vB[1]);

        if (ang < 0) {
            ang += 2*pi;
        }

        return ang;
    }

    void GetVisiblePoly (double pts[][2], const int num, std::vector<Vert2> &poly, int ind) {
        // pts ist ccw

        poly.clear();
        std::vector<Vert> _poly;

        double _ei[] = {-1, 0},
            x[] = {pts[ind][0], pts[ind][1]};

        // startpunkt
        std::vector<Cut> cuts, cuts2;

        for (int i = 0; i < num; i++) {
            double *ptI = pts[i],
                *ptJ = pts[(i+1)%num];

            if (i == ind) {
                cuts.push_back(Cut(ind, x, 0, 0));
            } else {
                double s[2], t, u;
                if (Intersect(x, _ei, ptI, ptJ, s, &t, &u)) {
                    Cut cut(i, s, t, u);
                    cuts.push_back(std::move(cut));
                }
            }
        }

        std::vector<Cut>::const_iterator itr;

        /*
        std::cout << "cuts=[";
        for (itr = cuts.begin(); itr != cuts.end(); ++itr) {
            std::cout << "(" << itr->i << ", " << itr->t << "), ";
        }
        std::cout << "]" << std::endl;
        */

        std::sort(cuts.begin(), cuts.end(), [](const Cut &a, const Cut &b) {
            if (std::abs(b.t-a.t) < 1e-6) {
                return false;
            } else {
                return a.t < b.t;
            }
        });

        /*
        std::cout << "cuts'=[";
        for (itr = cuts.begin(); itr != cuts.end(); ++itr) {
            std::cout << "(" << itr->i << ", " << itr->t << "), ";
        }
        std::cout << "]" << std::endl;
        */

        auto res2 = std::find_if(cuts.begin(), cuts.end(), [ind](const Cut &c) { return c.i == ind; });

        cuts2.insert(cuts2.begin(), res2, cuts.end());

        /*
        std::cout << "cuts2=[";
        for (itr = cuts2.begin(); itr != cuts2.end(); ++itr) {
            std::cout << "(" << itr->i << ", " << itr->t << "), ";
        }
        std::cout << "]" << std::endl;
        */

        std::remove_if(cuts2.begin()+1, cuts2.end(), [&](const Cut &c) {
            if (c.u < 1e-6) {
                double a = pts[(c.i+1)%num][1]-x[1],
                    b = pts[(c.i+num-1)%num][1]-x[1];
                return (a < -1e-6 && b < -1e-6) || (a > 1e-6 && b > 1e-6);
            } else {
                return false;
            }
        });

        /*
        for (int i= 0; i < num; i++) {
            std::cout << (i == 0 ? "M" : "L") << pts[i][0] << "," << pts[i][1] << " ";
        }
        std::cout << "Z" << std::endl;
        */

        assert(cuts2[0].i == ind);

        double *xA = pts[(ind+1)%num],
            *xB = pts[(ind+num-1)%num];

        double a[] = {xA[0]-x[0], xA[1]-x[1]},
            b[] = {xB[0]-x[0], xB[1]-x[1]};

        double angA = GetAng(a, b),
            angB = GetAng(a, _ei);

        /*
        std::cout << (angA*180/pi) << ", " << (angB*180/pi) << " => " << (angB < angA) << std::endl;
        */

        int f = 0;
        if (cuts2.size() > 1 && angB < angA) {
           f = 1;
        }

        Cut &cutA = cuts2[f];

        int o = 0;
        if (cutA.u > 1e-6) {
            _poly.push_back(Vert(x, cutA.s, 1));
            o++;
        }

        for (int i = 0; i < num; i++) {
            int k = (cutA.i+o+i)%num;
            Vert vert(x, pts[k], (i+o+1)%(num+o), k);

            /*
            std::cout << i << ": " << vert.ind << ", " << (vert.phi*180/pi) << ", r=["
                << vert.r[0] << ", " << vert.r[1] << "]" << std::endl;
            */

            _poly.push_back(std::move(vert));
        }

        /*
        std::cout << "f=" << f << ", o=" << o << std::endl;

        int i = 0;
        do {
            std::cout << (i == 0 ? "M" : "L") << _poly[i].pt[0] << "," << _poly[i].pt[1] << " ";
        } while ((i = _poly[i].succ) > 0);
        std::cout << "Z" << std::endl;
        */

        std::vector<int> vp;
        vp.push_back(0);
        vp.push_back(1);

        std::stack<Bag> lBags, rBags;

        if (f == 1) {
            double p[] = {-9999999, x[1]};
            Vert vP(x, p, NOUSE);
            _poly.push_back(vP);
            rBags.push(Bag(0, _poly.size()-1, 0));
        }

        std::vector<int>::iterator itr2;

        double tol = 1e-6;

        int t = 0, u, v;

        for (;;) {

            Vert &tVert = _poly[t],
                &uVert = _poly[(u = tVert.succ)],
                &vVert = _poly[(v = uVert.succ)];

            /*
            std::cout << "-- " << t << ": (" << u << ", " << v << ")" << std::endl;
            */

            double tPt[3], uPt[3], vPt[3];
            CPY(tPt, tVert.pt)
            CPY(uPt, uVert.pt)
            CPY(vPt, vVert.pt)

            double dT = GetD(tPt, vPt);

            double cA = (uPt[1]-tPt[1])*(vPt[0]-tPt[0])-(uPt[0]-tPt[0])*(vPt[1]-tPt[1]),
                cB = (uPt[1]-x[1])*(vPt[0]-x[0])-(uPt[0]-x[0])*(vPt[1]-x[1]);

            /*
            std::cout << "cA=" << cA << " cB=" << cB << std::endl;
            */

            if (v == 0) {
                break;
            }

            if (cB < -tol) {
                Bag *bag = NULL;
                while (!rBags.empty() && rBags.top().phi < vVert.phi) {
                    Bag &top = rBags.top();
                    if (bag == NULL) {
                        bag = new Bag(top);
                    } else {
                        *bag = top;
                    }

                    /*
                    std::cout << "pop bag" << std::endl;
                    */
                    rBags.pop();
                }

                double k[2];
                if (bag != NULL && Intersect2(uPt, vPt, _poly[bag->f].pt, _poly[bag->g].pt, k)) {
                    /*
                    std::cout << "Case 1(a)" << std::endl;
                    */

                    int w = v;
                    for (;;) {
                        Vert &vA = _poly[w],
                            &vB = _poly[vA.succ];

                        double l[2];
                        if (Intersect2(_poly[bag->f].pt, k, vA.pt, vB.pt, l)) {
                            Vert vL(x, l, vA.succ);
                            _poly.push_back(std::move(vL));
                            int indL = _poly.size()-1;

                            Vert vK(x, k, indL);
                            _poly.push_back(std::move(vK));
                            int indK = _poly.size()-1;

                            //uVert.succ = indK;
                            _poly[u].succ = indK;

                            vp.push_back(indK);
                            vp.push_back(indL);

                            rBags.push(Bag(bag->f, indL, bag->phi));

                            /*
                            std::cout << "add M" << _poly[bag->f].pt[0]
                                << "," << _poly[bag->f].pt[1]
                                << " L" << l[0]
                                << "," << l[1] << std::endl;
                            */

                            t = indK;

                            break;
                        }

                        w = vA.succ;
                        if (w == 0) {
                            break;
                        }

                    }

                } else {
                    /*
                    std::cout << "Case 1(b)" << std::endl;
                    */

                    if (bag != NULL) {
                        /*
                        std::cout << "re-push bag" << std::endl;
                        */
                        rBags.push(*bag);
                    }

                    vp.push_back(v);
                    t = u;
                }

                if (bag != NULL) {
                    delete bag;
                }

            } else if ((cA > -tol && cB > tol) || dT < tol) {
                /*
                std::cout << "Case 2" << std::endl;
                */

                int w = v;
                for (;;) {
                    Vert &vA = _poly[w],
                        &vB = _poly[vA.succ];

                    /*
                    std::cout << "edge=[" << w << "," << vA.succ << "]" << std::endl;
                    */

                    double k[2], q;
                    if (Intersect(x, uVert.r, vA.pt, vB.pt, k, &q) && q > tol) {
                        Vert vK(x, k, vA.succ);
                        _poly.push_back(std::move(vK));
                        int indK = _poly.size()-1;

                        //uVert.succ = indK;
                        // durch push_back ist die ref. uVert nicht mehr gültig
                        _poly[u].succ = indK;

                        lBags.push(Bag(u, indK, _poly[u].phi));

                        /*
                        std::cout << "add M" << _poly[u].pt[0]
                            << "," << _poly[u].pt[1]
                            << " L" << k[0]
                            << "," << k[1] << std::endl;
                        */

                        vp.push_back(indK);

                        t = u;

                        break;
                    }

                    w = vA.succ;
                    if (w == 0) {
                        // schnitt ist genau auf 0
                        uVert.succ = 0;
                        break;
                    }
                }

            } else if (cA < 0 && cB > tol) {
                Bag *bag = NULL;
                while (!lBags.empty() && lBags.top().phi > vVert.phi) {
                    Bag &top = lBags.top();
                    if (bag == NULL) {
                        bag = new Bag(top);
                    } else {
                        *bag = top;
                    }

                    /*
                    std::cout << "pop bag" << std::endl;
                    */
                    lBags.pop();
                }

                double k[2];
                if (bag != NULL && Intersect2(uPt, vPt, _poly[bag->f].pt, _poly[bag->g].pt, k)) {
                    /*
                    std::cout << "Case 3(a)" << std::endl;
                    */

                    int w = v;
                    for (;;) {
                        Vert &vA = _poly[w],
                            &vB = _poly[vA.succ];

                        double l[2];
                        if (Intersect2(_poly[bag->f].pt, k, vA.pt, vB.pt, l)) {
                            int m;
                            do {
                                m = vp.back();
                                vp.pop_back();
                            } while (m != bag->g);

                            Vert vL(x, l, vA.succ);
                            _poly.push_back(std::move(vL));
                            int indL = _poly.size()-1;

                            _poly[bag->f].succ = indL;

                            vp.push_back(indL);
                            lBags.push(Bag(bag->f, indL, bag->phi));

                            /*
                            std::cout << "add M" << _poly[bag->f].pt[0]
                                << "," << _poly[bag->f].pt[1]
                                << " L" << l[0]
                                << "," << l[1] << std::endl;
                            */

                            t = bag->f;

                            break;
                        }

                        w = vA.succ;
                        if (w == 0) {
                            break;
                        }

                    }

                } else {
                    /*
                    std::cout << "Case 3(b)" << std::endl;
                    */

                    if (bag != NULL) {
                        /*
                        std::cout << "re-push bag" << std::endl;
                        */
                        lBags.push(*bag);
                    }

                    while (vp.size() > 1) {
                        int pre = *(vp.end()-2);
                        
                        Vert &vA = _poly[pre],
                            &vB = _poly[vA.succ];

                        /*
                        std::cout << "pop vert " << vp.back() << std::endl;
                        */
                        vp.pop_back();

                        double k[2], q, r;
                        if (Intersect(x, vVert.r, vA.pt, vB.pt, k, &q, &r) && q > tol) {
                            int indK;

                            if (r > tol) {
                                Vert vK(x, k, v);
                                _poly.push_back(std::move(vK));
                                indK = _poly.size()-1;

                                _poly[pre].succ = indK;
                                vp.push_back(indK);

                                /*
                                std::cout << "(1)" << std::endl;
                                */

                            } else {
                                // k ist vA
                                vA.succ = v;
                                indK = vp.back();

                                /*
                                std::cout << "(2)" << std::endl;
                                */
                            }

                            t = indK;

                            // zweiter teil

                            Vert &_vVert = _poly[v];
                            Vert &wVert = _poly[_vVert.succ];

                            double *wPt = wVert.pt;

                            double cC = (vPt[1]-x[1])*(wPt[0]-x[0])-(vPt[0]-x[0])*(wPt[1]-x[1]),
                                cD = (vPt[1]-uPt[1])*(wPt[0]-uPt[0])-(vPt[0]-uPt[0])*(wPt[1]-uPt[1]);

                            /*
                            std::cout << "cC=" << cC << " cD=" << cD << std::endl;
                            */

                            if (cC > -tol) {
                                /*
                                std::cout << "Reg a" << std::endl;
                                */

                                vp.push_back(v);
                            } else if (cC < 0 && cD > -tol) {
                                /*
                                std::cout << "Reg b" << std::endl;
                                */

                                vp.push_back(v);
                                lBags.push(Bag(v, indK, _vVert.phi));

                                /*
                                std::cout << "add M" << _poly[v].pt[0]
                                    << "," << _poly[v].pt[1]
                                    << " L" << k[0]
                                    << "," << k[1] << std::endl;
                                */

                            } else if (cC < 0 && cD < 0) {
                                /*
                                std::cout << "Reg c" << std::endl;
                                */

                                int w = _vVert.succ;
                                for (;;) {
                                    Vert &vA = _poly[w],
                                        &vB = _poly[vA.succ];

                                    double l[2];
                                    if (Intersect2(k, vPt, vA.pt, vB.pt, l)) {
                                        Vert vL(x, l, vA.succ);
                                        _poly.push_back(std::move(vL));
                                        int indL = _poly.size()-1;

                                        (_poly.end()-2)->succ = indL;

                                        vp.push_back(indL);
                                        rBags.push(Bag(v, indL, _vVert.phi));

                                        /*
                                        std::cout << "add M" << _poly[v].pt[0]
                                            << "," << _poly[v].pt[1]
                                            << " L" << l[0]
                                            << "," << l[1] << std::endl;
                                        */

                                        break;
                                    }

                                    w = vA.succ;
                                    if (w == 0) {
                                        break;
                                    }

                                }

                            }

                            break;
                        }

                    }

                    // sonderfall

                    if (vp.size() == 1) {
                        vp.push_back(_poly[0].succ);
                        _poly[vp.back()].succ = v;
                        t = 0;
                    }

                }

                if (bag != NULL) {
                    delete bag;
                }

            } else {
                // die sonderfälle

                double dU = GetD(uPt, x),
                    dV = GetD(vPt, x);

                if (vVert.ind == ind) {
                    /*
                    std::cout << "Case 4(a)" << std::endl;
                    */

                    vp.push_back(v);
                    t = u;

                } else if (uVert.ind == ind) {
                    /*
                    std::cout << "Case 4(b)" << std::endl;
                    */

                    vp.push_back(v);
                    t = u;

                } else if (dV < tol) {
                    /*
                    std::cout << "Case 4(c)" << std::endl;
                    */

                    double ref[] = {uPt[0]-vPt[0], uPt[1]-vPt[1]};

                    double minAng = DBL_MAX;

                    int last = vVert.succ, m;

                    for (;;) {
                        Vert &next = _poly[last];

                        if (next.ind == ind) {
                            break;
                        }

                        assert(next.ind != -1);

                        double vec[] = {next.pt[0]-vPt[0], next.pt[1]-vPt[1]},
                            ang = GetAng(ref, vec);

                        /*
                        std::cout << next.ind << ": " << (ang*180/pi) << std::endl;
                        */

                        if (ang < minAng) {
                            minAng = ang;
                            m = last;
                        }

                        last = next.succ;
                    }

                    /*
                    std::cout << (minAng*180/pi) << " -> " << m << ", " << _poly[m].ind << std::endl;
                    */

                    while (vp.size() > 1) {
                        Vert &vA = _poly[*(vp.end()-2)],
                            &vB = _poly[vA.succ];

                        /*
                        std::cout << "pop vert " << vp.back() << std::endl;
                        */
                        vp.pop_back();

                        double k[2], q, r;
                        if (Intersect(x, _poly[m].r, vA.pt, vB.pt, k, &q, &r) && q > tol) {
                            int indK;

                            if (r > tol) {
                                Vert vK(x, k, m);
                                _poly.push_back(std::move(vK));
                                indK = _poly.size()-1;

                                vA.succ = indK;
                                vp.push_back(indK);

                            } else {
                                // k ist vA
                                vA.succ = m;
                                indK = vp.back();
                            }

                            vp.push_back(m);

                            lBags.push(Bag(m, indK, _poly[m].phi));

                            t = indK;

                            break;
                        }

                    }

                } else if (dU < tol) {
                    /*
                    std::cout << "Case 4(d)" << std::endl;
                    */

                    assert(true);
                } else {
                    // lin. abh.
                    /*
                    std::cout << "Case 4(e)" << std::endl;
                    */

                    vp.push_back(v);
                    t = u;
                }
            }

            /*
            for (itr2 = vp.begin(); itr2 != vp.end(); ++itr2) {                
                std::cout << (itr2-vp.begin() == 0 ? "M" : "L") << _poly[*itr2].pt[0] << "," << _poly[*itr2].pt[1] << " ";
            }
            std::cout << "Z" << std::endl;
            */

        }

        for (itr2 = vp.begin(); itr2 != vp.end(); ++itr2) {
            Vert2 v(_poly[*itr2]);
            poly.push_back(std::move(v));
        }

    }

    bool check (double pts[][2], int num) {
        int i = 0;
        while (i < num) {
            if (GetD(pts[i], pts[(i+1)%num]) < 1e-6) {
                break;
            }
            i++;
        }

        return i != num;
    }

    void SubD::update (int a, int b, int _w) {
        /*
        std::cout << "updating with w=" << _w
                  << " and pair " << a << ", " << b << std::endl;
        */

        if (_w <= w) {
            if (_w < w) {
                S.clear();
                w = _w;
            }

            if (!S.empty() && a <= S.getB().f) {
                /*
                std::cout << "skipped" << std::endl;
                */
                return;
            }

            while (!S.empty() && S.getB().g >= b) {
                /*
                std::cout << "popping (" << S.getB().f << ", " << S.getB().g << ")" << std::endl;
                */
                S.popB();
            }

            S.push(Pair(a, b));

            /*
            std::cout << "S=" << S << std::endl;
            */
        } else {
            /*
            std::cout << "current w is better" << std::endl;
            */
        }
    }

    bool Decompose::IsReflex (int i, int j, int k) {
        double vA[] = {pts[j][0]-pts[i][0], pts[j][1]-pts[i][1], 0},
            vB[] = {pts[k][0]-pts[i][0], pts[k][1]-pts[i][1], 0};

        return (vA[0]*vB[1]-vB[0]*vA[1] < -1e-6) || GetD(pts[j], pts[k]) < 1e-6;
    }

    void Decompose::TypeA (int i, int j, int k) {
        assert(refls[i]);

        /*
        std::cout << "TypeA(" << i << ", " << j << ", " << k << ")";
        */

        Pair pr = Pair(i, j);
        if (visVerts.count(pr) > 0) {
            /*
            std::cout << std::endl
                      << "i, j are visible" << std::endl;
            */

            SubD &sd = problems[pr];
            Deq &s = sd.S;
            int w = sd.w,
                a = j;

            /*
            std::cout << "w=" << w << std::endl;
            */

            // j=k-1 und j=i+1 wird ausgelassen
            // diese kanten hätten w=-1

            if (k-j > 1) {
                Pair pr_ = Pair(j, k);
                if (visVerts.count(pr_) == 0) {
                    /*
                    std::cout << "j, k not visible" << std::endl;
                    */
                    return;
                }
                w += problems[pr_].w+1;
            }

            if (j-i > 1) {
                // wenn die ersten beiden keine reflexe sind,
                // dann den letzten löschen und den anderen auf s lassen

                if (!IsReflex(j, k, s.getF().g)) {
                    while (s.size() > 1 && !IsReflex(j, k, s.getPreF().g)) {
                        /*
                        std::cout << "pop " << s.getF() << std::endl;
                        */
                        s.popF();
                    }

                    if (!s.empty() && !IsReflex(i, s.getF().f, k)) {
                        /*
                        std::cout << ">1" << std::endl;
                        */
                        a = s.getF().f;
                    } else {
                        /*
                        std::cout << ">2" << std::endl;
                        */
                        w++;
                    }
                } else {
                    // gleich der erste auf s war ein reflex
                    /*
                    std::cout << ">3" << std::endl;
                    */
                    w++;
                }
            }

            problems[Pair(i, k)].update(a, j, w);
        } else {
            /*
            std::cout << " -> not visible" << std::endl;
            */
        }

    }

    void Decompose::TypeB (int i, int j, int k) {
        assert(refls[k]);

        /*
        std::cout << "TypeB(" << i << ", " << j << ", " << k << ")";
        */

        Pair pr = Pair(j, k);
        if (visVerts.count(pr) > 0) {
            /*
            std::cout << std::endl
                      << "j, k are visible" << std::endl;
            */

            SubD &sd = problems[pr];
            Deq &s = sd.S;
            int w = sd.w,
                a = j;

            /*
            std::cout << "w=" << w << std::endl;
            */

            // j=k-1 und j=i+1 wird ausgelassen
            // diese kanten hätten w=-1

            if (j-i > 1) {
                Pair pr_ = Pair(i, j);
                if (visVerts.count(pr_) == 0) {
                    /*
                    std::cout << "j, k not visible" << std::endl;
                    */
                    return;
                }
                w += problems[pr_].w+1;
            }

            if (k-j > 1) {
                // wenn die ersten beiden keine reflexe sind,
                // dann den letzten löschen und den anderen auf s lassen

                if (!IsReflex(j, s.getB().f, i)) {
                    while (s.size() > 1 && !IsReflex(j, s.getPreB().f, i)) {
                        /*
                        std::cout << "pop " << s.getB() << std::endl;
                        */
                        s.popB();
                    }

                    if (!s.empty() && !IsReflex(k, i, s.getB().g)) {
                        /*
                        std::cout << ">1" << std::endl;
                        */
                        a = s.getB().g;
                    } else {
                        /*
                        std::cout << ">2" << std::endl;
                        */
                        w++;
                    }
                } else {
                    // gleich der erste auf s war ein reflex
                    /*
                    std::cout << ">3" << std::endl;
                    */
                    w++;
                }
            }

            problems[Pair(i, k)].update(j, a, w);
        } else {
            /*
            std::cout << " -> not visible" << std::endl;
            */
        }

    }

    void Decompose::run () {
        /*
        std::cout << "<path class=\"A\" d=\"";
        for (int i = 0; i < num; i++) {
            std::cout << (i == 0 ? "M" : "L") << pts[i][0] << "," << pts[i][1] << " ";
        }
        std::cout << "Z\"/>" << std::endl;
        */

        assert(!check(pts, num));


        std::vector<Vert2> poly;
        std::vector<Vert2>::const_iterator itr;

        for (int i = 0; i < num; i++) {
            bool r = IsReflex(i, (i+1)%num, (i+num-1)%num);
            refls.push_back(r);
        }

        isConcave = std::count(refls.begin(), refls.end(), true) > 0;

        if (isConcave) {
            // index 0 ist per konvention ein reflex
            refls[0] = true;

            for (int i = 0; i < num; i++) {
                if (refls[i]) {
                    /*
                    std::cout << "GetVisiblePoly(" << i << ")" << std::endl;
                    */

                    GetVisiblePoly(pts, num, poly, i);

                    for (itr = poly.begin()+1; itr != poly.end(); ++itr) {
                        if (itr->ind != NOUSE) {
                            if (i < itr->ind) {
                                visVerts.insert(Pair(i, itr->ind));
                            } else if (i > itr->ind) {
                                visVerts.insert(Pair(itr->ind, i));

                                // das entspricht der forderung: "if f is not reflex, then g must be"
                            }
                        }
                    }

                    int j = 0;

                    /*
                    std::cout << "<path class=\"B\" d=\"";
                    */
                    for (itr = poly.begin(); itr != poly.end(); ++itr) {
                        if (itr->ind != NOUSE) {
                            /*
                            std::cout << (j == 0 ? "M" : "L") << itr->pt[0] << "," << itr->pt[1] << " ";
                            */
                            j++;
                        }
                    }
                    /*
                    std::cout << "Z\"/>" << std::endl;
                    */
                }
            }

            /*
            std::cout << "init" << std::endl;
            */

            std::set<Pair>::iterator itr2;
            for (itr2 = visVerts.begin(); itr2 != visVerts.end(); ++itr2) {
                /*
                std::cout << *itr2 << std::endl;
                */

                if (!refls[itr2->f]) {
                    assert(refls[itr2->g]);
                }

                // w ist auf 9999 gesetzt

                SubD &sd = problems[*itr2];

                if (itr2->f == itr2->g-1) {
                    /*
                    std::cout << "edge" << std::endl;
                    */
                    sd.w = 0;
                } else if (itr2->f == itr2->g-2) {
                    sd.w = 0;
                    /*
                    std::cout << "triangle" << std::endl;
                    */
                    sd.S.push(Pair(itr2->f+1, itr2->f+1));
                }

            }

            /*
            std::cout << "decompose" << std::endl;
            */

            // prinzipieller ablauf

            for (int l = 3; l < num; l++) {
                /*
                std::cout << "--- l=" << l << std::endl;

                std::cout << "for(1)" << std::endl;
                */

                for (int i = 0; i < num-l; i++) {
                    if (refls[i]) {
                        int k = i+l;

                        if (visVerts.count(Pair(i, k)) > 0) {

                            if (refls[k]) {
                                for (int j = i+1; j < k; j++) {
                                    /*
                                    std::cout << "(1) ";
                                    */
                                    TypeA(i, j, k);
                                }
                            } else {
                                for (int j = i+1; j < k-1; j++) {
                                    if (refls[j]) {
                                        /*
                                        std::cout << "(2) ";
                                        */
                                        TypeA(i, j, k);
                                    } else {
                                        /*
                                        std::cout << "skip " << i << ", " << j << ", " << k << std::endl;
                                        */
                                    }
                                }

                                // in der schleife wird j = k-1 ausgelassen
                                /*
                                std::cout << "(3) ";
                                */
                                TypeA(i, k-1, k);
                            }
                        }
                    }
                }

                /*
                std::cout << "for(2)" << std::endl;
                */

                for (int k = l; k < num; k++) {
                    if (refls[k]) {
                        int i = k-l;

                        if (!refls[i] && visVerts.count(Pair(i, k)) > 0) {

                            // hier ist es genauso; in der nachstehenden schleife wird j = i+1 ausgelassen
                            /*
                            std::cout << "(1) ";
                            */
                            TypeB(i, i+1, k);

                            for (int j = i+2; j < k; j++) {
                                if (refls[j]) {
                                    /*
                                    std::cout << "(2) ";
                                    */
                                    TypeB(i, j, k);
                                } else {
                                    /*
                                    std::cout << "skip " << i << ", " << j << ", " << k << std::endl;
                                    */
                                }
                            }
                        }
                    }
                }

            }
        }  

    }

    void Decompose::Recover (int i, int k) {
        int j;
        
        if (k-i < 2) {
            return;
        }

        Deq &s = problems[Pair(i, k)].S;
        if (refls[i]) {
            j = s.getF().g;
            Recover(j, k);

            if (j-i > 1) {
                if (s.getF().f != s.getF().g) {
                    Deq &t = problems[Pair(i, j)].S;
                    t.restore();

                    while (!t.empty() && s.getF().f != t.getF().f) {
                        t.popF();
                    }

                    assert(!t.empty());
                }
                Recover(i, j);
            }
        } else {
            j = s.getB().f;
            Recover(i, j);

            if (k-j > 1) {
                if (s.getB().g != s.getB().f) {
                    Deq &t = problems[Pair(j, k)].S;
                    t.restore();

                    while (!t.empty() && s.getB().g != t.getB().g) {
                        t.popB();
                    }

                    assert(!t.empty());
                }
                Recover(j, k);
            }
        }
    }

    void Decompose::Draw (int i, int k) {
        int j; 
        bool realA, realB;
        
        if (k-i < 2) {
            return;
        }

        /*
        std::cout << "-- (" << i << ", " << k << ")" << std::endl;
        */

        Deq &s = problems[Pair(i, k)].S;
        
        if (refls[i]) {
            j = s.getF().g;
            realA = (s.getF().f == s.getF().g);
            realB = true;

            /*
            std::cout << "is reflex j=" << j << " " << s.getF() << std::endl;
            */
        } else {
            j = s.getB().f;
            realA = true;
            realB = (s.getB().g == s.getB().f);

            /*
            std::cout << "is not j=" << j << " " << s.getB() << std::endl;
            */
        }

        if (realA && j-i > 1) {
            /*
            std::cout << "draw diag i=" << i << ", j=" << j << std::endl;
            */
            diags.push_back(Pair(i, j));
        }

        if (realB && k-j > 1) {
            /*
            std::cout << "draw diag j=" << j << ", k=" << k << std::endl;
            */
            diags.push_back(Pair(j, k));
        }

        Draw(i, j);
        Draw(j, k);
    }

    void Decompose::collect (DecType &dec) {
        run();

        /*
        std::cout << "collect" << std::endl;
        */

        dec.clear();

        IdsType poly;
        dec.push_back(poly);

        if (isConcave) {

            PType::const_iterator itr;

            for (itr = problems.begin(); itr != problems.end(); ++itr) {
                const Pair &pr = itr->first;
                const SubD &sd = itr->second;

                /*
                std::cout << "Pair " << pr << " w=" << sd.w << " S=" << sd.S << std::endl;
                */
            }

            Recover(0, num-1);
            Draw(0, num-1);

            std::sort(diags.begin(), diags.end(), [](const Pair &a, const Pair &b) {
                const int aG = -a.g, bG = -b.g;
                return std::tie(a.f, aG) < std::tie(b.f, bG); });

            std::vector<Pair>::const_iterator itr2;
            for (itr2 = diags.begin(); itr2 != diags.end(); ++itr2) {
                /*
                std::cout << *itr2 << std::endl;
                */
            }

            // bildet einzelne polygone

            int i = 0, p = 0, q = 0;

            std::vector<int> ps, rs;

            while (i < num) {
                IdsType &c = dec[p];

                if (c.empty() || c.back() != i) {
                    c.push_back(i);
                }

                if (!rs.empty() && i == diags[rs.back()].g) {
                    if (c.front() != diags[rs.back()].f) {
                        c.push_back(diags[rs.back()].f);
                    }

                    rs.pop_back();

                    p = ps.back();
                    ps.pop_back();

                } else if (q < diags.size() && i == diags[q].f) {
                    if (c.front() != diags[q].g) {
                        c.push_back(diags[q].g);
                    }

                    dec.push_back(poly);

                    rs.push_back(q);
                    ps.push_back(p);

                    p = ++q;
                } else {
                    i++;
                }
            }

        } else {
            // polygon war bereits konvex

            for (int i = 0; i < num; i++) {
                dec[0].push_back(i);
            }
        }

        /*
        DecType::const_iterator itr3;
        IdsType::const_iterator itr4;
        for (itr3 = dec.begin(); itr3 != dec.end(); ++itr3) {
            std::cout << "[";
            for (itr4 = itr3->begin(); itr4 != itr->end(); ++itr4) {
                std::cout << *itr4 << ", ";
            }
            std::cout << "]" << std::endl;

            //assert(itr3->size() > 2);

            std::cout << "<path class=\"C\" d=\"";
            for (itr4 = itr3->begin(); itr4 != itr3->end(); ++itr4) {
                std::cout << (itr4-itr3->begin() == 0 ? "M" : "L") << pts[*itr4][0] << "," << pts[*itr4][1] << " ";
            }
            std::cout << "Z\"/>" << std::endl;
        }
        */

    }

}
