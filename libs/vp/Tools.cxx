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

#include <cmath>
#include <cassert>
#include <cfloat>

#include "Tools.h"

double Normalize (double *v, const int n) {
    double l;

    // pow(0, -.5) liefert inf
    // 1/inf -> 0

    if (n == 3) {
        l = std::pow(v[0]*v[0]+v[1]*v[1]+v[2]*v[2], -.5);
        v[0] *= l;
        v[1] *= l;
        v[2] *= l;

    } else {
        l = std::pow(v[0]*v[0]+v[1]*v[1], -.5);
        v[0] *= l;
        v[1] *= l;
    }

    return 1/l;
}

double GetAngle (double *vA, double *vB) {
    double ang = std::atan2(vA[0]*vB[1]-vB[0]*vA[1], vA[0]*vB[0]+vA[1]*vB[1]);
    if (ang < 0) {
        ang += 2*PI;
    }
    return ang;
}

bool Ld (double *a, double *b, double *c, double *r) {
    double vA[] = {b[0]-a[0], b[1]-a[1]},
        vB[] = {c[0]-a[0], c[1]-a[1]};

    if (r != nullptr) {
        *r = (vA[0]*vA[0]+vA[1]*vA[1])/(vB[0]*vB[0]+vB[1]*vB[1]);
    }

    return std::abs(vA[0]*vB[1]-vA[1]*vB[0]) < E;
}

double Cross (double *a, double *b, double *c) {
    // kreuzprodukt der vektoren ab und ac
    return (b[1]-a[1])*(c[0]-a[0])-(b[0]-a[0])*(c[1]-a[1]);
}

std::shared_ptr<D> Intersect (const double *o, const double *r, const double *pA, const double *pB) {
    double oB[] = {o[0]+r[0], o[1]+r[1]};

    double m11 = oB[0]-o[0],
        m12 = pA[0]-pB[0],
        m21 = oB[1]-o[1],
        m22 = pA[1]-pB[1],
        v1 = pA[0]-o[0],
        v2 = pA[1]-o[1];

    double det = m11*m22-m12*m21;

    if (std::abs(det) < E) {
        return nullptr;
    }

    double t1 = (v1*m22-m12*v2)/det,
        t2 = (m11*v2-v1*m21)/det;

    if (t2 > -E && t2 < 1-E) {
        double s[] = {pA[0]+t2*(pB[0]-pA[0]), pA[1]+t2*(pB[1]-pA[1])};
        return std::make_shared<D>(s, t1, t2);
    } else {
        return nullptr;
    }
}

std::shared_ptr<D> Intersect2 (const double *oA, const double *oB, const double *pA, const double *pB, const Bnds &bnds) {
    double m11 = oB[0]-oA[0],
        m12 = pA[0]-pB[0],
        m21 = oB[1]-oA[1],
        m22 = pA[1]-pB[1],
        v1 = pA[0]-oA[0],
        v2 = pA[1]-oA[1];

    double det = m11*m22-m12*m21;

    if (std::abs(det) < E) {
        return nullptr;
    }

    double t1 = (v1*m22-m12*v2)/det,
        t2 = (m11*v2-v1*m21)/det;

    // if ((t1 > -E && t1 < 1+E)
    //     && (t2 > E && t2 < 1+E)) {

    if (t1 > bnds.lA && t1 < bnds.uA
        && t2 > bnds.lB && t2 < bnds.uB) {

        double s[] = {pA[0]+t2*(pB[0]-pA[0]), pA[1]+t2*(pB[1]-pA[1])};
        return std::make_shared<D>(s, t1, t2);
    } else {
        return nullptr;
    }
}

std::shared_ptr<D> Intersect2 (const double *oA, const double *oB, const double *pA, const double *pB) {
    const Bnds bnds;
    return Intersect2(oA, oB, pA, pB, bnds);
}

bool IsFrontfaced (double *r, double *a, double *b) {
    double rAB[] = {a[1]-b[1], b[0]-a[0]}; // um pi/2 gedreht
    Normalize(rAB);

    double _x = r[0]*rAB[0]+r[1]*rAB[1];

    //return std::acos(r[0]*rAB[0]+r[1]*rAB[1]) > PI/2;
    return _x < 0;
}

bool IsNear (const double *a, const double *b) {
    return std::abs(b[0]-a[0]) < E
        && std::abs(b[1]-a[1]) < E;
}

double GetT (double *a, double *b, double *c) {
    if (IsNear(b, c)) {
        return 1;
    }

    double vA[] = {b[0]-a[0], b[1]-a[1]},
        vB[] = {c[0]-a[0], c[1]-a[1]},
        lA = Normalize(vA),
        lB = Normalize(vB);

    return (vA[0]*vB[0]+vA[1]*vB[1])*lB/lA;
}

bool TestCW (const PolyType &poly) {
    // http://mathworld.wolfram.com/PolygonArea.html

    int num = poly.size();

    double sum = 0;

    for (int i = 0; i < num; i++) {
        const Point &a = poly[i],
            &b = poly[(i+1)%num];
        sum += a.x*b.y-b.x*a.y;
    }

    return sum < 0;
}

std::string GetAbsolutePath (const PolyType &poly) {
    std::stringstream path;

    for (const Point& p : poly) {
        path << p.x << "," << p.y << " ";
    }

    std::string svg = "M" + path.str() + "Z";

    return svg;
}

bool TestPIP (PolyType &poly, Point &pt) {
    // Point-In-Polygon

    for (auto& p : poly) {
        if (IsNear(p.pt, pt.pt)) {
            return false;
        }
    }

    int num = poly.size();

    bool in = false;

    for (int i = 0; i < num; i++) {
        Point &a = poly[i],
            &b = poly[(i+1)%num];

        if ((a.x <= pt.x || b.x <= pt.x)
            && (a.y < pt.y && b.y >= pt.y
                || b.y < pt.y && a.y >= pt.y)) {

            // schnittpunkt mit bounding box und strahlensatz
            if (a.x+(pt.y-a.y)*(b.x-a.x)/(b.y-a.y) < pt.x) {
                in = !in;
            }
        }
    }

    return in;
}

void GetExt (const PolyType &poly, Ext &ext) {
    for (auto &p : poly) {
        if (p.x < ext.minX) {
            ext.minX = p.x;
        }
        if (p.x > ext.maxX) {
            ext.maxX = p.x;
        }
        if (p.y < ext.minY) {
            ext.minY = p.y;
        }
        if (p.y > ext.maxY) {
            ext.maxY = p.y;
        }
    }
}

double GetArea (const PolyType &poly) {
    int num = poly.size();

    double sum = 0;

    for (int i = 0; i < num; i++) {
        const Point &a = poly[i],
            &b = poly[(i+1)%num];
        sum += a.x*b.y-b.x*a.y;
    }

    return std::abs(sum);
}

double GetDis (const Point &pA, const Point &pB, const Point &pC, double &t, double *pro) {
    double n[] = {pA.y-pB.y, pB.x-pA.x},
        l = Normalize(n),
        d = n[0]*(pA.x-pC.x)+n[1]*(pA.y-pC.y),
        v[] = {pC.x+d*n[0], pC.y+d*n[1]},
        w[] = {v[0]-pA.x, v[1]-pA.y},
        m = Normalize(w);

    t = (n[1]*w[0]-n[0]*w[1])*m/l;

    if (pro != nullptr) {
        Cpy(pro, v);
    }

    return std::abs(d);
}

double GetSqDis (const Point &a, const Point &b) {
    double v[] = {
        b.x-a.x,
        b.y-a.y
    };
    return v[0]*v[0]+v[1]*v[1];
}

void GetSect (int tagA, int tagB, PolyType &poly) {
    std::rotate(poly.begin(), std::find_if(poly.begin(), poly.end(), [&tagA](const Point &p) {
        return p.tag == tagA;
    }), poly.end());

    poly.erase(std::find_if(poly.begin()+1, poly.end(), [&tagB](const Point &p) {
        return p.tag == tagB;
    })+1, poly.end());
}

#ifdef _WIN32
int Point::_tag;
#endif
