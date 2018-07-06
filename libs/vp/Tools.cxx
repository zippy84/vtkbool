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

void Move (double *a, double *b, double *c) {
    // damit a, b oder c nicht der ursprung ist
    // det(a,b,c) wäre dann auch 0

    double x = std::min({a[0], b[0], c[0]}),
        y = std::min({a[1], b[1], c[1]});

    double m[] = {1-x, 2-y};

    a[0] += m[0]; a[1] += m[1];
    b[0] += m[0]; b[1] += m[1];
    c[0] += m[0]; c[1] += m[1];
}

bool Ld (double *a, double *b, double *c) {
    double vA[] = {b[0]-a[0], b[1]-a[1]},
        vB[] = {c[0]-a[0], c[1]-a[1]};

    return std::abs(vA[0]*vB[1]-vA[1]*vB[0]) < E;
}

double Cross (double *a, double *b, double *c) {
    // kreuzprodukt der vektoren ab und ac
    return (b[1]-a[1])*(c[0]-a[0])-(b[0]-a[0])*(c[1]-a[1]);
}

std::shared_ptr<D> Intersect (double *o, double *r, double *pA, double *pB) {
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

std::shared_ptr<D> Intersect2 (double *oA, double *oB, double *pA, double *pB) {
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

    if ((t1 > -E && t1 < 1+E)
        && (t2 > E && t2 < 1+E)) {

        double s[] = {pA[0]+t2*(pB[0]-pA[0]), pA[1]+t2*(pB[1]-pA[1])};
        return std::make_shared<D>(s, t1, t2);
    } else {
        return nullptr;
    }
}

bool IsFrontfaced (double *r, double *a, double *b) {
    double rAB[] = {a[1]-b[1], b[0]-a[0]}; // um pi/2 gedreht
    Normalize(rAB);

    double _x = r[0]*rAB[0]+r[1]*rAB[1];

    //return std::acos(r[0]*rAB[0]+r[1]*rAB[1]) > PI/2;
    return _x < 0;
}

bool IsNear (double *a, double *b) {
    return std::abs(b[0]-a[0]) < E
        && std::abs(b[1]-a[1]) < E;
}

double GetT (double *a, double *b, double *c) {
    if (IsNear(b, c)) {
        return 1;
    }

    double _a[2], _b[2], _c[2];
    Cpy(_a, a);
    Cpy(_b, b);
    Cpy(_c, c);

    Move(_a, _b, _c);

    double m11 = _a[0],
        m12 = _b[0]-_a[0],
        m21 = _a[1],
        m22 = _b[1]-_a[1],
        v1 = _c[0],
        v2 = _c[1];

    double det = m11*m22-m12*m21;

    assert(std::abs(det) > E);

    return (m11*v2-v1*m21)/det;
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
        path << "L" << p.x << "," << p.y << " ";
    }

    std::string svg = "M" + path.str().substr(1) + "Z";

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
