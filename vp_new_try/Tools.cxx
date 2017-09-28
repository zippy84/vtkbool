#include <cmath>

#include "Tools.h"

double Normalize (double *v) {
    double l = std::pow(v[0]*v[0]+v[1]*v[1], -.5);
    v[0] *= l;
    v[1] *= l;

    return 1/l;
}

double GetAngle (double *vA, double *vB) {
    double ang = std::atan2(vA[0]*vB[1]-vB[0]*vA[1], vA[0]*vB[0]+vA[1]*vB[1]);
    if (ang < 0) {
        ang += 2*PI;
    }
    return ang;
}

double Ld (double *a, double *b, double *c) {
    return std::abs(a[0]*(b[1]-c[1])-a[1]*(b[0]-c[0])+(b[0]*c[1]-c[0]*b[1]));
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
    return std::acos(r[0]*rAB[0]+r[1]*rAB[1]) > PI/2;
}

bool IsNear (double *a, double *b) {
    return std::abs(b[0]-a[0]) < E
        && std::abs(b[1]-a[1]) < E;
}

double GetT (double *a, double *b, double *c) {
    double m11 = a[0],
        m12 = b[0]-a[0],
        m21 = a[1],
        m22 = b[1]-a[1],
        v1 = c[0],
        v2 = c[1];

    double det = m11*m22-m12*m21;

    return (m11*v2-v1*m21)/det;
}

bool IsOnSeg (double *a, double *b, double *c) {
    if (IsNear(a, c) || IsNear(b, c)) {
        return false;
    }

    if ((c[0] < a[0] && c[0] < b[0])
        || (c[0] > a[0] && c[0] > b[0])) {
        return false;
    }

    if ((c[1] < a[1] && c[1] < b[1])
        || (c[1] > a[1] && c[1] > b[1])) {
        return false;
    }

    return std::abs(Cross(a, b, c)) < E;
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

bool TestPIP (const PolyType &poly, const Point &pt) {
    // Point-In-Polygon

    int num = poly.size();

    bool in = false;

    for (int i = 0; i < num; i++) {
        const Point &a = poly[i],
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
