#ifndef __VisPoly_h
#define __VisPoly_h

#include <vector>
#include <ostream>

#include "Tools.h"

class Point {
public:
    Point (double _x, double _y, int _id = NO_USE) : id(_id) {
        pt[0] = _x;
        pt[1] = _y;
    }
    Point (double *_pt, int _id = NO_USE) : Point(_pt[0], _pt[1], _id) {}

    Point (const Point& p) : id(p.id) {
        pt[0] = p.pt[0];
        pt[1] = p.pt[1];
    }

    Point& operator= (const Point &p) {
        pt[0] = p.pt[0];
        pt[1] = p.pt[1];
        id = p.id;

        return *this;
    }

    double pt[2];

    const double &x = pt[0],
        &y = pt[1];
    int id;

    friend std::ostream& operator<< (std::ostream &out, const Point &p) {
        out << "id: " << p.id << ", pt: [" << p.x << "," << p.y << "]";
        return out;
    }
};

typedef std::vector<Point> PolyType;

class Vert {
public:
    double pt[2], r[2], phi;
    int id, nxt;

    Vert (double *_x, double *_pt, int _id) : id(_id), nxt(NO_USE) {
        pt[0] = _pt[0];
        pt[1] = _pt[1];

        r[0] = pt[0]-_x[0];
        r[1] = pt[1]-_x[1];

        Normalize(r);
    }
    Vert (double *_x, double *_pt, double *ref, int _nxt = NO_USE) : Vert(_x, _pt, NO_USE) {
        nxt = _nxt;
        phi = GetAngle(ref, r);
    }

};

typedef std::vector<Vert> VertsType;

class Bag {
public:
    Bag (int _f, int _g, double _phi) : f(_f), g(_g), phi(_phi) {

    }

    int f, g;
    double phi;

    friend std::ostream& operator<< (std::ostream &out, const Bag &b) {
        out << "phi: " << b.phi
            << ", g: " << b.g
            << ", f: " << b.f;
        return out;
    }
};

void GetVisPoly (PolyType &poly, PolyType &res, int ind = 0);

#endif
