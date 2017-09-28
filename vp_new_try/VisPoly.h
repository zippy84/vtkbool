#ifndef __VisPoly_h
#define __VisPoly_h

#include <vector>
#include <ostream>

#include "Tools.h"

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

// diese darf nicht direkt verwendet werden
void GetVisPoly (PolyType &poly, PolyType &res, int ind = 0);

void GetVisPoly_wrapper (PolyType &poly, PolyType &res, int ind);

#endif
