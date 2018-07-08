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

#ifndef __VisPoly_h
#define __VisPoly_h

#include <vector>
#include <ostream>
#include <exception>

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

class vp_error : public std::exception {
public:
    const char* what() const throw() {
        return "Too many pop's.";
    }
};

// diese darf nicht direkt verwendet werden
void GetVisPoly (PolyType &poly, PolyType &res, int ind = 0);

bool GetVisPoly_wrapper (PolyType &poly, PolyType &res, int ind);

#endif
