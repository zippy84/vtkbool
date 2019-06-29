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

#ifndef __VisPoly_h
#define __VisPoly_h

#include <vector>
#include <map>
#include <set>
#include <string>
#include <ostream>
#include <cassert>
#include <exception>

#include "Tools.h"

class Vert : public Point {
public:
    double r[2], phi;
    int nxt;

    Vert (double *_x, Point &p) : Point(p), nxt(NO_USE) {
        r[0] = pt[0]-_x[0];
        r[1] = pt[1]-_x[1];

        Normalize(r);
    }
    Vert (double *_x, Point &p, double *ref, int _nxt = NO_USE) : Vert(_x, p) {
        nxt = _nxt;
        phi = GetAngle(ref, r);
    }

    friend std::ostream& operator<< (std::ostream &out, const Vert &v) {
        out << (Point) v
            << ", nxt=" << v.nxt;
        return out;
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
            << ", f: " << (b.f+1)
            << ", g: " << (b.g+1);
        return out;
    }
};

class Vert2 {
public:
    Vert2 (int _i, double _l) : i(_i), l(_l) {}
    double l;
    int i;
};

typedef std::vector<Vert2> VertsType2;

class Pos {
public:
    int edA, edB;
    double t;
    Pos (int edA, int edB, double t) : edA(edA), edB(edB), t(t) {}
    Pos () {}

    friend std::ostream& operator<< (std::ostream &out, const Pos &p) {
        out << "t: " << p.t
            << ", edA: " << p.edA
            << ", edB: " << p.edB;
        return out;
    }

    bool operator== (const Pos &other) const {
        return edA == other.edA && edB == other.edB;
    }
};

class Tracker {
public:
    std::map<int, Pos> locs;
    Tracker (const PolyType &poly) {
        PolyType::const_reverse_iterator itr;
        for (itr = poly.rbegin(); itr != poly.rend(); ++itr) {
            locs[itr->tag] = { itr->tag, (itr+1 != poly.rend() ? itr+1 : poly.rbegin())->tag, 0 };
        }
    }
    void Track (const Point &before, const Point &after, const Point &p, double t) {
        const Pos &posB = locs[before.tag],
            &posA = locs[after.tag];

        double tA = posB == posA ? posA.t : 1;

        locs[p.tag] = { posB.edA, posB.edB, posB.t+(tA-posB.t)*t };
    }

};

class Vert4 : public Point {
public:
    double t;
    Vert4 (const Point &p, double t) : Point(p), t(t) {}

    friend std::ostream& operator<< (std::ostream &out, const Vert4 &v) {
        out << (Point) v
            << ", t: " << v.t;
        return out;
    }
};

typedef std::vector<Vert4> VertsType4;

typedef std::map<Pair, VertsType4> SavedPtsType;
typedef std::shared_ptr<SavedPtsType> SavedPtsPtr;

typedef std::set<int> SpecTagsType;
typedef std::shared_ptr<SpecTagsType> SpecTagsPtr;

void Simplify (const PolyType &poly, SavedPtsPtr &savedPts, SpecTagsPtr &specTags, PolyType &res, int skip, bool rev);

void Align (PolyType &poly, const Point &p);

void Restore (const PolyType &poly, const Tracker &tr, const SavedPtsType &savedPts, PolyType &res);
void Restore2 (const PolyType &poly, PolyType &res);

void SimpleRestore (const PolyType &poly, const SavedPtsType &savedPts, PolyType &res);

// diese darf nicht direkt verwendet werden
void GetVisPoly (PolyType &poly, Tracker &tr, PolyType &res, SavedPtsType &savedPts, int ind = 0);

void GetVisPoly_wrapper (PolyType &poly, PolyType &res, int ind);

#endif
