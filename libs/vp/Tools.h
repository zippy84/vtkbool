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

#ifndef __Tools_h
#define __Tools_h

#include <vector>
#include <cmath>
#include <memory>
#include <iostream>
#include <string>
#include <sstream>
#include <tuple>
#include <algorithm>
#include <cfloat>
#include <exception>

#define NO_USE -1
#define PI std::acos(-1)

#define E 1e-5

inline void Cpy (double *a, const double *b, const int n = 2) {
    std::copy_n(b, n, a);
}

class Pair {
public:
    int f, g;
    Pair () {}
    Pair (int _f, int _g) : f(_f), g(_g) {}
    bool operator< (const Pair &other) const {
        return std::tie(f, g) < std::tie(other.f, other.g);
    }
    bool operator== (const Pair &other) const {
        return f == other.f && g == other.g;
    }
    friend std::ostream& operator<< (std::ostream &out, const Pair &p) {
        out << "(" << p.f << ", " << p.g << ")";
        return out;
    }
};

typedef std::vector<int> IdsType;

class Point {
    static int _tag;
public:
    int tag;

    Point (double _x, double _y, int _id = NO_USE) : id(_id), tag(_tag++) {
        pt[0] = _x;
        pt[1] = _y;
    }
    Point (double *_pt, int _id = NO_USE) : Point(_pt[0], _pt[1], _id) {}

    Point (const Point& p) : id(p.id), tag(p.tag) {
        pt[0] = p.pt[0];
        pt[1] = p.pt[1];
    }

    Point& operator= (const Point &p) {
        pt[0] = p.pt[0];
        pt[1] = p.pt[1];
        id = p.id;
        tag = p.tag;

        return *this;
    }

    double pt[2];

    const double &x = pt[0],
        &y = pt[1];
    int id;

    friend std::ostream& operator<< (std::ostream &out, const Point &p) {
        out << "id: " << p.id << ", pt: [" << p.x << "," << p.y << "], tag: " << p.tag;
        return out;
    }

    bool operator< (const Point &other) const {
        const int x1 = static_cast<int>(x*1e5),
            y1 = static_cast<int>(y*1e5),
            x2 = static_cast<int>(other.x*1e5),
            y2 = static_cast<int>(other.y*1e5);

        return std::tie(x1, y1) < std::tie(x2, y2);
    }

};

typedef std::vector<Point> PolyType;

double Normalize (double *v, const int n = 2);
double GetAngle (double *vA, double *vB);
bool Ld (double *a, double *b, double *c, double *r = nullptr);
double Cross (double *a, double *b, double *c);

class D {
public:
    D (double *_s) {
        s[0] = _s[0];
        s[1] = _s[1];
    }
    D (double *_s, double _t1, double _t2) : D(_s) {
        t1 = _t1;
        t2 = _t2;
    }
    double s[2], t1, t2;

    friend std::ostream& operator<< (std::ostream &out, const D &d) {
        out << "t1: " << d.t1
            << ", t2: " << d.t2;
        return out;
    }
};

class Bnds {
public:
    Bnds () : lA(-E), uA(1+E), lB(E), uB(1+E) {}
    Bnds (double lTa, double uTa, double lTb, double uTb) : lA(lTa), uA(1+uTa), lB(lTb), uB(1+uTb) {}
    const double lA, uA, lB, uB;
};

std::shared_ptr<D> Intersect (const double *o, const double *r, const double *pA, const double *pB);
std::shared_ptr<D> Intersect2 (const double *oA, const double *oB, const double *pA, const double *pB, const Bnds &bnds);

std::shared_ptr<D> Intersect2 (const double *oA, const double *oB, const double *pA, const double *pB);

bool IsFrontfaced (double *r, double *a, double *b);
bool IsNear (const double *a, const double *b);
double GetT (double *a, double *b, double *c);

bool TestCW (const PolyType &poly);

std::string GetAbsolutePath (const PolyType &poly);

bool TestPIP (PolyType &poly, Point &pt);

template<typename T>
std::ostream& operator<< (typename std::enable_if<std::is_enum<T>::value, std::ostream>::type& stream, const T& e) {
    return stream << static_cast<typename std::underlying_type<T>::type>(e);
}

class Ext {
public:
    double minX, maxX, minY, maxY;
    Ext () : minX(DBL_MAX), maxX(DBL_MIN), minY(DBL_MAX), maxY(DBL_MIN) {}

    friend std::ostream& operator<< (std::ostream &out, const Ext &e) {
        out << "(" << e.minX
            << ", " << e.maxX
            << ", " << e.minY
            << ", " << e.maxY << ")";
        return out;
    }

    double GetDiag () {
        return (maxX-minX)*(maxX-minX)+(maxY-minY)*(maxY-minY);
    }
};

void GetExt (const PolyType &poly, Ext &ext);

double GetArea (const PolyType &poly);

double GetDis (const Point &pA, const Point &pB, const Point &pC, double &t, double *pro = nullptr);

double GetSqDis (const Point &a, const Point &b);

void GetSect (int tagA, int tagB, PolyType &poly);

inline void vtkbool_throw (bool test, const std::string &where, const std::string &msg) {
    if (!test) {
        throw std::runtime_error("Exception (" + where + ", " + msg + ")");
    }
}

#endif
