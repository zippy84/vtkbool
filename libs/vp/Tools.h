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

double Normalize (double *v, const int n = 2);
double GetAngle (double *vA, double *vB);
void Move (double *a, double *b, double *c);
double Ld (double *a, double *b, double *c);
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
};

std::shared_ptr<D> Intersect (double *o, double *r, double *pA, double *pB);
std::shared_ptr<D> Intersect2 (double *oA, double *oB, double *pA, double *pB);

bool IsFrontfaced (double *r, double *a, double *b);
bool IsNear (double *a, double *b);
double GetT (double *a, double *b, double *c);

bool IsOnSeg (double *a, double *b, double *c);

bool TestCW (const PolyType &poly);

std::string GetAbsolutePath (const PolyType &poly);

bool TestPIP (PolyType &poly, Point &pt);

template<typename T>
std::ostream& operator<< (typename std::enable_if<std::is_enum<T>::value, std::ostream>::type& stream, const T& e) {
    return stream << static_cast<typename std::underlying_type<T>::type>(e);
}

#endif
