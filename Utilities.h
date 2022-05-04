/*
Copyright 2012-2022 Ronald Römer

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

#ifndef __Utilities_h
#define __Utilities_h

#include <iostream>
#include <type_traits>

#include <vtkPolyData.h>
#include <vtkKdTreePointLocator.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkMath.h>

#define NOTSET -1

double GetAngle (double *vA, double *vB, double *n);

/* VTK */
void ComputeNormal (vtkPoints *pts, double *n, vtkIdType num, const vtkIdType *poly);
void FindPoints (vtkKdTreePointLocator *pl, const double *pt, vtkIdList *pts, double tol = 1e-6);
void WriteVTK (const char *name, vtkPolyData *pd);

class Point3d {
public:
    const double x, y, z;
    vtkIdType id;

    Point3d () = delete;
    Point3d (const double _x, const double _y, const double _z, vtkIdType _id = NOTSET) : x(_x), y(_y), z(_z), id(_id) {}
    bool operator< (const Point3d &other) const {
        const long x1 = std::lround(x*1e5),
            y1 = std::lround(y*1e5),
            z1 = std::lround(z*1e5),
            x2 = std::lround(other.x*1e5),
            y2 = std::lround(other.y*1e5),
            z2 = std::lround(other.z*1e5);

        return std::tie(x1, y1, z1) < std::tie(x2, y2, z2);
    }
    bool operator== (const Point3d &other) const {
        return !(*this < other) && !(other < *this);
    }
};

typedef std::vector<vtkIdType> IdsType;

template<typename T,
    typename std::enable_if<std::is_integral<T>::value, bool>::type = true>
class _Pair {
public:
    T f, g;
    _Pair () = delete;
    _Pair (T _f, T _g) : f(_f), g(_g) {}
    bool operator< (const _Pair &other) const {
        return std::tie(f, g) < std::tie(other.f, other.g);
    }
    bool operator== (const _Pair &other) const {
        return f == other.f && g == other.g;
    }
    friend std::ostream& operator<< (std::ostream &out, const _Pair &p) {
        out << "(" << p.f << ", " << p.g << ")";
        return out;
    }
};

// typedef _Pair<vtkIdType> Pair;
using Pair = _Pair<vtkIdType>;

template<typename T,
    typename U,
    typename std::enable_if<std::is_integral<T>::value
        && std::is_integral<U>::value
        && std::is_signed<T>::value, bool>::type = true>
T Mod (T a, U b) {
    T _b = static_cast<T>(b);

    return ((a%_b)+_b)%_b;
}

class Base {
public:
    Base () {}
    Base (vtkPoints *pts, vtkIdType num, const vtkIdType *poly);
    double n[3], ei[3], ej[3], d;
};

void Transform (const double *in, double *out, const Base &base);
void BackTransform (const double *in, double *out, const Base &base);

inline void Cpy (double *a, const double *b, const int n = 2) {
    std::copy_n(b, n, a);
}

template<typename T>
std::ostream& operator<< (typename std::enable_if<std::is_enum<T>::value, std::ostream>::type& stream, const T& e) {
    return stream << static_cast<typename std::underlying_type<T>::type>(e);
}

typedef std::vector<Point3d> Poly;
typedef std::vector<Poly> PolysType;

#endif
