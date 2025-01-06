/*
Copyright 2012-2024 Ronald Römer

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
#include <functional>
#include <vector>
#include <deque>
#include <map>
#include <set>
#include <memory>

#include <vtkPolyData.h>
#include <vtkKdTreePointLocator.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkMath.h>

#define NOTSET -1

double GetAngle (const double *vA, const double *vB, const double *n);

/* VTK */
double ComputeNormal (vtkPoints *pts, double *n, vtkIdType num, const vtkIdType *poly);
bool CheckNormal (vtkPoints *pts, vtkIdType num, const vtkIdType *poly, const double *n, double d);

void FindPoints (vtkKdTreePointLocator *pl, const double *pt, vtkIdList *pts, double tol = 1e-6);
void WriteVTK (const char *name, vtkPolyData *pd);

class Point3d {
public:
    const double x, y, z;
    vtkIdType id;

    Point3d () = delete;
    Point3d (const double x, const double y, const double z, vtkIdType id = NOTSET) : x(x), y(y), z(z), id(id) {}
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

    friend std::ostream& operator<< (std::ostream &out, const Point3d &p) {
        out << "Point3d(x=" << p.x
            << ", y=" << p.y
            << ", z=" << p.z
            << ", id=" << p.id << ")";
        return out;
    }

    static double GetVec (const Point3d &a, const Point3d &b, double *v) {
        v[0] = b.x-a.x;
        v[1] = b.y-a.y;
        v[2] = b.z-a.z;

        return vtkMath::Normalize(v);
    }

    static double GetDist (const Point3d &a, const Point3d &b) {
        double dx = b.x-a.x;
        double dy = b.y-a.y;
        double dz = b.z-a.z;

        return dx*dx+dy*dy+dz*dz;
    }
};

typedef std::vector<vtkIdType> IdsType;
typedef std::set<vtkIdType> _IdsType;

class Pair {
public:
    vtkIdType f, g;
    Pair () = delete;
    Pair (vtkIdType f, vtkIdType g) : f(f), g(g) {}
    bool operator< (const Pair &other) const {
        return std::tie(f, g) < std::tie(other.f, other.g);
    }
    bool operator== (const Pair &other) const {
        return f == other.f && g == other.g;
    }
    vtkIdType operator& (const Pair &other) const {
        if (f == other.f || f == other.g) {
            return f;
        }
        return g;
    }
    friend std::ostream& operator<< (std::ostream &out, const Pair &p) {
        out << "(" << p.f << ", " << p.g << ")";
        return out;
    }
};

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

typedef std::vector<Point3d> Poly, Points;
typedef std::vector<Poly> PolysType;

double ComputeNormal (const Poly &poly, double *n);
bool PointInPoly (const Poly &poly, const Point3d &p);

void WritePolys (const char *name, const PolysType &polys);

typedef std::deque<vtkIdType> IndexedPoly;
typedef std::vector<IndexedPoly> IndexedPolysType;

typedef std::map<vtkIdType, std::reference_wrapper<const Point3d>> ReferencedPointsType;

void GetPolys (const ReferencedPointsType &pts, const IndexedPolysType &indexedPolys, PolysType &polys);

void GetPoly (vtkPoints *pts, vtkIdType num, const vtkIdType *poly, Poly &out);
void FlattenPoly (const Poly &poly, Poly &out, const Base &base);

class Proj {
public:
    Proj (vtkIdType a, vtkIdType b, const Point3d &proj, double d) : edge(a, b), proj(proj), d(d) {}

    Pair edge;
    Point3d proj;
    double d;
};

std::shared_ptr<Proj> GetEdgeProj (const Poly &poly, const Point3d &p);

void ProjOnLine (const Point3d &a, const Point3d &b, const Point3d &p, double *d, double *t, std::shared_ptr<Point3d> &proj);
void ProjOnLine (vtkPolyData *pd, const Pair &line, const Point3d &p, std::shared_ptr<Point3d> &proj);

#endif
