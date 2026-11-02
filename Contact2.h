/*
Copyright 2012-2026 Ronald Römer

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

#ifndef __Contact2_h
#define __Contact2_h

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include <optional>
#include <vector>
#include <set>
#include <iostream>

struct Vector3 {
    double x, y, z;
};

struct Face {
    std::vector<vtkIdType> vertIds;
    Vector3 normal{};
    Vector3 centroid{};
    double d = 0;
    std::optional<size_t> globalPlaneId;
    vtkIdType origCellId;
};

struct Mesh {
    std::vector<Vector3> verts;
    std::vector<Face> faces;
    std::vector<std::vector<vtkIdType>> pointToFaces;
};

struct GlobalPlane {
    Vector3 normal;
    double d;
    std::vector<size_t> faceIdsA, faceIdsB;
};

struct Pair2 {
    vtkIdType a, b;
};

enum class Src {
    A,
    B
};

enum class End {
    None,
    A,
    B
};

enum class PointSrc {
    Calculated,
    Copied
};

/*template<typename T>
std::ostream& operator<< (typename std::enable_if<std::is_enum<T>::value, std::ostream>::type& stream, const T& e) {
    return stream << static_cast<typename std::underlying_type<T>::type>(e);
}*/

class InterPt {
public:
    InterPt () = delete;

    InterPt (double x, double y, double z, double t, vtkIdType a, vtkIdType b, End end, Src src, PointSrc pointSrc = PointSrc::Calculated) : t(t), edge{a, b}, end(end), src(src), srcA(std::nullopt), srcB(std::nullopt), pointSrc(pointSrc) {
        pt[0] = x;
        pt[1] = y;
        pt[2] = z;
    }

    double pt[3], t;
    Pair2 edge;
    End end;
    Src src;
    std::optional<vtkIdType> srcA, srcB;
    PointSrc pointSrc;

    friend std::ostream& operator<< (std::ostream &out, const InterPt &s) {
        out << "pt [" << s.pt[0] << ", " << s.pt[1] << ", " << s.pt[2] << "]"
            << ", t " << s.t
            << ", edge (" << s.edge.a << ", " << s.edge.b << ")"
            /*<< ", end " << s.end
            << ", src " << s.src
            << ", pointSrc " << s.pointSrc*/;

        return out;
    }

    void Merge (const InterPt &other) {
        assert(src != other.src);

        if (src == Src::A) {
            srcA = GetEnd();
        } else {
            srcB = GetEnd();
        }

        if (std::abs(other.t-t) < 1e-5) {
            if (other.src == Src::A) {
                srcA = other.GetEnd();
            } else {
                srcB = other.GetEnd();
            }
        }
    }

    inline std::optional<vtkIdType> GetEnd () const {
        if (end == End::A) {
            return edge.a;
        }

        if (end == End::B) {
            return edge.b;
        }

        return std::nullopt;
    }

};

typedef std::vector<InterPt> InterPts;

typedef std::set<std::pair<vtkIdType, vtkIdType>> NonManifoldEdges;

class Contact {
    Mesh meshA, meshB;
    vtkSmartPointer<vtkPolyData> pdA, pdB, lines;

    double tol;

    vtkSmartPointer<vtkPolyData> CreateLines ();
public:
    Contact (vtkPolyData *origPdA, vtkPolyData *origPdB);
    vtkSmartPointer<vtkPolyData> GetPolyDataA ();
    vtkSmartPointer<vtkPolyData> GetPolyDataB ();
    vtkSmartPointer<vtkPolyData> GetLines ();

    void Save (vtkPolyData *pd, const std::string &filename);
    void SaveAll (const std::string prefix = "");
};

#endif
