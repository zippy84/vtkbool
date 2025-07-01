/*
Copyright 2012-2025 Ronald Römer

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

#ifndef __Contact_h
#define __Contact_h

#include "Utilities.h"

#include <vtkOBBTree.h>
#include <vtkMatrix4x4.h>

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

class InterPt {
public:
    InterPt () = delete;

    InterPt (double x, double y, double z, double t, vtkIdType a, vtkIdType b, End end, Src src, PointSrc pointSrc) : t(t), edge(a, b), end(end), src(src), srcA(NOTSET), srcB(NOTSET), pointSrc(pointSrc) {
        pt[0] = x;
        pt[1] = y;
        pt[2] = z;
    }

    double pt[3], t;
    Pair edge;
    End end;
    Src src;
    vtkIdType srcA, srcB;

    PointSrc pointSrc;

    friend std::ostream& operator<< (std::ostream &out, const InterPt &s) {
        out << "pt [" << s.pt[0] << ", " << s.pt[1] << ", " << s.pt[2] << "]"
            << ", t " << s.t
            << ", edge " << s.edge
            << ", end " << s.end
            << ", src " << s.src
            << ", pointSrc " << s.pointSrc;

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

    inline vtkIdType GetEnd () const {
        if (end == End::A) {
            return edge.f;
        }

        if (end == End::B) {
            return edge.g;
        }

        return NOTSET;
    }

};

typedef std::vector<InterPt> InterPtsType;
typedef std::vector<std::tuple<InterPt, InterPt>> OverlapsType;

typedef std::set<Pair> NonManifoldEdgesType;

vtkSmartPointer<vtkPolyData> Clean (vtkPolyData *pd);

class Contact {
public:
    Contact () = delete;
    Contact (vtkPolyData *newPdA, vtkPolyData *newPdB);

    vtkPolyData *newPdA, *newPdB;

    vtkSmartPointer<vtkPoints> pts;
    vtkSmartPointer<vtkPolyData> lines;
    vtkSmartPointer<vtkIdTypeArray> contA, contB, sourcesA, sourcesB;

    bool touchesEdgesA, touchesEdgesB;

    NonManifoldEdgesType edgesA, edgesB;

    vtkSmartPointer<vtkPolyData> GetLines ();

    void GetNonManifoldEdges (vtkPolyData *pd, NonManifoldEdgesType &edges);

    void InterEdgeLine (InterPtsType &interPts, const Point3d &pA, const Point3d &pB, Src src);

    bool InterPolyLine (InterPtsType &interPts, const Base2 &base, const Poly &poly, Src src);

    void InterPolys (vtkIdType idA, vtkIdType idB);

    bool CheckInters (const InterPtsType &interPts, vtkPolyData *pd);

    void OverlapLines (OverlapsType &overlaps, InterPtsType &intersA, InterPtsType &intersB);

    void AddContactLines (InterPtsType &intersA, InterPtsType &intersB, vtkIdType idA, vtkIdType idB);

    static int InterNodes (vtkOBBNode *nodeA, vtkOBBNode *nodeB, vtkMatrix4x4 *vtkNotUsed(matrix), void *ptr);

    std::vector<Pair> pairs;

    std::map<vtkIdType, IdsType> replsA, replsB;

    void IntersectReplacements ();
};

#endif
