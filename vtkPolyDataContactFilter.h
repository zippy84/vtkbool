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

#ifndef __vtkPolyDataContactFilter_h
#define __vtkPolyDataContactFilter_h

#include <map>
#include <tuple>
#include <set>
#include <string>

#include <vtkPolyDataAlgorithm.h>
#include <vtkIdTypeArray.h>

#include "Utilities.h"

class vtkOBBNode;
class vtkMatrix4x4;

enum class Src {
    A,
    B
};

enum class End {
    NONE,
    BEGIN,
    END
};

enum class PointSrc {
    CALCULATED,
    COPIED
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
            srcA = getEnd();
        } else {
            srcB = getEnd();
        }

        if (std::abs(other.t-t) < 1e-5) {
            if (other.src == Src::A) {
                srcA = other.getEnd();
            } else {
                srcB = other.getEnd();
            }
        }
    }

    inline vtkIdType getEnd() const {
        if (end == End::BEGIN) {
            return edge.f;
        }

        if (end == End::END) {
            return edge.g;
        }

        return NOTSET;
    }

};

typedef std::vector<InterPt> InterPtsType;

typedef std::vector<std::tuple<InterPt, InterPt>> OverlapsType;

typedef std::set<Pair> InvalidEdgesType;

class VTK_EXPORT vtkPolyDataContactFilter : public vtkPolyDataAlgorithm {

    void CopyPolyData (vtkPolyData *pd, vtkPolyData *newPd);

    static void InterEdgeLine (InterPtsType &interPts, vtkPolyData *pd, vtkIdType idA, vtkIdType idB, const double *r, const double *pt, Src src);
    static bool InterPolyLine (InterPtsType &interPts, vtkPolyData *pd, vtkIdType num, const vtkIdType *poly, const double *r, const double *pt, Src src, const double *n);
    void InterPolys (vtkIdType idA, vtkIdType idB);
    void OverlapLines (OverlapsType &ols, InterPtsType &intersA, InterPtsType &intersB, vtkIdType idA, vtkIdType idB);
    void AddContactLines (InterPtsType &intersA, InterPtsType &intersB, vtkIdType idA, vtkIdType idB);

    static bool CheckInters (InterPtsType &interPts, vtkPolyData *pd);

    vtkSmartPointer<vtkPolyData> newPdA, newPdB;

    vtkSmartPointer<vtkPolyData> contLines;
    vtkSmartPointer<vtkIdTypeArray> contA, contB, sourcesA, sourcesB;

    bool invalidA, invalidB;
    InvalidEdgesType edgesA, edgesB;

    void GetInvalidEdges (vtkPolyData *pd, InvalidEdgesType &edges);

    bool aborted;

public:
    vtkTypeMacro(vtkPolyDataContactFilter, vtkPolyDataAlgorithm);

    static vtkPolyDataContactFilter* New();

    static int InterOBBNodes (vtkOBBNode *nodeA, vtkOBBNode *nodeB, vtkMatrix4x4 *mat, void *caller);

protected:
    vtkPolyDataContactFilter ();
    ~vtkPolyDataContactFilter ();

    int RequestData (vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;

    void PrintSelf (ostream&, vtkIndent) override {};

private:
    vtkPolyDataContactFilter (const vtkPolyDataContactFilter&) = delete;
    void operator= (const vtkPolyDataContactFilter&) = delete;

};


using Edges = std::map<Pair, Points>;

class SnapPoint {
public:
    SnapPoint () = delete;
    SnapPoint (vtkIdType cellId, const Pair &line, const Point3d &point, const Point3d &inter, double d) : cellId(cellId), line(line), point(point), inter(inter), d(d) {}

    vtkIdType cellId;
    Pair line;
    Point3d point;
    Point3d inter;
    double d;

    friend std::ostream& operator<< (std::ostream &out, const SnapPoint &s) {
        out << "cellId " << s.cellId
            << ", line " << s.line
            << ", point " << s.point
            << ", inter " << s.inter
            << ", d " << s.d;

        return out;
    }
};

class SnapEdge {
public:
    SnapEdge () = delete;
    SnapEdge (vtkIdType cellId, const Pair &line, const Pair &edge, const Point3d &proj, const Point3d &inter, double d) : cellId(cellId), line(line), edge(edge), proj(proj), inter(inter), d(d) {}

    vtkIdType cellId;
    Pair line;
    Pair edge;
    Point3d proj;
    Point3d inter;
    double d;

    friend std::ostream& operator<< (std::ostream &out, const SnapEdge &s) {
        out << "cellId " << s.cellId
            << ", line " << s.line
            << ", edge " << s.edge
            << ", proj " << s.proj
            << ", inter " << s.inter
            << ", d " << s.d;

        return out;
    }
};

class PreventEqualCaptPoints {
    vtkPolyData *pdA, *pdB;
public:
    static void TriangluteCell (vtkPolyData *pd, vtkIdType cellId, const Edges &edges);
    static void MovePoint (vtkPolyData *pd, vtkIdType ind, const Point3d &p);

    PreventEqualCaptPoints (vtkPolyData *pdA, vtkPolyData *pdB);
    void Run ();
private:
    void Find (vtkPolyData *pd, vtkPolyData *other, const std::string &name);
};

#endif
