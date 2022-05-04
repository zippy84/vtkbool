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

#ifndef __vtkPolyDataContactFilter_h
#define __vtkPolyDataContactFilter_h

#include <map>
#include <tuple>
#include <set>

#include <vtkPolyDataAlgorithm.h>
#include <vtkIdTypeArray.h>

#include "Utilities.h"

class vtkOBBNode;
class vtkMatrix4x4;

enum class Src {
    A = 1,
    B = 2
};

class InterPt {
public:
    InterPt () : onEdge(false), end(NOTSET), srcA(NOTSET), srcB(NOTSET) {}

    InterPt (double _t, vtkIdType _end, double x, double y, double z) : InterPt() {
        t = _t;
        onEdge = true;
        end = _end;
        pt[0] = x;
        pt[1] = y;
        pt[2] = z;
    }

    InterPt (double _t, vtkIdType _end, double x, double y, double z, vtkIdType a, vtkIdType b, Src _src) : InterPt(_t, _end, x, y, z) {
        assert(_end == a || _end == b);

        edge[0] = a;
        edge[1] = b;
        src = _src;
    }

    double t;
    bool onEdge;

    vtkIdType edge[2], end, srcA, srcB;

    double pt[3];

    friend std::ostream& operator<< (std::ostream &out, const InterPt &s) {
        out << "pt [" << s.pt[0] << ", " << s.pt[1] << ", " << s.pt[2] << "]"
            << ", t " << s.t
            << ", edge [" << s.edge[0] << ", " << s.edge[1] << "]"
            << ", end " << s.end
            << ", src " << s.src;

        return out;
    }

    void Merge (const InterPt &other) {
        assert(src != other.src);

        if (src == Src::A) {
            srcA = end == NOTSET ? edge[0] : end;
        } else {
            srcB = end == NOTSET ? edge[0] : end;
        }

        if (std::abs(other.t-t) < 1e-5) {
            if (other.src == Src::A) {
                srcA = other.end == NOTSET ? other.edge[0] : other.end;
            } else {
                srcB = other.end == NOTSET ? other.edge[0] : other.end;
            }
        }
    }

    Src src;
};

typedef std::vector<InterPt> InterPtsType;

typedef std::vector<std::tuple<InterPt, InterPt, vtkIdType, vtkIdType>> OverlapsType;

typedef std::set<Pair> SpecialEdgesType;

class VTK_EXPORT vtkPolyDataContactFilter : public vtkPolyDataAlgorithm {

    void PreparePolyData (vtkPolyData *pd, SpecialEdgesType &edges);

    static void InterEdgeLine (InterPtsType &interPts, const double *eA, const double *eB, const double *r, const double *pt);
    static void InterPolyLine (InterPtsType &interPts, vtkPolyData *pd, vtkIdType num, const vtkIdType *poly, const double *r, const double *pt, Src src, const double *n);
    void InterPolys (vtkIdType idA, vtkIdType idB);
    void OverlapLines (OverlapsType &ols, InterPtsType &intersA, InterPtsType &intersB, const vtkIdType *polyA, const vtkIdType *polyB, vtkIdType idA, vtkIdType idB);
    void AddContactLines (InterPtsType &intersA, InterPtsType &intersB, const vtkIdType *polyA, const vtkIdType *polyB, vtkIdType idA, vtkIdType idB);

    vtkIdTypeArray *contA, *contB;

    vtkPolyData *contLines;
    vtkPoints *contPts;

    vtkPolyData *pdA, *pdB;

    vtkIdTypeArray *sourcesA, *sourcesB;

    vtkIdTypeArray *neigsA, *neigsB;

    SpecialEdgesType edgesA, edgesB;

    bool invalidA, invalidB;

public:
    vtkTypeMacro(vtkPolyDataContactFilter, vtkPolyDataAlgorithm);

    static vtkPolyDataContactFilter* New();

    static int InterOBBNodes (vtkOBBNode *nodeA, vtkOBBNode *nodeB, vtkMatrix4x4 *mat, void *caller);

protected:
    vtkPolyDataContactFilter ();
    ~vtkPolyDataContactFilter ();

    int ProcessRequest (vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector);

    void PrintSelf (ostream&, vtkIndent) override {};

private:
    vtkPolyDataContactFilter (const vtkPolyDataContactFilter&) = delete;
    void operator= (const vtkPolyDataContactFilter&) = delete;

};

#endif
