/*
Copyright 2012-2020 Ronald Römer

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

#include <vtkPolyDataAlgorithm.h>

#include "Utilities.h"

class vtkOBBNode;
class vtkMatrix4x4;

enum class _Src {
    A = 1,
    B = 2
};

class InterPtType {
public:
    InterPtType () : onEdge(false), end(NO_USE), srcA(NO_USE), srcB(NO_USE), count(1) {}

    double pt[3];
    double t;
    bool onEdge;

    vtkIdType edge[2], end, srcA, srcB;

    int count;

    bool operator< (const InterPtType &other) const {
        return t < other.t;
    }

    friend std::ostream& operator<< (std::ostream &out, const InterPtType &s) {
        out << "pt [" << s.pt[0] << ", " << s.pt[1] << ", " << s.pt[2] << "]"
            << ", t " << s.t
            << ", edge [" << s.edge[0] << ", " << s.edge[1] << "]"
            << ", end " << s.end
            << ", src " << s.src
            << ", count " << s.count;

        return out;
    }

    void Merge (const InterPtType &other) {
        assert(src != other.src);

        if (src == _Src::A) {
            srcA = end == NO_USE ? edge[0] : end;
        } else {
            srcB = end == NO_USE ? edge[0] : end;
        }

        if (std::abs(other.t-t) < 1e-5) {
            if (other.src == _Src::A) {
                srcA = other.end == NO_USE ? other.edge[0] : other.end;
            } else {
                srcB = other.end == NO_USE ? other.edge[0] : other.end;
            }
        }
    }

    _Src src;

};

typedef std::vector<InterPtType> InterPtsType;

typedef std::vector<std::pair<InterPtType, InterPtType> > OverlapsType;

class VTK_EXPORT vtkPolyDataContactFilter : public vtkPolyDataAlgorithm {

    void PreparePolyData (vtkPolyData *pd);

    static void InterEdgeLine (InterPtType &inter, const double *eA, const double *eB, const double *r, const double *pt, vtkIdType _pid);
    static void InterPolyLine (InterPtsType &interPts, vtkPolyData *pd, vtkIdType num, const vtkIdType *poly, const double *r, const double *pt, _Src src, double _pt[3], double _n[3], double _d, vtkIdType _pid);
    void InterPolys (vtkIdType idA, vtkIdType idB);
    static void OverlapLines (OverlapsType &ols, InterPtsType &intersA, InterPtsType &intersB);

    vtkIntArray *contA, *contB;

    vtkPolyData *contLines;
    vtkPoints *contPts;

    vtkPolyData *pdA, *pdB;

    vtkIntArray *sourcesA, *sourcesB;

public:
    vtkTypeMacro(vtkPolyDataContactFilter, vtkPolyDataAlgorithm);

    static vtkPolyDataContactFilter* New();

    static int InterOBBNodes (vtkOBBNode *nodeA, vtkOBBNode *nodeB, vtkMatrix4x4 *mat, void *caller);

protected:
    vtkPolyDataContactFilter ();
    ~vtkPolyDataContactFilter ();

    int ProcessRequest (vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector);

private:
    vtkPolyDataContactFilter (const vtkPolyDataContactFilter&) = delete;
    void operator= (const vtkPolyDataContactFilter&) = delete;

};

#endif
