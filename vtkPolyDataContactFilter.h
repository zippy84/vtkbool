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
    InterPtType () : ind(0), onEdge(false), end(NO_USE), count(1), srcA(NO_USE), srcB(NO_USE) {}

    int ind;

    double pt[3];
    double t;
    bool onEdge;

    int edge[2];

    int end, count;

    int srcA, srcB;

    bool operator< (const InterPtType &other) const {
        return t < other.t;
    }

    friend std::ostream& operator<< (std::ostream &out, const InterPtType &s) {
        out << "ind " << s.ind
            << ", pt [" << s.pt[0] << ", " << s.pt[1] << ", " << s.pt[2] << "]"
            << ", t " << s.t
            << ", edge [" << s.edge[0] << ", " << s.edge[1] << "]"
            << ", end " << s.end
            << ", src " << s.src;

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

    InterPtType InterEdgeLine (double *edgePtA, double *edgePtB, double *r, double *pt, int _pid);
    InterPtsType InterPolyLine (vtkPoints *pts, vtkIdList *poly, double *r, double *pt, _Src src, int _pid);
    void InterPolys (vtkIdType idA, vtkIdType idB);
    OverlapsType OverlapLines (InterPtsType &intersA, InterPtsType &intersB);

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
