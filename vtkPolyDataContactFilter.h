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

#ifndef __vtkPolyDataContactFilter_h
#define __vtkPolyDataContactFilter_h

#include <vtkPolyDataAlgorithm.h>

class vtkOBBNode;
class vtkMatrix4x4;

class InterPtType {
public:
    InterPtType () : onEdge(false), end(-1), count(1) {}

    double pt[3];
    double t;
    bool onEdge;

    int edge[2];

    int end, count;

    bool operator< (const InterPtType &other) const {
        return t < other.t;
    }

#ifdef DEBUG
    int ind;
#endif
};

typedef std::vector<InterPtType> InterPtsType;

typedef std::vector<std::pair<InterPtType, InterPtType> > OverlapsType;

class VTK_EXPORT vtkPolyDataContactFilter : public vtkPolyDataAlgorithm {

    void PreparePolyData (vtkPolyData *pd);

    InterPtType InterEdgeLine (double *edgePtA, double *edgePtB, double *r, double *pt);
    InterPtsType InterPolyLine (vtkPoints *pts, vtkIdList *poly, double *r, double *pt);
    void InterPolys (vtkIdType idA, vtkIdType idB);
    OverlapsType OverlapLines (InterPtsType &intersA, InterPtsType &intersB);

    vtkIntArray *contA, *contB;

    vtkPolyData *contLines;
    vtkPoints *contPts;

    vtkPolyData *pdA, *pdB;

public:
    vtkTypeMacro(vtkPolyDataContactFilter, vtkPolyDataAlgorithm);

    static vtkPolyDataContactFilter* New();

    static int InterOBBNodes (vtkOBBNode *nodeA, vtkOBBNode *nodeB, vtkMatrix4x4 *mat, void *caller);

protected:
    vtkPolyDataContactFilter ();
    ~vtkPolyDataContactFilter ();

    int ProcessRequest (vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector);

};

#endif
