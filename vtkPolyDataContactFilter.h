/*
   Copyright 2012, 2013 Ronald Römer

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
    InterPtType () : isOnEdge(false), onEnd(-1) {}

    bool operator< (const InterPtType &otherPt) const {
        return t < otherPt.t;
    }

    double t;
    bool isOnEdge;
    int onEnd;
    double pt[3];
};

typedef std::vector<InterPtType> InterPtsType;

typedef std::vector<std::pair<InterPtType, InterPtType> > OverlapsType;

class VTK_EXPORT vtkPolyDataContactFilter : public vtkPolyDataAlgorithm {

    InterPtType IntersectEdgeAndLine (double *edgePtA, double *edgePtB, double *r, double *pt);
    InterPtsType IntersectPolyAndLine (vtkPoints *pts, vtkIdList *poly, double *r, double *pt);
    void IntersectPolys (vtkIdType idA, vtkIdType idB);
    OverlapsType OverlapLines (InterPtsType &intersA, InterPtsType &intersB);

    vtkIntArray *contA, *contB;

    vtkPolyData *contLines;
    vtkPoints *contPts;

    vtkPolyData *pdA, *pdB;

    bool MergeLines;

public:
    vtkTypeMacro(vtkPolyDataContactFilter, vtkPolyDataAlgorithm);

    static vtkPolyDataContactFilter* New();

    static int IntersectOBBNodes (vtkOBBNode *nodeA, vtkOBBNode *nodeB, vtkMatrix4x4 *mat, void *caller);

    vtkSetMacro(MergeLines, bool);
    vtkGetMacro(MergeLines, bool);
    vtkBooleanMacro(MergeLines, bool);

protected:
    vtkPolyDataContactFilter ();
    ~vtkPolyDataContactFilter ();

    int ProcessRequest (vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector);

};

#endif
