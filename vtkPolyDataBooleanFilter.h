/*
   Copyright 2012-2014 Ronald Römer

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

#ifndef __vtkPolyDataBooleanFilter_h
#define __vtkPolyDataBooleanFilter_h

#include <vector>
#include <deque>
#include <map>
#include <utility>

#include <vtkPolyDataAlgorithm.h>
#include <vtkKdTreePointLocator.h>

#define LOC_NONE 0
#define LOC_INSIDE 1
#define LOC_OUTSIDE 2

#define OPER_UNION 0
#define OPER_INTERSECTION 1
#define OPER_DIFFERENCE 2
#define OPER_DIFFERENCE2 3

#define CAPT_NOT 0
#define CAPT_EDGE 1
#define CAPT_A 2
#define CAPT_B 3

class StripPtType {
public:
    StripPtType () : T(0), t(0), edgeA(-1), edgeB(-1), onEnd(false), capt(CAPT_NOT) {}

    double T, t;
    int ind;
    double pt[3];
    int edgeA, edgeB;
    bool onEnd;
    int capt;
    double captPt[3];

    double cutPt[3];

    bool operator< (const StripPtType &other) const {
        return t < other.t;
    }

#ifdef DEBUG
    // zur überprüfung der sortierung
    int pos;
#endif
};

typedef std::deque<StripPtType> StripType;

typedef std::vector<StripType> StripsType;

typedef std::map<int, StripsType> PolyStripsType;

typedef std::map<int, StripPtType> StripPointsType;

typedef std::vector<std::pair<double, int> > ReplsType;

class MergeVertType {
public:
    MergeVertType () : used(false) {}

    int polyInd;

    // dies ist die ecke die ausgetauscht werden muss
    int vert, vertInd;

    // dies ist der punkt für die paarbildung
    double pt[3];

    bool used;

#ifdef DEBUG
    int ind;
#endif

};


class VTK_EXPORT vtkPolyDataBooleanFilter : public vtkPolyDataAlgorithm {
    vtkPolyData *resultA, *resultB, *contLines;

    vtkPolyData *modPdA, *modPdB;

    vtkCellData *cellDataA, *cellDataB;

    unsigned long timePdA, timePdB;

    std::vector<int> locked;

    StripPointsType GetStripPoints(vtkPolyData *pd, int polyInd, std::deque<int> &lines);
    PolyStripsType GetPolyStrips (vtkPolyData *pd, vtkIntArray *conts);

    StripsType GetAllStrips (PolyStripsType &polyStrips);

    void CutCells (vtkPolyData *pd, PolyStripsType &polyStrips);
    void RemoveDuplicates (vtkPolyData *pd, std::deque<int> &lines);

    void CompleteStrips (StripsType &strips);

    bool HasArea (StripType &strip);

    void RestoreOrigPoints (vtkPolyData *pd, StripsType &strips);
    void RestoreOrigPt (vtkPolyData *pd, vtkKdTreePointLocator *loc, StripPtType &stripPt);

    void ResolveOverlaps (vtkPolyData *pd, StripsType &strips, vtkIntArray *conts);

    void AddAdjacentPoints (vtkPolyData *pd, StripsType &strips);
    void DisjoinPolys (vtkPolyData *pd, StripsType &strips);
    void MergePoints (vtkPolyData *pd, StripsType &strips);

    void CombineRegions ();

    int OperMode;

public:
    vtkTypeMacro(vtkPolyDataBooleanFilter, vtkPolyDataAlgorithm);
    static vtkPolyDataBooleanFilter* New ();

    static bool SortFct (const StripType &stripA, const StripType &stripB, double *n);

    vtkSetClampMacro(OperMode, int, OPER_UNION, OPER_DIFFERENCE2);
    vtkGetMacro(OperMode, int);

    void SetOperModeToUnion () { OperMode = OPER_UNION; };
    void SetOperModeToIntersection () { OperMode = OPER_INTERSECTION; };
    void SetOperModeToDifference () { OperMode = OPER_DIFFERENCE; };
    void SetOperModeToDifference2 () { OperMode = OPER_DIFFERENCE2; };

protected:
    vtkPolyDataBooleanFilter ();
    ~vtkPolyDataBooleanFilter ();

    int ProcessRequest (vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector);

};

#endif
