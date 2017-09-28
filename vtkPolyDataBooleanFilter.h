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

#ifndef __vtkPolyDataBooleanFilter_h
#define __vtkPolyDataBooleanFilter_h

#include <vector>
#include <deque>
#include <map>
#include <set>
#include <utility>
#include <iostream>

#include <vtkPolyDataAlgorithm.h>
#include <vtkKdTreePointLocator.h>

#include "Utilities.h"

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

#define SIDE_START 0
#define SIDE_END 1

class StripPt {
public:
    StripPt () : t(0), capt(CAPT_NOT) {
        edge[0] = NO_USE;
        edge[1] = NO_USE;
    }

    double t;
    int ind;
    double pt[3];
    int edge[2];
    int capt;
    double captPt[3];

    double cutPt[3];

    friend std::ostream& operator<< (std::ostream &out, const StripPt &s) {
        out << "ind=" << s.ind
            << ", edge=[" << s.edge[0] << ", " << s.edge[1] << "]"
            << ", t=" << s.t
            << ", capt=" << s.capt;
        return out;
    }
};

class StripPtR {
public:
    StripPtR (int _ind) : ind(_ind) {
        strip = NO_USE;
        side = NO_USE;
        ref = NO_USE;
    }

    int ind;
    int desc[2];

    // nicht gesetzt bei CAPT_NOT
    int strip;
    int side;
    int ref;

    friend std::ostream& operator<< (std::ostream &out, const StripPtR &s) {
        out << "ind=" << s.ind
            << ", desc=[" << s.desc[0] << ", " << s.desc[1] << "]"
            << ", strip=" << s.strip
            << ", side=" << s.side
            << ", ref=" << s.ref;
        return out;
    }
};

typedef std::map<int, StripPt> StripPtsType;
typedef std::deque<StripPtR> StripType;
typedef std::vector<StripType> StripsType;

class PStrips {
public:
    PStrips () {}
    double n[3];
    IdsType poly;
    StripPtsType pts;
    StripsType strips;
};

typedef std::map<int, PStrips> PolyStripsType;

typedef std::vector<std::reference_wrapper<StripPtR>> RefsType;

class StripPtL {
public:
    StripPtL (const StripPt &sp) : ind(sp.ind) {
        CPY(pt, sp.pt)
        CPY(cutPt, sp.cutPt)
    }

    int ind;
    double pt[3];
    double cutPt[3];

    bool operator< (const StripPtL &other) const {
        return ind < other.ind;
    }
};

class StripPtL2 {
public:
    StripPtL2 (const StripPt &sp) : ind(sp.ind), t(sp.t) {
        CPY(pt, sp.pt)
        edge[0] = sp.edge[0];
        edge[1] = sp.edge[1];
    }

    int ind;
    int edge[2];
    double pt[3];
    double t;

    bool operator< (const StripPtL2 &other) const {
        return ind < other.ind;
    }
};

class StripPtL3 {
public:
    StripPtL3 (const double *_pt, double _t, int _ind = NO_USE) : t(_t), ind(_ind) {
        CPY(pt, _pt)
    }
    int ind;
    double pt[3];
    double t;

    bool operator< (const StripPtL3 &other) const {
        return t < other.t;
    }

    friend std::ostream& operator<< (std::ostream &out, const StripPtL3 &s) {
        out << "ind=" << s.ind
            << ", t=" << s.t
            << ", pt=[" << s.pt[0] << ", " << s.pt[1] << ", " << s.pt[2] << "]";
        return out;
    }
};

class MergePt {
public:
    MergePt (int _polyInd, int _ind, double *_pt) : polyInd(_polyInd), ind(_ind) {
        CPY(pt, _pt)
    }
    int polyInd;
    int ind;
    double pt[3];
};


class VTK_EXPORT vtkPolyDataBooleanFilter : public vtkPolyDataAlgorithm {
    vtkPolyData *resultA, *resultB, *contLines;
    vtkPolyData *modPdA, *modPdB;
    vtkCellData *cellDataA, *cellDataB;

    unsigned long timePdA, timePdB;

    PolyStripsType polyStripsA, polyStripsB;

    std::set<int> involvedA, involvedB;

    void GetStripPoints (vtkPolyData *pd, PStrips &pStrips, IdsType &lines);
    void GetPolyStrips (vtkPolyData *pd, vtkIntArray *conts, PolyStripsType &polyStrips);
    void RemoveDuplicates (IdsType &lines);
    void CompleteStrips (PStrips &pStrips);
    bool HasArea (StripType &strip);
    void CutCells (vtkPolyData *pd, PolyStripsType &polyStrips);
    void RestoreOrigPoints (vtkPolyData *pd, PolyStripsType &polyStrips);
    void DisjoinPolys (vtkPolyData *pd, PolyStripsType &polyStrips);
    void ResolveOverlaps (vtkPolyData *pd, vtkIntArray *conts, PolyStripsType &polyStrips);
    void AddAdjacentPoints (vtkPolyData *pd, vtkIntArray *conts, PolyStripsType &polyStrips);
    void MergePoints (vtkPolyData *pd, PolyStripsType &polyStrips);
    void CombineRegions ();
    void MergeRegions ();

    int OperMode;
    bool MergeAll;

public:
    vtkTypeMacro(vtkPolyDataBooleanFilter, vtkPolyDataAlgorithm);
    static vtkPolyDataBooleanFilter* New ();

    vtkSetClampMacro(OperMode, int, OPER_UNION, OPER_DIFFERENCE2);
    vtkGetMacro(OperMode, int);

    void SetOperModeToUnion () { OperMode = OPER_UNION; };
    void SetOperModeToIntersection () { OperMode = OPER_INTERSECTION; };
    void SetOperModeToDifference () { OperMode = OPER_DIFFERENCE; };
    void SetOperModeToDifference2 () { OperMode = OPER_DIFFERENCE2; };

    vtkSetMacro(MergeAll, bool);
    vtkGetMacro(MergeAll, bool);
    vtkBooleanMacro(MergeAll, bool);

protected:
    vtkPolyDataBooleanFilter ();
    ~vtkPolyDataBooleanFilter ();

    int ProcessRequest (vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector);

};

#endif
