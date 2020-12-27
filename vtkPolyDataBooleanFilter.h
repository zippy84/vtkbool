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

#define OPER_UNION 0
#define OPER_INTERSECTION 1
#define OPER_DIFFERENCE 2
#define OPER_DIFFERENCE2 3

enum class Capt {
    NOT,
    EDGE,
    A,
    B
};

enum class Side {
    NONE,
    START,
    END
};

enum class Loc {
    NONE,
    INSIDE,
    OUTSIDE
};

class StripPt {
public:
    StripPt () : t(0), capt(Capt::NOT), catched(true) {
        edge[0] = NO_USE;
        edge[1] = NO_USE;
    }

    double t;
    Capt capt;
    double captPt[3];

    vtkIdType ind, edge[2];

    double pt[3];
    double cutPt[3];

    friend std::ostream& operator<< (std::ostream &out, const StripPt &s) {
        out << "ind " << s.ind
            << ", edge [" << s.edge[0] << ", " << s.edge[1] << "]"
            << ", t " << s.t
            << ", capt " << s.capt
            << ", polyId " << s.polyId;
        return out;
    }

    std::vector<Pair> history;

    vtkIdType polyId;

    vtkIdType src;
    bool catched;
};

class StripPtR {
public:
    StripPtR (vtkIdType _ind) : ind(_ind), strip(NO_USE), ref(NO_USE), side(Side::NONE) {
        desc[0] = NO_USE;
        desc[1] = NO_USE;
    }

    vtkIdType ind, desc[2];

    // nicht gesetzt bei CAPT_NOT
    vtkIdType strip, ref;

    Side side;

    friend std::ostream& operator<< (std::ostream &out, const StripPtR &s) {
        out << "ind " << s.ind
            << ", desc [" << s.desc[0] << ", " << s.desc[1] << "]"
            << ", strip " << s.strip
            << ", side " << s.side
            << ", ref " << s.ref;
        return out;
    }
};

typedef std::map<vtkIdType, StripPt> StripPtsType;
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

typedef std::map<vtkIdType, PStrips> PolyStripsType;

typedef std::vector<std::reference_wrapper<StripPtR>> RefsType;

class MergePt {
public:
    MergePt (vtkIdType _polyId, vtkIdType _ind) : polyId(_polyId), ind(_ind) {}
    vtkIdType polyId, ind;
    bool operator== (const MergePt &other) const {
        return polyId == other.polyId && ind == other.ind;
    };

    friend std::ostream& operator<< (std::ostream &out, const MergePt &p) {
        out << "ind " << p.ind
            << ", polyId " << p.polyId;
        return out;
    }
};

// typedef std::vector<IdsType> HolesType;

// class _Wrapper {
//     vtkPolyData *pd;
//     IdsType descIds;
//     vtkIdType origId;

//     Base base;
//     HolesType holes;
// public:
//     _Wrapper (vtkPolyData* _pd, IdsType& _descIds, vtkIdType _origId)
//         : pd(_pd), descIds(_descIds), origId(_origId) {}

//     void MergeAll ();
//     void Add (IdsType &hole) {
//         holes.push_back(hole);
//     }
// };

typedef std::set<vtkIdType> InvolvedType;

enum class Rel {
    ORIG = 1,
    DEC = 2
};

typedef std::map<vtkIdType, Rel> RelationsType;

class VTK_EXPORT vtkPolyDataBooleanFilter : public vtkPolyDataAlgorithm {
    vtkPolyData *resultA, *resultB, *contLines;
    vtkPolyData *modPdA, *modPdB;
    vtkCellData *cellDataA, *cellDataB;
    vtkIdTypeArray *cellIdsA, *cellIdsB;

    unsigned long timePdA, timePdB;

    PolyStripsType polyStripsA, polyStripsB;

    InvolvedType involvedA, involvedB;

    RelationsType relsA, relsB;

    void GetStripPoints (vtkPolyData *pd, vtkIdTypeArray *sources, PStrips &pStrips, IdsType &lines);
    bool GetPolyStrips (vtkPolyData *pd, vtkIdTypeArray *conts, vtkIdTypeArray *sources, PolyStripsType &polyStrips);
    void RemoveDuplicates (IdsType &lines);
    void CompleteStrips (PStrips &pStrips);
    bool HasArea (const StripType &strip) const;
    void CutCells (vtkPolyData *pd, PolyStripsType &polyStrips);
    void RestoreOrigPoints (vtkPolyData *pd, PolyStripsType &polyStrips);
    void DisjoinPolys (vtkPolyData *pd, PolyStripsType &polyStrips);
    void ResolveOverlaps (vtkPolyData *pd, vtkIdTypeArray *conts, PolyStripsType &polyStrips);
    void AddAdjacentPoints (vtkPolyData *pd, vtkIdTypeArray *conts, PolyStripsType &polyStrips);
    void MergePoints (vtkPolyData *pd, PolyStripsType &polyStrips);
    // void DecPolys_ (vtkPolyData *pd, InvolvedType &involved, RelationsType &rels);
    bool CombineRegions ();
    void MergeRegions ();

    int OperMode;
    bool MergeRegs, DecPolys;

public:
    vtkTypeMacro(vtkPolyDataBooleanFilter, vtkPolyDataAlgorithm);
    static vtkPolyDataBooleanFilter* New ();

    vtkSetClampMacro(OperMode, int, OPER_UNION, OPER_DIFFERENCE2);
    vtkGetMacro(OperMode, int);

    void SetOperModeToUnion () { OperMode = OPER_UNION; }
    void SetOperModeToIntersection () { OperMode = OPER_INTERSECTION; }
    void SetOperModeToDifference () { OperMode = OPER_DIFFERENCE; }
    void SetOperModeToDifference2 () { OperMode = OPER_DIFFERENCE2; }

    vtkSetMacro(MergeRegs, bool);
    vtkGetMacro(MergeRegs, bool);
    vtkBooleanMacro(MergeRegs, bool);

    vtkSetMacro(DecPolys, bool);
    vtkGetMacro(DecPolys, bool);
    vtkBooleanMacro(DecPolys, bool);

protected:
    vtkPolyDataBooleanFilter ();
    ~vtkPolyDataBooleanFilter ();

    int ProcessRequest (vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector);

private:
    vtkPolyDataBooleanFilter (const vtkPolyDataBooleanFilter&) = delete;
    void operator= (const vtkPolyDataBooleanFilter&) = delete;

};

#endif
