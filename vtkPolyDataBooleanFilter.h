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

#ifndef __vtkPolyDataBooleanFilter_h
#define __vtkPolyDataBooleanFilter_h

#include <vector>
#include <deque>
#include <map>
#include <set>
#include <utility>
#include <iostream>
#include <type_traits>

#include <vtkPolyDataAlgorithm.h>
#include <vtkKdTree.h>
#include <vtkModifiedBSPTree.h>

#include "Utilities.h"

enum OperMode {
    OPER_NONE = 0,
    OPER_UNION,
    OPER_INTERSECTION,
    OPER_DIFFERENCE,
    OPER_DIFFERENCE2
};

enum class Capt {
    Not = 1 << 0,
    Edge = 1 << 1,
    A = 1 << 2,
    B = 1 << 3,
    Branched = 1 << 4,
    Boundary = 0xe
};

// inline std::underlying_type_t<Capt> operator| (Capt lhs, Capt rhs) {
//     return static_cast<std::underlying_type_t<Capt>>(lhs) |
//         static_cast<std::underlying_type_t<Capt>>(rhs);
// }

inline std::underlying_type_t<Capt> operator& (Capt lhs, Capt rhs) {
    return static_cast<std::underlying_type_t<Capt>>(lhs) &
        static_cast<std::underlying_type_t<Capt>>(rhs);
}

enum class Side {
    None,
    Start,
    End
};

enum class Loc {
    None,
    Inside,
    Outside
};

class StripPt {
public:
    StripPt () : t(0), capt(Capt::Not), catched(true) {
        edge[0] = NOTSET;
        edge[1] = NOTSET;
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

    vtkIdType polyId;

    bool catched;
};

class StripPtR {
public:
    StripPtR () = delete;

    StripPtR (vtkIdType ind, std::size_t strip) : ind(ind), strip(strip), ref(NOTSET), side(Side::None) {
        desc[0] = NOTSET;
        desc[1] = NOTSET;
    }

    vtkIdType ind;
    std::size_t strip;
    vtkIdType ref, desc[2];
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

typedef std::vector<std::reference_wrapper<StripType>> _StripsType;

class PStrips {
public:
    PStrips (vtkPolyData *pd, vtkIdType cellId) {
        const vtkIdType *cell;
        vtkIdType numPts;

        pd->GetCellPoints(cellId, numPts, cell);

        for (vtkIdType i = 0; i < numPts; i++) {
            poly.push_back(cell[i]);
        }

        ComputeNormal(pd->GetPoints(), n, numPts, cell);

        base = Base(pd->GetPoints(), numPts, cell);
    }

    StripPtsType pts;
    StripsType strips;

    double n[3];
    IdsType poly;
    Base base;
};

typedef std::map<vtkIdType, PStrips> PolyStripsType;

typedef std::vector<std::reference_wrapper<StripPtR>> RefsType;
typedef std::vector<std::reference_wrapper<const StripPtR>> ConstRefsType;

class VTK_EXPORT vtkPolyDataBooleanFilter : public vtkPolyDataAlgorithm {
    vtkPolyData *resultA, *resultB, *resultC;

    vtkSmartPointer<vtkPolyData> modPdA, modPdB, contLines;

    vtkSmartPointer<vtkCellData> cellDataA, cellDataB;
    vtkSmartPointer<vtkIdTypeArray> cellIdsA, cellIdsB;

    vtkIdTypeArray *contsA, *contsB;

    unsigned long timePdA, timePdB;

    PolyStripsType polyStripsA, polyStripsB;

    void GetStripPoints (vtkPolyData *pd, vtkIdTypeArray *sources, PStrips &pStrips, IdsType &lines);
    bool GetPolyStrips (vtkPolyData *pd, vtkIdTypeArray *conts, vtkIdTypeArray *sources, PolyStripsType &polyStrips);
    bool CleanStrips ();
    void RemoveDuplicates (IdsType &lines);
    void CompleteStrips (PStrips &pStrips);
    bool HasArea (const StripType &strip) const;
    bool CutCells (vtkPolyData *pd, PolyStripsType &polyStrips);
    void RestoreOrigPoints (vtkPolyData *pd, PolyStripsType &polyStrips);
    void DisjoinPolys (vtkPolyData *pd, PolyStripsType &polyStrips);
    void ResolveOverlaps (vtkPolyData *pd, PolyStripsType &polyStrips);
    void AddAdjacentPoints (vtkPolyData *pd, vtkIdTypeArray *conts, PolyStripsType &polyStrips);
    void MergePoints (vtkPolyData *pd, PolyStripsType &polyStrips);
    bool CombineRegions ();

    int OperMode;

public:
    vtkTypeMacro(vtkPolyDataBooleanFilter, vtkPolyDataAlgorithm);
    static vtkPolyDataBooleanFilter* New ();

    vtkSetClampMacro(OperMode, int, OPER_NONE, OPER_DIFFERENCE2);
    vtkGetMacro(OperMode, int);

    void SetOperModeToNone () { OperMode = OPER_NONE; Modified(); }
    void SetOperModeToUnion () { OperMode = OPER_UNION; Modified(); }
    void SetOperModeToIntersection () { OperMode = OPER_INTERSECTION; Modified(); }
    void SetOperModeToDifference () { OperMode = OPER_DIFFERENCE; Modified(); }
    void SetOperModeToDifference2 () { OperMode = OPER_DIFFERENCE2; Modified(); }

protected:
    vtkPolyDataBooleanFilter ();
    ~vtkPolyDataBooleanFilter ();

    int RequestData (vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;

    void PrintSelf (ostream&, vtkIndent) override {};

private:
    vtkPolyDataBooleanFilter (const vtkPolyDataBooleanFilter&) = delete;
    void operator= (const vtkPolyDataBooleanFilter&) = delete;

};

#endif
