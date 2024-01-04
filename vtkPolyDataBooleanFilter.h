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

#define OPER_NONE 0
#define OPER_UNION 1
#define OPER_INTERSECTION 2
#define OPER_DIFFERENCE 3
#define OPER_DIFFERENCE2 4

enum class Capt {
    NOT = 1 << 0,
    EDGE = 1 << 1,
    A = 1 << 2,
    B = 1 << 3,
    BRANCHED = 1 << 4,
    BOUNDARY = 0xe
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

    StripPtR (vtkIdType ind, std::size_t strip) : ind(ind), strip(strip), ref(NOTSET), side(Side::NONE) {
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

// Merger

typedef std::vector<std::size_t> GroupType;

typedef std::map<vtkIdType, std::size_t> SourcesType;

class Conn {
public:
    Conn () = delete;
    Conn (double d, vtkIdType i, vtkIdType j) : d(d), i(i), j(j) {}

    double d;
    vtkIdType i, j;

    bool operator< (const Conn &other) const {
        return d < other.d;
    }

    friend std::ostream& operator<< (std::ostream &out, const Conn &c) {
        out << "Conn(d=" << c.d
            << ", i=" << c.i
            << ", j=" << c.j
            << ")";
        return out;
    }

};

struct ConnCmp {
    bool operator() (const Conn &a, const Conn &b) const {
        return std::tie(a.i, a.j) < std::tie(b.i, b.j);
    }
};

typedef std::vector<Conn> ConnsType;
typedef std::map<std::size_t, ConnsType> PolyConnsType;

typedef std::set<Conn, ConnCmp> ConnsType2;

inline std::ostream& operator<< (std::ostream &out, const PolyConnsType& polyConns) {
    PolyConnsType::const_iterator itr;

    for (itr = polyConns.begin(); itr != polyConns.end(); ++itr) {
        out << itr->first << ": [";
        for (auto &conn : itr->second) {
            out << conn << ", ";
        }
        out << "]" << std::endl;
    }

    return out;
}

class Prio {
public:
    Prio () = delete;
    Prio (const Conn &conn, const std::set<std::size_t> &solvable, double d) : conn(conn), solvable(solvable), d(d) {}

    Conn conn;
    std::set<std::size_t> solvable;
    double d;

    friend std::ostream& operator<< (std::ostream &out, const Prio &p) {
        out << "Prio(conn=" << p.conn
            << ", d=" << p.d
            << ")";
        return out;
    }
};

struct Cmp {
    bool operator() (const Prio &a, const Prio &b) const {
        const auto _a = a.solvable.size(),
            _b = b.solvable.size();
        return std::tie(_a, a.d) < std::tie(_b, b.d);
    }
};

typedef std::set<Prio, Cmp> PriosType;

typedef std::map<std::size_t, Prio> PolyPriosType;

class Merger {
    vtkPolyData *pd;
    const PStrips &pStrips;
    vtkIdType origId;

    PolysType polys;
    std::vector<std::size_t> innerIds;
public:
    Merger (vtkPolyData *pd, const PStrips &pStrips, const StripsType &strips, const IdsType &descIds, vtkIdType origId);
    void Run ();

private:
    void MergeGroup (const GroupType &group, PolysType &merged);
    bool FindConns (vtkPolyData *lines, vtkSmartPointer<vtkKdTree> kdTree, vtkSmartPointer<vtkModifiedBSPTree> bspTree, PolyConnsType &polyConns, const IndexedPolysType &indexedPolys, const SourcesType &sources, int &n);

    void MergeStage1 (const IndexedPolysType &indexedPolys, const ReferencedPointsType &refPts, const SourcesType &sources, const ConnsType &conns, IndexedPoly &polyA);
    void MergeStage2 (const ConnsType2 &conns, const ReferencedPointsType &refPts, const ConnsType2 &usedConns, IndexedPolysType &splitted);
};

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
