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

#ifndef __Merger_h
#define __Merger_h

#include "vtkPolyDataBooleanFilter.h"

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
    Merger () = delete;
    Merger (vtkPolyData *pd, const PStrips &pStrips, const StripsType &strips, const IdsType &descIds, vtkIdType origId);
    void Run ();

private:
    void MergeGroup (const GroupType &group, PolysType &merged);
    bool FindConns (vtkPolyData *lines, vtkSmartPointer<vtkKdTree> kdTree, vtkSmartPointer<vtkModifiedBSPTree> bspTree, PolyConnsType &polyConns, const IndexedPolysType &indexedPolys, const SourcesType &sources, int &n);

    void MergeStage1 (const IndexedPolysType &indexedPolys, const ReferencedPointsType &refPts, const SourcesType &sources, const ConnsType &conns, IndexedPoly &polyA);
    void MergeStage2 (const ConnsType2 &conns, const ReferencedPointsType &refPts, const ConnsType2 &usedConns, IndexedPolysType &splitted);
};

#endif
