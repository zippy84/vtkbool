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

#ifndef __Optimize_h
#define __Optimize_h

#include "Utilities.h"

typedef std::map<Pair, Points> Edges;

class SnapPoint {
public:
    SnapPoint () = delete;
    SnapPoint (vtkIdType cellId, const Pair &line, const Point3d &point, const Point3d &inter, double d) : cellId(cellId), line(line), point(point), inter(inter), d(d) {}

    vtkIdType cellId;
    Pair line;
    Point3d point;
    Point3d inter;
    double d;

    friend std::ostream& operator<< (std::ostream &out, const SnapPoint &s) {
        out << "cellId " << s.cellId
            << ", line " << s.line
            << ", point " << s.point
            << ", inter " << s.inter
            << ", d " << s.d;

        return out;
    }
};

class SnapEdge {
public:
    SnapEdge () = delete;
    SnapEdge (vtkIdType cellId, const Pair &line, const Pair &edge, const Point3d &proj, const Point3d &inter, double d) : cellId(cellId), line(line), edge(edge), proj(proj), inter(inter), d(d) {}

    vtkIdType cellId;
    Pair line;
    Pair edge;
    Point3d proj;
    Point3d inter;
    double d;

    friend std::ostream& operator<< (std::ostream &out, const SnapEdge &s) {
        out << "cellId " << s.cellId
            << ", line " << s.line
            << ", edge " << s.edge
            << ", proj " << s.proj
            << ", inter " << s.inter
            << ", d " << s.d;

        return out;
    }
};

class PreventEqualCaptPoints {
    vtkPolyData *pdA, *pdB;
public:
    static IdsType TriangulateCell (vtkPolyData *pd, vtkIdType cellId, const Edges &edges);
    static void MovePoint (vtkPolyData *pd, vtkIdType ind, const Point3d &p);

    PreventEqualCaptPoints () = delete;
    PreventEqualCaptPoints (vtkPolyData *pdA, vtkPolyData *pdB);
    void Run ();
private:
    void Find (vtkPolyData *pd, vtkPolyData *other, const std::string &name);
};

#endif
