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

#include "Optimize.h"

#include <vtkCellArrayIterator.h>
#include <vtkModifiedBSPTree.h>
#include <vtkCellData.h>
#include <vtkPolygon.h>
#include <vtkSmartPointer.h>

// #define _DEBUG

PreventEqualCaptPoints::PreventEqualCaptPoints (vtkPolyData *pdA, vtkPolyData *pdB) : pdA(pdA), pdB(pdB) {}

void PreventEqualCaptPoints::Run () {
#ifdef _DEBUG
    WriteVTK("captA.vtk", pdA);
    WriteVTK("captB.vtk", pdB);
#endif

    pdA->BuildLinks();
    pdB->BuildLinks();

    Find(pdA, pdB, "A");

#ifdef _DEBUG
    WriteVTK("modB.vtk", pdB);
#endif

    Find(pdB, pdA, "B");

#ifdef _DEBUG
    WriteVTK("modA.vtk", pdA);
#endif
}

void PreventEqualCaptPoints::Find (vtkPolyData *pd, vtkPolyData *other, [[maybe_unused]] const std::string &name) {
#ifdef _DEBUG
    std::cout << "Find(" << name << ")" << std::endl;
#endif

    vtkIdType num;
    const vtkIdType *poly;

    vtkIdType i, j;

    std::set<Pair> lines;

    auto polyItr = vtk::TakeSmartPointer(pd->GetPolys()->NewIterator());

    for (polyItr->GoToFirstCell(); !polyItr->IsDoneWithTraversal(); polyItr->GoToNextCell()) {
        polyItr->GetCurrentCell(num, poly);

        for (i = 0; i < num; i++) {
            j = i+1;

            if (j == num) {
                j = 0;
            }

            if (poly[i] < poly[j]) {
                lines.emplace(poly[i], poly[j]);
            } else {
                lines.emplace(poly[j], poly[i]);
            }
        }
    }

    auto tree = vtkSmartPointer<vtkModifiedBSPTree>::New();
    tree->SetDataSet(other);
    tree->BuildLocator();

    auto pts = vtkSmartPointer<vtkPoints>::New();
    auto cells = vtkSmartPointer<vtkIdList>::New();

    double pA[3], pB[3];

    vtkIdType cellId;

    double tr[2];

#ifdef _DEBUG
    auto pdVerts = vtkSmartPointer<vtkPolyData>::New();
    pdVerts->Allocate(1);

    auto ptsVerts = vtkSmartPointer<vtkPoints>::New();
#endif

    std::map<vtkIdType, std::vector<SnapPoint>> pointSnaps;
    std::map<Point3d, std::vector<SnapEdge>> edgeSnaps;

    for (auto &line : lines) {
        pd->GetPoint(line.f, pA);
        pd->GetPoint(line.g, pB);

        if (tree->IntersectWithLine(pA, pB, 1e-5, pts, cells) == 0) {
            continue;
        }

        for (i = 0; i < pts->GetNumberOfPoints(); i++) {
            const double *pt = pts->GetPoint(i);

            Point3d sA(pt[0], pt[1], pt[2]);

            cellId = cells->GetId(i);

            other->GetCellPoints(cellId, num, poly);

            Base base(other->GetPoints(), num, poly);

            Poly polyA, polyB;

            GetPoly(other->GetPoints(), num, poly, polyA);

            FlattenPoly(polyA, polyB, base);

            Transform(pt, tr, base);

            Point3d sB(tr[0], tr[1], 0);

            if (PointInPoly(polyB, sB)) {
#ifdef _DEBUG
                auto vert = vtkSmartPointer<vtkIdList>::New();
                vert->InsertNextId(ptsVerts->InsertNextPoint(pt));

                pdVerts->InsertNextCell(VTK_VERTEX, vert);
#endif

                // snap auf ecke oder kante?

                auto snap = std::find_if(polyA.begin(), polyA.end(), [&](const Point3d &p) { return Point3d::GetDist(p, sA) < 1e-10; });

                if (snap != polyA.end()) {
                    double d = Point3d::GetDist(*snap, sA);

                    pointSnaps[snap->id].emplace_back(cellId, line, *snap, sA, d);

                } else {
                    // projektion auf kante

                    auto edgeProj = GetEdgeProj(polyA, sA);

                    if (edgeProj != nullptr) {
                        edgeSnaps[edgeProj->proj].emplace_back(cellId, line, edgeProj->edge, edgeProj->proj, sA, edgeProj->d);
                    }
                }
            }
        }
    }

    for (const auto& [id, snaps] : pointSnaps) {
        if (snaps.size() > 1) {
#ifdef _DEBUG
            std::cout << "id " << id << ", snaps" << std::endl;

            for (const auto &s : snaps) {
                std::cout << s << std::endl;
            }
#endif

            std::set<Pair> allLines;

            std::transform(snaps.begin(), snaps.end(), std::inserter(allLines, allLines.end()), [](const SnapPoint &p) { return p.line; });

#ifdef _DEBUG
            std::cout << "allLines" << std::endl;

            for (const auto &line : allLines) {
                std::cout << line << std::endl;
            }
#endif

            auto itrA = snaps.begin();
            decltype(itrA) itrB;

            double d;

            bool collapse = true;

            for (; itrA != snaps.end()-1 && collapse; ++itrA) {
                for (itrB = itrA+1; itrB < snaps.end(); ++itrB) {
                    d = Point3d::GetDist(itrA->inter, itrB->inter);

#ifdef _DEBUG
                    std::cout << d << std::endl;
#endif

                    if (d > 1e-10) {
                        collapse = false;
                        break;
                    }
                }
            }

            if (collapse) {
                continue;
            }

            if (allLines.size() == 1) {
                std::shared_ptr<Point3d> proj;
                ProjOnLine(pd, snaps.back().line, snaps.back().point, proj);

                PreventEqualCaptPoints::MovePoint(other, id, *proj);

            } else {
                // gemeinsamen punkt ermitteln

                std::set<vtkIdType> allIds;

                for (auto &line : allLines) {
                    allIds.insert(line.f);
                    allIds.insert(line.g);
                }

                if (allIds.size() == allLines.size()+1) {
                    std::vector<Pair> _allLines(allLines.begin(), allLines.end());

                    vtkIdType s = _allLines[0] & _allLines[1];

#ifdef _DEBUG
                    std::cout << "line "
                        << _allLines[0]
                        << " & "
                        << _allLines[1]
                        << " -> "
                        << s
                        << std::endl;
#endif

                    double pt[3];
                    pd->GetPoint(s, pt);

                    Point3d p(pt[0], pt[1], pt[2]);

                    PreventEqualCaptPoints::MovePoint(other, id, p);
                } else {
                    throw std::runtime_error("");
                }
            }
        }
    }

    using Snap = std::tuple<SnapEdge, Point3d, double>;

    struct Cmp {
        bool operator() (const Snap &l, const Snap &r) const {
            return std::get<2>(l) < std::get<2>(r);
        }
    };

    std::map<Pair, std::set<Snap, Cmp>> allEdgeSnaps;

    double pt[3], t;

    for (const auto& [proj, snaps] : edgeSnaps) {
        if (snaps.size() > 1) {
            if (snaps.size() > 2) {
#ifdef _DEBUG
                std::cout << proj << std::endl;

                for (const auto &s : snaps) {
                    std::cout << s << std::endl;

                    auto &line = s.line;

                    pd->GetPoint(line.f, pA);
                    pd->GetPoint(line.g, pB);

                    std::cout << Point3d(pA[0], pA[1], pA[2]) << std::endl;
                    std::cout << Point3d(pB[0], pB[1], pB[2]) << std::endl;
                }
#endif

                continue;
            }

            const auto &snapA = snaps[0];
            const auto &snapB = snaps[1];

            if (Point3d::GetDist(snapA.inter, snapB.inter) > 1e-10) {
                Pair edge(snapA.edge);

                if (edge.f > edge.g) {
                    std::swap(edge.f, edge.g);
                }

                std::shared_ptr<Point3d> p;

                if (snapA.line == snapB.line) {
                    ProjOnLine(pd, snapA.line, snapA.proj, p);

                } else {
                    vtkIdType s = snapA.line & snapB.line;

                    pd->GetPoint(s, pt);

                    p = std::make_shared<Point3d>(pt[0], pt[1], pt[2]);
                }

                other->GetPoint(edge.f, pt);

                Point3d q(pt[0], pt[1], pt[2]);

                t = Point3d::GetDist(q, snapA.proj);

                allEdgeSnaps[edge].emplace(snapA, *p, t);
            }
        }
    }

    std::map<vtkIdType, Edges> newCells;

    for (const auto& [edge, data] : allEdgeSnaps) {
        Points pts;

        for (const auto &d : data) {
            const auto &p = std::get<1>(d);
            pts.emplace_back(p.x, p.y, p.z, other->GetPoints()->InsertNextPoint(p.x, p.y, p.z));
        }

        const auto &first = std::get<0>(*(data.begin()));

        auto neigs = vtkSmartPointer<vtkIdList>::New();

        other->GetCellEdgeNeighbors(first.cellId, first.edge.f, first.edge.g, neigs);

        if (neigs->GetNumberOfIds() != 1) {
            throw std::runtime_error("");
        }

        vtkIdType neig = neigs->GetId(0);

        auto _edge = Pair(first.edge.g, first.edge.f);

        if (edge == first.edge) {
            std::for_each(pts.begin(), pts.end(), [&](const Point3d &p) { newCells[first.cellId][first.edge].emplace_back(p); });
            std::for_each(pts.rbegin(), pts.rend(), [&](const Point3d &p) { newCells[neig][_edge].emplace_back(p); });
        } else {
            std::for_each(pts.rbegin(), pts.rend(), [&](const Point3d &p) { newCells[first.cellId][first.edge].emplace_back(p); });
            std::for_each(pts.begin(), pts.end(), [&](const Point3d &p) { newCells[neig][_edge].emplace_back(p); });
        }
    }

    for (const auto& [cellId, edges] : newCells) {
        PreventEqualCaptPoints::TriangulateCell(other, cellId, edges);
        other->DeleteCell(cellId);
    }

    other->RemoveDeletedCells();

#ifdef _DEBUG
    pdVerts->SetPoints(ptsVerts);

    auto fileName = "verts" + name + ".vtk";

    WriteVTK(fileName.c_str(), pdVerts);
#endif

}

IdsType PreventEqualCaptPoints::TriangulateCell (vtkPolyData *pd, vtkIdType cellId, const Edges &edges) {
    vtkIdTypeArray *origCellIds = vtkIdTypeArray::SafeDownCast(pd->GetCellData()->GetScalars("OrigCellIds"));

    IdsType newCellIds;

    vtkIdType num;
    const vtkIdType *poly;

    pd->GetCellPoints(cellId, num, poly);

    Poly _poly;

    double pt[3];

    vtkIdType i, j;

    vtkIdType newNum = num;

    for (i = 0; i < num; i++) {
        j = i+1;

        if (j == num) {
            j = 0;
        }

        pd->GetPoint(poly[i], pt);

        _poly.emplace_back(pt[0], pt[1], pt[2], poly[i]);

        Pair edge(poly[i], poly[j]);

        auto itr = edges.find(edge);

        if (itr != edges.end()) {
            auto &pts = itr->second;

            for (const auto &p : pts) {
                _poly.emplace_back(p);

                newNum++;
            }
        }
    }

    auto vtkPoly = vtkSmartPointer<vtkPolygon>::New();

    vtkPoly->GetPointIds()->SetNumberOfIds(newNum);
    vtkPoly->GetPoints()->SetNumberOfPoints(newNum);

    Base base(pd->GetPoints(), num, poly);

    Poly flattened;

    FlattenPoly(_poly, flattened, base);

    for (const auto &p : flattened) {
        vtkPoly->GetPointIds()->SetId(p.id, p.id);
        vtkPoly->GetPoints()->SetPoint(p.id, p.x, p.y, p.z);
    }

    auto triangles = vtkSmartPointer<vtkIdList>::New();

#if (VTK_MAJOR_VERSION >= 9 && VTK_MINOR_VERSION > 3)
    if (vtkPoly->TriangulateLocalIds(0, triangles) != 1) {
        throw std::runtime_error("");
    }
#else
    if (vtkPoly->Triangulate(triangles) != 1) {
        throw std::runtime_error("");
    }
#endif

    auto ids = vtkSmartPointer<vtkIdList>::New();

    for (const auto &p : _poly) {
        ids->InsertNextId(p.id);
    }

    vtkIdType origId = origCellIds->GetValue(cellId);

    auto triangle = vtkSmartPointer<vtkIdList>::New();
    triangle->SetNumberOfIds(3);

    for (i = 0; i < triangles->GetNumberOfIds(); i += 3) {
        triangle->SetId(0, ids->GetId(triangles->GetId(i)));
        triangle->SetId(1, ids->GetId(triangles->GetId(i+1)));
        triangle->SetId(2, ids->GetId(triangles->GetId(i+2)));

        newCellIds.push_back(pd->InsertNextCell(VTK_TRIANGLE, triangle));

        origCellIds->InsertNextValue(origId);
    }

    return newCellIds;

}

void PreventEqualCaptPoints::MovePoint (vtkPolyData *pd, vtkIdType ind, const Point3d &p) {
    auto cells = vtkSmartPointer<vtkIdList>::New();
    pd->GetPointCells(ind, cells);

    vtkIdType i, cellId;

    for (i = 0; i < cells->GetNumberOfIds(); i++) {
        cellId = cells->GetId(i);

        if (pd->GetCellType(cellId) == VTK_POLYGON) {
            PreventEqualCaptPoints::TriangulateCell(pd, cellId, {});
            pd->DeleteCell(cellId);
        }
    }

#ifdef _DEBUG
    double pt[3];
    pd->GetPoints()->GetPoint(ind, pt);

    Point3d q(pt[0], pt[1], pt[2]);

    std::cout << q
        << " -> "
        << p
        << std::endl;
#endif

    pd->GetPoints()->SetPoint(ind, p.x, p.y, p.z);
}
