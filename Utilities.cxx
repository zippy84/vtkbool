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

#include "Utilities.h"

#include <cmath>

#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkDataWriter.h>
#include <vtkSmartPointer.h>

#include <vtkPolyDataWriter.h>

double ComputeNormal (vtkPoints *pts, double *n, vtkIdType num, const vtkIdType *poly) {
    n[0] = 0; n[1] = 0; n[2] = 0;

    double pA[3], pB[3];

    vtkIdType i, a, b;

    for (i = 0; i < num; i++) {
        a = poly[i];
        b = poly[i+1 == num ? 0 : i+1];

        pts->GetPoint(a, pA);
        pts->GetPoint(b, pB);

        n[0] += (pA[1]-pB[1])*(pA[2]+pB[2]);
        n[1] += (pA[2]-pB[2])*(pA[0]+pB[0]);
        n[2] += (pA[0]-pB[0])*(pA[1]+pB[1]);
    }

    return vtkMath::Normalize(n);
}

void FindPoints (vtkKdTreePointLocator *pl, const double *pt, vtkIdList *pts, double tol) {
    pts->Reset();

    vtkPolyData *pd = vtkPolyData::SafeDownCast(pl->GetDataSet());

    vtkIdList *closest = vtkIdList::New();

    pl->FindPointsWithinRadius(std::max(1e-3, tol), pt, closest);

    vtkIdType i, numPts = closest->GetNumberOfIds();

    double c[3], v[3];

    for (i = 0; i < numPts; i++) {
        pd->GetPoint(closest->GetId(i), c);
        vtkMath::Subtract(pt, c, v);

        if (vtkMath::Norm(v) < tol) {
            pts->InsertNextId(closest->GetId(i));
        }
    }

    closest->Delete();
}

#ifdef DEBUG
void WriteVTK (const char *name, vtkPolyData *pd) {
    vtkPolyDataWriter *w = vtkPolyDataWriter::New();
    w->SetInputData(pd);
    w->SetFileName(name);
    w->Update();
    w->Delete();
}
#endif

double GetAngle (const double *vA, const double *vB, const double *n) {
    // http://math.stackexchange.com/questions/878785/how-to-find-an-angle-in-range0-360-between-2-vectors

    double _vA[3];

    vtkMath::Cross(n, vA, _vA);
    double ang = std::atan2(vtkMath::Dot(_vA, vB), vtkMath::Dot(vA, vB));

    if (ang < 0) {
        ang += 2*M_PI;
    }

    return ang;
}

Base::Base (vtkPoints *pts, vtkIdType num, const vtkIdType *poly) {
    ComputeNormal(pts, n, num, poly);

    double ptA[3],
        ptB[3];

    pts->GetPoint(poly[0], ptA);
    pts->GetPoint(poly[1], ptB);

    ei[0] = ptB[0]-ptA[0];
    ei[1] = ptB[1]-ptA[1];
    ei[2] = ptB[2]-ptA[2];

    vtkMath::Normalize(ei);

    ej[0] = n[1]*ei[2]-n[2]*ei[1];
    ej[1] = -n[0]*ei[2]+n[2]*ei[0];
    ej[2] = n[0]*ei[1]-n[1]*ei[0];

    vtkMath::Normalize(ej);

    d = n[0]*ptA[0]+n[1]*ptA[1]+n[2]*ptA[2];
}

void Transform (const double *in, double *out, const Base &base) {
    double x = base.ei[0]*in[0]+base.ei[1]*in[1]+base.ei[2]*in[2],
        y = base.ej[0]*in[0]+base.ej[1]*in[1]+base.ej[2]*in[2];

    out[0] = x;
    out[1] = y;
}

void BackTransform (const double *in, double *out, const Base &base) {
    double x = in[0]*base.ei[0]+in[1]*base.ej[0]+base.d*base.n[0],
        y = in[0]*base.ei[1]+in[1]*base.ej[1]+base.d*base.n[1],
        z = in[0]*base.ei[2]+in[1]*base.ej[2]+base.d*base.n[2];

    out[0] = x;
    out[1] = y;
    out[2] = z;
}

double ComputeNormal (const Poly &poly, double *n) {
    n[0] = 0; n[1] = 0; n[2] = 0;

    Poly::const_iterator itrA, itrB;

    for (itrA = poly.begin(); itrA != poly.end(); ++itrA) {
        itrB = itrA+1;

        if (itrB == poly.end()) {
            itrB = poly.begin();
        }

        const Point3d &ptA = *itrA,
            &ptB = *itrB;

        n[0] += (ptA.y-ptB.y)*(ptA.z+ptB.z);
        n[1] += (ptA.z-ptB.z)*(ptA.x+ptB.x);
        n[2] += (ptA.x-ptB.x)*(ptA.y+ptB.y);
    }

    return vtkMath::Normalize(n);
}

bool PointInPoly (const Poly &poly, const Point3d &p) {
    bool in = false;

    Poly::const_iterator itrA, itrB;

    for (itrA = poly.begin(); itrA != poly.end(); ++itrA) {
        itrB = itrA+1;

        if (itrB == poly.end()) {
            itrB = poly.begin();
        }

        const Point3d &ptA = *itrA,
            &ptB = *itrB;

        if ((ptA.x <= p.x || ptB.x <= p.x)
            && ((ptA.y < p.y && ptB.y >= p.y)
                || (ptB.y < p.y && ptA.y >= p.y))) {

            // schnittpunkt mit bounding box und strahlensatz
            if (ptA.x+(p.y-ptA.y)*(ptB.x-ptA.x)/(ptB.y-ptA.y) < p.x) {
                in = !in;
            }
        }
    }

    return in;
}

#ifdef DEBUG
void WritePolys (const char *name, const PolysType &polys) {
    auto pts = vtkSmartPointer<vtkPoints>::New();

    auto pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(pts);
    pd->Allocate(1);

    for (auto &poly : polys) {
        auto cell = vtkSmartPointer<vtkIdList>::New();

        for (auto &p : poly) {
            cell->InsertNextId(pts->InsertNextPoint(p.x, p.y, p.z));
        }

        pd->InsertNextCell(VTK_POLYGON, cell);
    }

    WriteVTK(name, pd);
}
#endif

void GetPolys (const ReferencedPointsType &pts, const IndexedPolysType &indexedPolys, PolysType &polys) {
    for (const auto &poly : indexedPolys) {
        Poly newPoly;

        for (auto &id : poly) {
            newPoly.push_back(pts.at(id));
        }

        polys.push_back(std::move(newPoly));
    }
}

void GetPoly (vtkPoints *pts, vtkIdType num, const vtkIdType *poly, Poly &out) {
    vtkIdType i;

    double p[3];

    for (i = 0; i < num; i++) {
        pts->GetPoint(poly[i], p);

        out.emplace_back(p[0], p[1], p[2], poly[i]);
    }
}

void FlattenPoly (const Poly &poly, Poly &out, const Base &base) {
    double pt[3], tr[2];

    vtkIdType i = 0;

    for (auto &p : poly) {
        pt[0] = p.x;
        pt[1] = p.y;
        pt[2] = p.z;

        Transform(pt, tr, base);

        out.emplace_back(tr[0], tr[1], 0, i++);
    }
}

void FlattenPoly2 (const Poly &poly, Poly &out, const Base2 &base) {
    double pt[3], tr[3];

    for (auto &p : poly) {
        pt[0] = p.x;
        pt[1] = p.y;
        pt[2] = p.z;

        base.Transform(pt, tr);

        out.emplace_back(tr[0], tr[1], tr[2], p.id);
    }
}

std::shared_ptr<Proj> GetEdgeProj (const Poly &poly, const Point3d &p) {
    Poly::const_iterator itrA, itrB;

    double d, t;

    std::shared_ptr<Point3d> proj;

    for (itrA = poly.begin(); itrA != poly.end(); ++itrA) {
        itrB = itrA+1;

        if (itrB == poly.end()) {
            itrB = poly.begin();
        }

        ProjOnLine(*itrA, *itrB, p, &d, &t, proj);

        if (d > 0 && d < 1e-5 && t > 0 && t < 1) {
            return std::make_shared<Proj>(itrA->id, itrB->id, *proj, d);
        }

    }

    return nullptr;
}

void ProjOnLine (const Point3d &a, const Point3d &b, const Point3d &p, double *d, double *t, std::shared_ptr<Point3d> &proj) {
    double v[3], w[3];

    double v_ = Point3d::GetVec(a, b, v);
    double w_ = Point3d::GetVec(a, p, w);

    double pr = std::min(std::max(vtkMath::Dot(v, w), -1.), 1.);

    *d = std::sin(std::acos(pr))*w_;

    vtkMath::MultiplyScalar(v, std::sqrt(w_*w_-*d**d));

    proj = std::make_shared<Point3d>(a.x+v[0], a.y+v[1], a.z+v[2]);

    *t = pr*w_/v_;
}

void ProjOnLine (vtkPolyData *pd, const Pair &line, const Point3d &p, std::shared_ptr<Point3d> &proj) {
    double pA[3], pB[3];

    pd->GetPoint(line.f, pA);
    pd->GetPoint(line.g, pB);

    Point3d a(pA[0], pA[1], pA[2]);
    Point3d b(pB[0], pB[1], pB[2]);

    double d, t;

    ProjOnLine(a, b, p, &d, &t, proj);
}

vtkSmartPointer<vtkPolyData> CreatePolyData (const PolysType &polys) {
    auto pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetDataTypeToDouble();

    auto pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(pts);
    pd->Allocate(100);

    for (const auto &poly : polys) {
        auto cell = vtkSmartPointer<vtkIdList>::New();

        for (const auto &p : poly) {
            cell->InsertNextId(pts->InsertNextPoint(p.x, p.y, p.z));
        }

        pd->InsertNextCell(VTK_POLYGON, cell);
    }

    pd->Squeeze();

    return pd;
}
