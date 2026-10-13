/*
Copyright 2012-2026 Ronald Römer

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

#include "Optimize2.h"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>

#include <vector>
#include <algorithm>
#include <execution>
#include <cmath>

struct Vector3 {
    double x, y, z;
};

struct Plane {
    Vector3 normal;
    Vector3 centroid;
};

inline bool Solve3x3 (const double A[3][3], const double b[3], Vector3 &out) {
    double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
       - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
       + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

    // lösbar/singulär?
    if (std::abs(det) < 1e-8) {
        return false;
    }

    double invDet = 1.0 / det;

    out.x = (b[0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
        - A[0][1] * (b[1] * A[2][2] - A[1][2] * b[2])
        + A[0][2] * (b[1] * A[2][1] - A[1][1] * b[2])) * invDet;

    out.y = (A[0][0] * (b[1] * A[2][2] - A[1][2] * b[2])
        - b[0] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
        + A[0][2] * (A[1][0] * b[2] - b[1] * A[2][0])) * invDet;

    out.z = (A[0][0] * (A[1][1] * b[2] - b[1] * A[2][1])
        - A[0][1] * (A[1][0] * b[2] - b[1] * A[2][0])
        + b[0] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])) * invDet;

    return true;
}

void PlanarizeMesh (vtkPolyData *pd, int maxIter) {
    if (!pd) {
        return;
    }

    vtkPoints *points = pd->GetPoints();
    if (!points) {
        return;
    }

    const vtkIdType numPoints = points->GetNumberOfPoints();
    const vtkIdType numFaces = pd->GetNumberOfPolys();

    if (numPoints == 0 || numFaces == 0) {
        return;
    }

    std::vector<std::vector<vtkIdType>> faces(numFaces);
    std::vector<std::vector<vtkIdType>> pointToFaces(numPoints);

    std::vector<bool> polygonFaces(numFaces, false);

    vtkCellArray *polys = pd->GetPolys();
    auto ids = vtkSmartPointer<vtkIdList>::New();

    polys->InitTraversal();
    for (vtkIdType i = 0; i < numFaces; i++) {
        polys->GetNextCell(ids);
        vtkIdType npts = ids->GetNumberOfIds();
        faces[i].reserve(npts);

        for (vtkIdType j = 0; j < npts; j++) {
            vtkIdType id = ids->GetId(j);
            faces[i].push_back(id);
            pointToFaces[id].push_back(i);
        }

        if (npts > 3) {
            polygonFaces[i] = true;
        }
    }

    std::vector<Vector3> oldVerts(numPoints);
    std::vector<Vector3> newVerts(numPoints);

    double p[3];

    for (vtkIdType i = 0; i < numPoints; i++) {
        points->GetPoint(i, p);
        oldVerts[i] = {p[0], p[1], p[2]};
    }

    newVerts = oldVerts;

    std::vector<size_t> indsVerts(numPoints);
    for (vtkIdType i = 0; i < numPoints; i++) {
        indsVerts[i] = i;
    }

    std::vector<size_t> indsFaces(numFaces);
    for (vtkIdType i = 0; i < numFaces; i++) {
        indsFaces[i] = i;
    }

    std::vector<Plane> planes(numFaces);

    for (int iter = 0; iter < maxIter; iter++) {
        std::for_each(std::execution::par, indsFaces.begin(), indsFaces.end(),
            [&](size_t f) {
                const auto &face = faces[f];
                size_t n = face.size();

                Vector3 centroid{0, 0, 0};
                Vector3 normal{0, 0, 0};

                for (size_t i = 0; i < n; i++) {
                    const Vector3 &curr = oldVerts[face[i]];
                    const Vector3 &next = oldVerts[face[(i + 1) % n]];

                    centroid.x += curr.x;
                    centroid.y += curr.y;
                    centroid.z += curr.z;

                    normal.x += (curr.y - next.y) * (curr.z + next.z);
                    normal.y += (curr.z - next.z) * (curr.x + next.x);
                    normal.z += (curr.x - next.x) * (curr.y + next.y);
                }

                centroid.x /= static_cast<double>(n);
                centroid.y /= static_cast<double>(n);
                centroid.z /= static_cast<double>(n);

                double len = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

                if (len > 1e-12) {
                    normal.x /= len; normal.y /= len; normal.z /= len;
                } else {
                    // kann aber nicht auftreten
                    normal = {0, 0, 1};
                }

                planes[f] = {normal, centroid};
            });

        std::for_each(std::execution::par, indsVerts.begin(), indsVerts.end(),
            [&](size_t v) {
                const auto &neigs = pointToFaces[v];
                if (neigs.empty()) {
                    return;
                }

                bool hasPolygon = false;
                for (vtkIdType neig : neigs) {
                    if (polygonFaces[neig]) {
                        hasPolygon = true;
                        break;
                    }
                }

                if (!hasPolygon) {
                    newVerts[v] = oldVerts[v];
                    return;
                }

                double ATA[3][3] = {0};
                double ATb[3] = {0};

                for (vtkIdType neig : neigs) {
                    const Plane &plane = planes[neig];
                    const double nx = plane.normal.x;
                    const double ny = plane.normal.y;
                    const double nz = plane.normal.z;
                    const double d  = nx * plane.centroid.x + ny * plane.centroid.y + nz * plane.centroid.z;

                    ATA[0][0] += nx * nx; ATA[0][1] += nx * ny; ATA[0][2] += nx * nz;
                    ATA[1][0] += ny * nx; ATA[1][1] += ny * ny; ATA[1][2] += ny * nz;
                    ATA[2][0] += nz * nx; ATA[2][1] += nz * ny; ATA[2][2] += nz * nz;

                    ATb[0] += nx * d;
                    ATb[1] += ny * d;
                    ATb[2] += nz * d;
                }

                // tikhonov-regularisierung
                // sanfte dämpfung auf diagonale hinzufügen, um rangverlust bei flachen/parallelen flächen zu verhindern

                const double lambda = 1e-4;
                const Vector3 &currVert = oldVerts[v];

                ATA[0][0] += lambda;
                ATA[1][1] += lambda;
                ATA[2][2] += lambda;

                ATb[0] += lambda * currVert.x;
                ATb[1] += lambda * currVert.y;
                ATb[2] += lambda * currVert.z;

                Vector3 newVert;
                if (Solve3x3(ATA, ATb, newVert)) {
                    // dämpfung für stetige konvergenz

                    newVerts[v].x = currVert.x + 0.5 * (newVert.x - currVert.x);
                    newVerts[v].y = currVert.y + 0.5 * (newVert.y - currVert.y);
                    newVerts[v].z = currVert.z + 0.5 * (newVert.z - currVert.z);
                } else {
                    newVerts[v] = currVert;
                }
            });

        oldVerts.swap(newVerts);
    }

    for (vtkIdType i = 0; i < numPoints; ++i) {
        points->SetPoint(i, oldVerts[i].x, oldVerts[i].y, oldVerts[i].z);
    }

    points->Modified();
    pd->Modified();
}

double ComputeMaxPlanarityError (vtkPolyData *pd) {
    double maxError = 0;

    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *polys = pd->GetPolys();

    auto ids = vtkSmartPointer<vtkIdList>::New();
    polys->InitTraversal();

    while (polys->GetNextCell(ids)) {
        vtkIdType npts = ids->GetNumberOfIds();

        if (npts < 4) {
            continue;
        }

        Vector3 centroid{0, 0, 0};
        Vector3 normal{0, 0, 0};

        double p1[3], p2[3];

        for (vtkIdType i = 0; i < npts; i++) {
            pts->GetPoint(ids->GetId(i), p1);
            pts->GetPoint(ids->GetId((i + 1) % npts), p2);

            centroid.x += p1[0]; centroid.y += p1[1]; centroid.z += p1[2];

            normal.x += (p1[1] - p2[1]) * (p1[2] + p2[2]);
            normal.y += (p1[2] - p2[2]) * (p1[0] + p2[0]);
            normal.z += (p1[0] - p2[0]) * (p1[1] + p2[1]);
        }

        centroid.x /= static_cast<double>(npts);
        centroid.y /= static_cast<double>(npts);
        centroid.z /= static_cast<double>(npts);

        double len = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

        if (len > 1e-12) {
            normal.x /= len; normal.y /= len; normal.z /= len;
        }

        double p[3];

        for (vtkIdType i = 0; i < npts; i++) {
            pts->GetPoint(ids->GetId(i), p);

            double dist = std::abs(normal.x * (p[0] - centroid.x) + normal.y * (p[1] - centroid.y) + normal.z * (p[2] - centroid.z));

            maxError = std::max(maxError, dist);
        }
    }

    return maxError;
}
