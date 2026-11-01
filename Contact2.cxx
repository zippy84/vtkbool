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

#include "Contact2.h"

#include <iostream>
#include <vector>
#include <execution>
#include <map>
#include <set>
#include <stack>
#include <optional>
#include <array>
#include <tuple>

#include <vtkSmartPointer.h>
#include <vtkOBBTree.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkCellData.h>
#include <vtkCleanPolyData.h>

#include <vtkSphereSource.h>
#include <vtkTransformFilter.h>
#include <vtkTransform.h>
#include <vtkStripper.h>
#include <vtkCellIterator.h>
#include <vtkVersionMacros.h>

std::ostream& operator<< (std::ostream &out, const Vector3 &vec) {
    out << "[" << vec.x << ", " << vec.y << ", " << vec.z << "]";
    return out;
}

int Intersect (vtkOBBNode *nodeA, vtkOBBNode *nodeB, vtkMatrix4x4 *, void *data) {
    auto* intersections = static_cast<std::vector<Pair2>*>(data);

    vtkIdList *cellsA = nodeA->Cells;
    vtkIdList *cellsB = nodeB->Cells;

    vtkIdType numCellsA = cellsA->GetNumberOfIds();
    vtkIdType numCellsB = cellsB->GetNumberOfIds();

    vtkIdType i, j, cellA, cellB;

    for (i = 0; i < numCellsA; i++) {
        cellA = cellsA->GetId(i);

        for (j = 0; j < numCellsB; j++) {
            cellB = cellsB->GetId(j);

            intersections->emplace_back(Pair2{cellA, cellB});
        }
    }

    return 1;
}

Mesh GetMesh (vtkPolyData *pd) {
    auto clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetOutputPointsPrecision(vtkAlgorithm::DOUBLE_PRECISION);
    clean->SetToleranceIsAbsolute(true);
    clean->SetAbsoluteTolerance(1e-6);
    clean->ConvertLinesToPointsOff();
    clean->ConvertPolysToLinesOff();
    clean->ConvertStripsToPolysOff();
    clean->SetInputData(pd);
    clean->Update();

    vtkPolyData *cleanedPd = clean->GetOutput();

    vtkPoints *points = cleanedPd->GetPoints();

    if (!points) {
        return {};
    }

    const vtkIdType numPoints = points->GetNumberOfPoints();

    if (numPoints == 0) {
        return {};
    }

    Mesh mesh;
    mesh.pointToFaces.resize(numPoints);

    vtkIdType faceId{0};

    vtkIdType cellId;
    vtkIdList *pts;
    vtkIdType npts;

    auto AddTriangle = [&](vtkIdType id0, vtkIdType id1, vtkIdType id2) {
        Face face;
        face.origCellId = cellId;
        face.vertIds = {id0, id1, id2};

        mesh.pointToFaces[id0].push_back(faceId);
        mesh.pointToFaces[id1].push_back(faceId);
        mesh.pointToFaces[id2].push_back(faceId);

        mesh.faces.push_back(std::move(face));
        faceId++;
    };

    auto GetDistance = [](const double a[3], const double b[3]) {
        double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
        return dx * dx + dy * dy + dz * dz;
    };

    auto IsNonConvex = [](const double a[3], const double b[3], const double c[3], const double nx, const double ny, const double nz) {
        double ax = b[0] - a[0], ay = b[1] - a[1], az = b[2] - a[2];
        double bx = c[0] - b[0], by = c[1] - b[1], bz = c[2] - b[2];

        double cx = ay * bz - az * by;
        double cy = az * bx - ax * bz;
        double cz = ax * by - ay * bx;

        return cx * nx + cy * ny + cz * nz < 0;
    };

    vtkCellIterator *cellItr = cleanedPd->NewCellIterator();

    for (cellItr->InitTraversal(); !cellItr->IsDoneWithTraversal(); cellItr->GoToNextCell()) {
        cellId = cellItr->GetCellId();
        pts = cellItr->GetPointIds();
        npts = pts->GetNumberOfIds();

        if (cellItr->GetCellType() == VTK_TRIANGLE || cellItr->GetCellType() == VTK_POLYGON) {
            Face face;
            face.origCellId = cellId;
            face.vertIds.resize(npts);

            for(vtkIdType i = 0; i < npts; i++) {
                face.vertIds[i] = pts->GetId(i);
                mesh.pointToFaces[pts->GetId(i)].push_back(faceId);
            }

            mesh.faces.push_back(std::move(face));

            faceId++;
        } else if (cellItr->GetCellType() == VTK_QUAD) {
            vtkIdType v0 = pts->GetId(0);
            vtkIdType v1 = pts->GetId(1);
            vtkIdType v2 = pts->GetId(2);
            vtkIdType v3 = pts->GetId(3);

            double p0[3], p1[3], p2[3], p3[3];
            points->GetPoint(v0, p0);
            points->GetPoint(v1, p1);
            points->GetPoint(v2, p2);
            points->GetPoint(v3, p3);

            double nx = (p0[1] - p1[1]) * (p0[2] + p1[2]) +
                (p1[1] - p2[1]) * (p1[2] + p2[2]) +
                (p2[1] - p3[1]) * (p2[2] + p3[2]) +
                (p3[1] - p0[1]) * (p3[2] + p0[2]);

            double ny = (p0[2] - p1[2]) * (p0[0] + p1[0]) +
                (p1[2] - p2[2]) * (p1[0] + p2[0]) +
                (p2[2] - p3[2]) * (p2[0] + p3[0]) +
                (p3[2] - p0[2]) * (p3[0] + p0[0]);

            double nz = (p0[0] - p1[0]) * (p0[1] + p1[1]) +
                (p1[0] - p2[0]) * (p1[1] + p2[1]) +
                (p2[0] - p3[0]) * (p2[1] + p3[1]) +
                (p3[0] - p0[0]) * (p3[1] + p0[1]);

            double lenSq = nx * nx + ny * ny + nz * nz;

            if (lenSq > 1e-12) {
                double len = std::sqrt(lenSq);

                nx /= len;
                ny /= len;
                nz /= len;
            }

            double cx = (p0[0] + p1[0] + p2[0] + p3[0]) / 4;
            double cy = (p0[1] + p1[1] + p2[1] + p3[1]) / 4;
            double cz = (p0[2] + p1[2] + p2[2] + p3[2]) / 4;

            double d = nx * cx + ny * cy + nz * cz;

            double d0 = nx * p0[0] + ny * p0[1] + nz * p0[2] - d;
            double d1 = nx * p1[0] + ny * p1[1] + nz * p1[2] - d;
            double d2 = nx * p2[0] + ny * p2[1] + nz * p2[2] - d;
            double d3 = nx * p3[0] + ny * p3[1] + nz * p3[2] - d;

            if (std::abs(d0) < 1e-6 && std::abs(d1) < 1e-6 && std::abs(d2) < 1e-6 && std::abs(d3) < 1e-6) {
                Face face;
                face.origCellId = cellId;
                face.vertIds = {v0, v1, v2, v3};

                mesh.pointToFaces[v0].push_back(faceId);
                mesh.pointToFaces[v1].push_back(faceId);
                mesh.pointToFaces[v2].push_back(faceId);
                mesh.pointToFaces[v3].push_back(faceId);

                mesh.faces.push_back(std::move(face));
                faceId++;
            } else {
                if (IsNonConvex(p3, p0, p1, nx, ny, nz)) {
                    AddTriangle(v0, v1, v2);
                    AddTriangle(v2, v3, v0);
                } else if (IsNonConvex(p0, p1, p2, nx, ny, nz)) {
                    AddTriangle(v1, v2, v3);
                    AddTriangle(v3, v0, v1);
                } else if (IsNonConvex(p1, p2, p3, nx, ny, nz)) {
                    AddTriangle(v2, v3, v0);
                    AddTriangle(v0, v1, v2);
                } else if (IsNonConvex(p2, p3, p0, nx, ny, nz)) {
                    AddTriangle(v3, v0, v1);
                    AddTriangle(v1, v2, v3);
                } else {
                    if (GetDistance(p0, p2) < GetDistance(p1, p3)) {
                        AddTriangle(v0, v1, v2);
                        AddTriangle(v2, v3, v0);
                    } else {
                        AddTriangle(v1, v2, v3);
                        AddTriangle(v3, v0, v1);
                    }
                }
            }
        } else if (cellItr->GetCellType() == VTK_TRIANGLE_STRIP) {
            for (vtkIdType i = 0; i < npts - 2; i++) {
                vtkIdType v0 = pts->GetId(i);
                vtkIdType v1 = pts->GetId(i + 1);
                vtkIdType v2 = pts->GetId(i + 2);

                if (v0 == v1 || v1 == v2 || v2 == v0) {
                    continue;
                }

                if (i % 2 == 0) {
                    AddTriangle(v0, v1, v2);
                } else {
                    AddTriangle(v0, v2, v1);
                }
            }
        }
    }

    cellItr->Delete();

    mesh.verts.resize(numPoints);

    double p[3];

    for (vtkIdType i = 0; i < numPoints; i++) {
        points->GetPoint(i, p);
        mesh.verts[i] = {p[0], p[1], p[2]};
    }

    return mesh;
}

void UpdatePlanes (Mesh &mesh) {
    std::for_each(std::execution::par, mesh.faces.begin(), mesh.faces.end(),
        [&](Face &face) {
            size_t numVerts = face.vertIds.size();

            face.normal.x = face.normal.y = face.normal.z = face.centroid.x = face.centroid.y = face.centroid.z = 0;

            for (size_t i = 0; i < numVerts; i++) {
                const auto &pA = mesh.verts[face.vertIds[i]];
                const auto &pB = mesh.verts[face.vertIds[(i + 1) % numVerts]];

                face.centroid.x += pA.x;
                face.centroid.y += pA.y;
                face.centroid.z += pA.z;

                face.normal.x += (pA.y - pB.y) * (pA.z + pB.z);
                face.normal.y += (pA.z - pB.z) * (pA.x + pB.x);
                face.normal.z += (pA.x - pB.x) * (pA.y + pB.y);
            }

            face.centroid.x /= static_cast<double>(numVerts);
            face.centroid.y /= static_cast<double>(numVerts);
            face.centroid.z /= static_cast<double>(numVerts);

            double lenSq = face.normal.x * face.normal.x + face.normal.y * face.normal.y + face.normal.z * face.normal.z;

            if (lenSq > 1e-12) {
                double len = std::sqrt(lenSq);
                face.normal.x /= len;
                face.normal.y /= len;
                face.normal.z /= len;
            }

            face.d = face.normal.x * face.centroid.x + face.normal.y * face.centroid.y + face.normal.z * face.centroid.z;
        });
}

inline bool Solve4x4 (double A[4][4], const double b[4], double x[4]) {
    double aug[4][5];

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            aug[i][j] = A[i][j];
        }
        aug[i][4] = b[i];
    }

    for (int i = 0; i < 4; i++) {
        int maxRow = i;

        for (int k = i + 1; k < 4; k++) {
            if (std::abs(aug[k][i]) > std::abs(aug[maxRow][i])) {
                maxRow = k;
            }
        }

        if (maxRow != i) {
            for (int k = i; k < 5; k++) {
                std::swap(aug[i][k], aug[maxRow][k]);
            }
        }

        if (std::abs(aug[i][i]) < 1e-12) {
            return false;
        }

        for (int k = i + 1; k < 4; k++) {
            double factor = aug[k][i] / aug[i][i];
            for (int j = i; j < 5; j++) {
                aug[k][j] -= factor * aug[i][j];
            }
        }
    }

    for (int i = 3; i >= 0; i--) {
        x[i] = aug[i][4];
        for (int j = i + 1; j < 4; j++) {
            x[i] -= aug[i][j] * x[j];
        }
        x[i] /= aug[i][i];
    }

    return true;
}

bool FitPlane (const std::vector<Vector3> &points, Vector3 &normal, double &d, int maxIters = 10) {
    if (points.size() < 3) {
        return false;
    }

    const double w = 100.0;

    for (int iter = 0; iter < maxIters; iter++) {
        double M[4][4] = {0.0};
        double rhs[4] = {0.0};

        for (const auto &p : points) {
            double rx = p.x;
            double ry = p.y;
            double rz = p.z;

            double r = normal.x * rx + normal.y * ry + normal.z * rz - d;

            M[0][0] += rx * rx; M[0][1] += rx * ry; M[0][2] += rx * rz; M[0][3] -= rx;
            M[1][0] += ry * rx; M[1][1] += ry * ry; M[1][2] += ry * rz; M[1][3] -= ry;
            M[2][0] += rz * rx; M[2][1] += rz * ry; M[2][2] += rz * rz; M[2][3] -= rz;
            M[3][0] -= rx; M[3][1] -= ry; M[3][2] -= rz; M[3][3] += 1.0;

            rhs[0] -= r * rx;
            rhs[1] -= r * ry;
            rhs[2] -= r * rz;
            rhs[3] -= r * (-1.0);
        }

        M[0][0] += w * normal.x * normal.x; M[0][1] += w * normal.x * normal.y; M[0][2] += w * normal.x * normal.z;
        M[1][0] += w * normal.y * normal.x; M[1][1] += w * normal.y * normal.y; M[1][2] += w * normal.y * normal.z;
        M[2][0] += w * normal.z * normal.x; M[2][1] += w * normal.z * normal.y; M[2][2] += w * normal.z * normal.z;

        double delta[4];

        if (!Solve4x4(M, rhs, delta)) {
            return false;
        }

        normal.x += delta[0];
        normal.y += delta[1];
        normal.z += delta[2];
        d += delta[3];

        double lenSq = normal.x * normal.x + normal.y * normal.y + normal.z * normal.z;

        if (lenSq > 1e-12) {
            double len = std::sqrt(lenSq);
            normal.x /= len;
            normal.y /= len;
            normal.z /= len;
        }

        if (std::sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2] + delta[3] * delta[3]) < 1e-8) {
            break;
        }
    }

    return true;
}

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

double FindNewVerts (Mesh &mesh, std::vector<GlobalPlane> &globalPlanes) {
    std::vector<size_t> indsVerts(mesh.verts.size());

    for (size_t i = 0; i < mesh.verts.size(); i++) {
        indsVerts[i] = i;
    }

    std::vector<Vector3> newVerts = mesh.verts;

    std::for_each(std::execution::par, indsVerts.begin(), indsVerts.end(),
        [&](size_t v) {
            const auto &neigs = mesh.pointToFaces[v];
            if (neigs.empty()) {
                return;
            }

            double ATA[3][3] = {0};
            double ATb[3] = {0};

            for (vtkIdType neig : neigs) {
                const Face &face = mesh.faces[neig];

                if (face.globalPlaneId) {
                    const auto &plane = globalPlanes[*face.globalPlaneId];

                    const double nx = plane.normal.x;
                    const double ny = plane.normal.y;
                    const double nz = plane.normal.z;
                    const double d  = plane.d;

                    ATA[0][0] += nx * nx; ATA[0][1] += nx * ny; ATA[0][2] += nx * nz;
                    ATA[1][0] += ny * nx; ATA[1][1] += ny * ny; ATA[1][2] += ny * nz;
                    ATA[2][0] += nz * nx; ATA[2][1] += nz * ny; ATA[2][2] += nz * nz;

                    ATb[0] += nx * d;
                    ATb[1] += ny * d;
                    ATb[2] += nz * d;
                } else {
                    const double nx = face.normal.x;
                    const double ny = face.normal.y;
                    const double nz = face.normal.z;
                    const double d  = face.d;

                    ATA[0][0] += nx * nx; ATA[0][1] += nx * ny; ATA[0][2] += nx * nz;
                    ATA[1][0] += ny * nx; ATA[1][1] += ny * ny; ATA[1][2] += ny * nz;
                    ATA[2][0] += nz * nx; ATA[2][1] += nz * ny; ATA[2][2] += nz * nz;

                    ATb[0] += nx * d;
                    ATb[1] += ny * d;
                    ATb[2] += nz * d;
                }
            }

            // tikhonov-regularisierung
            // sanfte dämpfung auf diagonale hinzufügen, um rangverlust bei flachen/parallelen flächen zu verhindern

            const double lambda = 1e-4;
            const Vector3 &currVert = mesh.verts[v];

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

    double maxDisp = std::transform_reduce(
        std::execution::par,
        indsVerts.begin(), indsVerts.end(),
        0.0,
        [](double a, double b) { return std::max(a, b); },
        [&](size_t v) {
            double dx = newVerts[v].x - mesh.verts[v].x;
            double dy = newVerts[v].y - mesh.verts[v].y;
            double dz = newVerts[v].z - mesh.verts[v].z;
            return std::sqrt(dx * dx + dy * dy + dz * dz);
        }
    );

    mesh.verts.swap(newVerts);

    return maxDisp;
}

void UpdatePoints (vtkPolyData *pd, const Mesh &mesh) {
    vtkPoints *points = pd->GetPoints();

    for (size_t i = 0; i < mesh.verts.size(); i++) {
        points->SetPoint(static_cast<vtkIdType>(i), mesh.verts[i].x, mesh.verts[i].y, mesh.verts[i].z);
    }

    points->Modified();
    pd->Modified();
}

class Base {
public:
    Base () {}
    Base (const double *t, const double *ei, const double *ek) : t{t[0], t[1], t[2]}, ei{ei[0], ei[1], ei[2]}, ek{ek[0], ek[1], ek[2]} {
        ej[0] = ek[1]*ei[2]-ek[2]*ei[1];
        ej[1] = -ek[0]*ei[2]+ek[2]*ei[0];
        ej[2] = ek[0]*ei[1]-ek[1]*ei[0];

        vtkMath::Normalize(ej.data());
    }
    void Transform (const double *in, double *out) const {
        double R[][3] = {
            { ei[0], ei[1], ei[2] },
            { ej[0], ej[1], ej[2] },
            { ek[0], ek[1], ek[2] }
        };

        double _t[] = { in[0]-t[0], in[1]-t[1], in[2]-t[2] };

        out[0] = R[0][0]*_t[0]+R[0][1]*_t[1]+R[0][2]*_t[2];
        out[1] = R[1][0]*_t[0]+R[1][1]*_t[1]+R[1][2]*_t[2];
        out[2] = R[2][0]*_t[0]+R[2][1]*_t[1]+R[2][2]*_t[2];
    }
    void BackTransform (const double *in, double *out) const {
        double R[][3] = {
            { ei[0], ej[0], ek[0] },
            { ei[1], ej[1], ek[1] },
            { ei[2], ej[2], ek[2] }
        };

        double _out[] = {
            R[0][0]*in[0]+R[0][1]*in[1]+R[0][2]*in[2],
            R[1][0]*in[0]+R[1][1]*in[1]+R[1][2]*in[2],
            R[2][0]*in[0]+R[2][1]*in[1]+R[2][2]*in[2]
        };

        out[0] = _out[0]+t[0];
        out[1] = _out[1]+t[1];
        out[2] = _out[2]+t[2];
    }

    std::array<double, 3> t, ei, ej, ek;
};

struct TransformedVert {
    double x, y;
    vtkIdType id;
};

std::vector<TransformedVert> TransformFace2d (const Mesh &mesh, const Face &face, const Base &base) {
    std::vector<TransformedVert> newVerts;

    newVerts.reserve(face.vertIds.size());

    double pt[3], tr[3];

    for (auto &id : face.vertIds) {
        const auto &vert = mesh.verts[id];

        pt[0] = vert.x; pt[1] = vert.y; pt[2] = vert.z;

        base.Transform(pt, tr);

        newVerts.push_back(TransformedVert{tr[0], tr[1], id});
    }

    return newVerts;
}

void IntersectEdgeAndLine (InterPts &interPts, const TransformedVert &pA, const TransformedVert &pB, Src src) {
    // schnitt mit x-achse

    // std::cout << "(" << pA.id << ", " << pB.id << ")" << std::endl;

    // std::cout << "pA.y " << pA.y << ", pB.y " << pB.y << std::endl;

    if (std::abs(pA.y) < 1e-5 && std::abs(pB.y) < 1e-5) {
        // std::cout << "-> 1" << std::endl;

        interPts.emplace_back(pA.x, 0, 0, pA.x, pA.id, pB.id, End::A, src, PointSrc::Copied);
        interPts.emplace_back(pB.x, 0, 0, pB.x, pA.id, pB.id, End::B, src, PointSrc::Copied);

        // std::cout << "pA.x " << pA.x << ", pB.x " << pB.x << std::endl;

        return;
    }

    double yA = pA.y;
    double yB = pB.y;

    if (std::abs(yA) < 1e-6) {
        yA = 0;
    }

    if (std::abs(yB) < 1e-6) {
        yB = 0;
    }

    if (std::min(yA, yB) < 1e-6 && std::max(yA, yB) > -1e-6) {
        double div = yB - yA;

        if (std::abs(div) < 1e-12) {
            return;
        }

        double t = -yA / div;

        // std::cout << "t " << t << std::endl;

        if (t > -1e-6 && t < 1+1e-6) {
            // std::cout << "-> 2" << std::endl;

            double x = pA.x + t * (pB.x - pA.x);

            // std::cout << "x " << x << std::endl;

            End end = End::None;

            if (std::abs(t) < 1e-6) {
                end = End::A;
            } else if (std::abs(t - 1) < 1e-6) {
                end = End::B;
            }

            interPts.emplace_back(x, 0, 0, x, pA.id, pB.id, end, src);
        }
    }
}

InterPts IntersectFaceAndLine (const Base &base, const std::vector<TransformedVert> &verts, Src src) {
    InterPts interPts;

    size_t numVerts = verts.size();

    for (size_t i = 0; i < numVerts; i++) {
        const auto &pA = verts[i];
        const auto &pB = verts[(i + 1) % numVerts];

        IntersectEdgeAndLine(interPts, pA, pB, src);
    }

    std::sort(interPts.begin(), interPts.end(), [](auto &a, auto &b) { return a.t < b.t; });

    struct Cmp {
        const double tol = 1e-5;

        bool operator() (const InterPt &lhs, const InterPt &rhs) const {
            auto l = lhs.GetEnd();
            auto r = rhs.GetEnd();

            if (lhs.pointSrc == PointSrc::Copied && rhs.pointSrc == PointSrc::Copied) {
                if (l && r && *l == *r) {
                    return false;
                }
            } else if (lhs.pointSrc == PointSrc::Copied) {
                if (l) {
                    if (*l == rhs.edge.a || *l == rhs.edge.b) {
                        return false;
                    }
                }
            } else if (rhs.pointSrc == PointSrc::Copied) {
                if (r) {
                    if (*r == lhs.edge.a || *r == lhs.edge.b) {
                        return false;
                    }
                }
            }

            const double d = lhs.t-rhs.t;

            if (std::abs(d) < tol) {
                return false;
            } else {
                return d < 0;
            }
        }
    };

    std::map<InterPt, InterPts, Cmp> groupedInterPts;

    for (auto &p : interPts) {
        groupedInterPts[p].push_back(p);
    }

    std::map<vtkIdType, double> positions;

    std::vector<InterPts> sortedInterPts;

    for (auto& [k, v] : groupedInterPts) {
        if (v.back().end == End::A) {
            positions.emplace(v.back().edge.a, v.back().t);
        } else if (v.back().end == End::B) {
            positions.emplace(v.back().edge.b, v.back().t);
        }

        sortedInterPts.push_back(std::move(v));
    }

    std::map<vtkIdType, TransformedVert> allVerts;

    for (auto &v : verts) {
        allVerts.emplace(v.id, v);
    }

    for (size_t i = 0; i < sortedInterPts.size(); i++) {
        auto &group = sortedInterPts[i];

        if (i == 0 || i == sortedInterPts.size() - 1) {
            if (group.size() == 2) {
                group.pop_back();
            }
        } else if (group.size() == 1 && group.back().end != End::None) {
            group.push_back(group.back());
        }

        if (group.back().end == End::None) {
            continue;
        }

        vtkIdType id, beforeId, afterId;

        if (group.back().end == End::A) {
            id = group.back().edge.a;

            afterId = group.back().edge.b;

            auto currItr = std::find_if(verts.begin(), verts.end(), [&id](auto &v) { return v.id == id; });

            if (currItr == verts.begin()) {
                beforeId = (verts.end() - 1)->id;
            } else {
                beforeId = (currItr - 1)->id;
            }

        } else {
            id = group.back().edge.b;

            beforeId = group.back().edge.a;

            auto currItr = std::find_if(verts.begin(), verts.end(), [&id](auto &v) { return v.id == id; });

            if (currItr == verts.end() - 1) {
                afterId = verts.begin()->id;
            } else {
                afterId = (currItr + 1)->id;
            }
        }

        const auto &searchBefore = positions.find(beforeId);
        const auto &searchAfter = positions.find(afterId);

        if (group.size() == 2) {
            if (searchBefore == positions.end() && searchAfter != positions.end()) {
                if ((searchAfter->second < group.back().t && allVerts[beforeId].y < 0) ||
                    (searchAfter->second > group.back().t && allVerts[beforeId].y > 0)) {
                    group.pop_back();
                }
            } else if (searchBefore != positions.end() && searchAfter == positions.end()) {
                if ((searchBefore->second < group.back().t && allVerts[afterId].y > 0) ||
                    (searchBefore->second > group.back().t && allVerts[afterId].y < 0)) {
                    group.pop_back();
                }
            }
        }

        if (searchBefore == positions.end() && searchAfter == positions.end()) {
            const auto &a = allVerts[beforeId];
            const auto &b = allVerts[afterId];

            if (std::signbit(a.y) != std::signbit(b.y)) {
                if (group.size() == 2) {
                    group.pop_back();
                }
            } else {
                if ((a.x > b.x) == std::signbit(a.y)) {
                    group.clear();
                }
            }
        }
    }

    double q[3];

    InterPts cleanedInterPts;

    for (auto &group : sortedInterPts) {
        for (auto &p : group) {
            base.BackTransform(p.pt, q);
            p.pt[0] = q[0]; p.pt[1] = q[1]; p.pt[2] = q[2];
        }

        cleanedInterPts.insert(
            cleanedInterPts.end(),
            std::make_move_iterator(group.begin()),
            std::make_move_iterator(group.end())
        );
    }

    return cleanedInterPts;
}

typedef std::vector<std::tuple<InterPt, InterPt>> OverlapedLines;

void OverlapLines (OverlapedLines &overlaps, InterPts &intersA, InterPts &intersB) {

    auto Add = [](InterPt &a, InterPt &b, InterPt &c, InterPt &d) {
        a.Merge(c);
        b.Merge(d);

        return std::make_tuple(a, b);
    };

    InterPts::iterator itr, itr2;

    for (itr = intersA.begin(); itr != intersA.end(); itr += 2) {
        for (itr2 = intersB.begin(); itr2 != intersB.end(); itr2 += 2) {
            if (itr->t <= itr2->t && (itr+1)->t > itr2->t) {
                if ((itr2+1)->t < (itr+1)->t) {
                    overlaps.push_back(Add(*itr2, *(itr2+1), *itr, *(itr+1)));
                } else {
                    overlaps.push_back(Add(*itr2, *(itr+1), *itr, *(itr2+1)));
                }
            } else if (itr2->t <= itr->t && (itr2+1)->t > itr->t) {
                if ((itr+1)->t < (itr2+1)->t) {
                    overlaps.push_back(Add(*itr, *(itr+1), *itr2, *(itr2+1)));
                } else {
                    overlaps.push_back(Add(*itr, *(itr2+1), *itr2, *(itr+1)));
                }
            }
        }
    }

}

NonManifoldEdges GetNonManifoldEdges (const Mesh &mesh) {
    NonManifoldEdges edges;
    NonManifoldEdges nonManifoldEdges;

    for (const auto &face : mesh.faces) {
        size_t numVerts = face.vertIds.size();

        for (size_t i = 0; i < numVerts; i++) {
            vtkIdType u = face.vertIds[i];
            vtkIdType v = face.vertIds[(i + 1) % numVerts];

            if (!edges.emplace(u, v).second) {
                nonManifoldEdges.emplace(u, v);
                nonManifoldEdges.emplace(v, u);
            }
        }
    }

    return nonManifoldEdges;
}

vtkSmartPointer<vtkPolyData> GetPolyData (const Mesh &mesh) {
    auto pd = vtkSmartPointer<vtkPolyData>::New();
    auto points = vtkSmartPointer<vtkPoints>::New();

    vtkIdType numCells = static_cast<vtkIdType>(mesh.faces.size());

    auto cellIds = vtkSmartPointer<vtkIdTypeArray>::New();
    cellIds->SetName("OrigCellIds");

#if VTK_VERSION_NUMBER >= VTK_VERSION_CHECK(9, 7, 0)
    cellIds->ReserveTuples(numCells);
#else
    cellIds->Allocate(numCells);
#endif

    const vtkIdType numVerts = static_cast<vtkIdType>(mesh.verts.size());
    points->SetNumberOfPoints(numVerts);
    for (vtkIdType i = 0; i < numVerts; i++) {
        points->SetPoint(i, mesh.verts[i].x, mesh.verts[i].y, mesh.verts[i].z);
    }
    pd->SetPoints(points);

    auto polys = vtkSmartPointer<vtkCellArray>::New();

    for (const auto &face : mesh.faces) {
        polys->InsertNextCell(
            static_cast<vtkIdType>(face.vertIds.size()),
            face.vertIds.data()
        );

        cellIds->InsertNextValue(face.origCellId);
    }
    pd->SetPolys(polys);

    pd->GetCellData()->SetScalars(cellIds);

    return pd;
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

        double lenSq = normal.x * normal.x + normal.y * normal.y + normal.z * normal.z;

        if (lenSq > 1e-12) {
            double len = std::sqrt(lenSq);
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

Contact::Contact (vtkPolyData *origPdA, vtkPolyData *origPdB) {
    meshA = GetMesh(origPdA);
    meshB = GetMesh(origPdB);

    std::cout << meshA.verts.size() << std::endl;
    std::cout << meshB.verts.size() << std::endl;

    pdA = GetPolyData(meshA);
    pdB = GetPolyData(meshB);

    std::cout << ComputeMaxPlanarityError(pdA) << std::endl;
    std::cout << ComputeMaxPlanarityError(pdB) << std::endl;

    lines = CreateLines();
}

vtkSmartPointer<vtkPolyData> Contact::GetPolyDataA () {
    return pdA;
}

vtkSmartPointer<vtkPolyData> Contact::GetPolyDataB () {
    return pdB;
}

vtkSmartPointer<vtkPolyData> Contact::GetLines () {
    return lines;
}

vtkSmartPointer<vtkPolyData> Contact::CreateLines () {
    std::vector<Pair2> intersections;

    auto treeA = vtkSmartPointer<vtkOBBTree>::New();
    treeA->SetDataSet(pdA);

    auto treeB = vtkSmartPointer<vtkOBBTree>::New();
    treeB->SetDataSet(pdB);

    treeA->BuildLocator();
    treeB->BuildLocator();

    treeA->IntersectWithOBBTree(treeB, nullptr, Intersect, &intersections);

    UpdatePlanes(meshA);
    UpdatePlanes(meshB);

    std::vector<Pair2> congruentPairs;
    std::vector<Pair2> nonCongruentPairs;

    for (auto pair : intersections) {
        auto faceA = meshA.faces[pair.a];
        auto faceB = meshB.faces[pair.b];

        // std::cout << pair.a
        //     << ", " << faceA.normal
        //     << ", " << faceA.d
        //     << ", " << pair.b
        //     << ", " << faceB.normal
        //     << ", " << faceB.d
        //     << std::endl;

        double dotprod = faceA.normal.x * faceB.normal.x + faceA.normal.y * faceB.normal.y + faceA.normal.z * faceB.normal.z;

        if (std::abs(dotprod) > .9) {
            // std::cout << dotprod << std::endl;

            bool congruent{true};

            for (auto &id : faceA.vertIds) {
                double distance = faceB.normal.x * meshA.verts[id].x + faceB.normal.y * meshA.verts[id].y + faceB.normal.z * meshA.verts[id].z - faceB.d;

                if (std::abs(distance) > 1e-5) {
                    congruent = false;
                    break;
                }
            }

            if (congruent) {
                std::cout << "(" << pair.a << ", " << pair.b << ")" << std::endl;

                congruentPairs.push_back(pair);
            } else {
                nonCongruentPairs.push_back(pair);
            }
        }   else {
            nonCongruentPairs.push_back(pair);
        }
    }

    std::cout << congruentPairs.size() << std::endl;

    struct Node {
        char src;
        size_t index;

        bool operator< (const Node& other) const {
            if (src != other.src) {
                return src < other.src;
            }
            return index < other.index;
        }
    };

    std::map<Node, std::set<Node>> graph;

    for (auto pair : congruentPairs) {
        graph[Node{'A', static_cast<size_t>(pair.a)}].emplace(Node{'B', static_cast<size_t>(pair.b)});
        graph[Node{'B', static_cast<size_t>(pair.b)}].emplace(Node{'A', static_cast<size_t>(pair.a)});
    }

    const double maxProd = 0.9998476951563913;

    std::set<Node> visited;

    std::vector<std::set<Node>> groups;

    for (auto& [k, v] : graph) {
        if (visited.find(k) == visited.end()) {
            std::set<Node> group;

            std::stack<Node> stack;

            stack.push(k);
            visited.insert(k);

            const auto &seedNormal = k.src == 'A' ? meshA.faces[k.index].normal : meshB.faces[k.index].normal;

            while (!stack.empty()) {
                Node node = stack.top();
                stack.pop();

                group.insert(node);

                auto itr = graph.find(node);

                if (itr != graph.end()) {
                    for (const auto &neighbor : itr->second) {
                        if (visited.find(neighbor) == visited.end()) {
                            const auto &normal = neighbor.src == 'A' ? meshA.faces[neighbor.index].normal : meshB.faces[neighbor.index].normal;

                            const double prod = normal.x * seedNormal.x + normal.y * seedNormal.y + normal.z * seedNormal.z;

                            if (std::abs(prod) > maxProd) {
                                visited.insert(neighbor);
                                stack.push(neighbor);
                            }
                        }
                    }
                }
            }

            groups.push_back(group);
        }
    }

    auto planesA = vtkSmartPointer<vtkIdTypeArray>::New();
    auto planesB = vtkSmartPointer<vtkIdTypeArray>::New();

    planesA->SetNumberOfTuples(pdA->GetNumberOfPolys());
    planesB->SetNumberOfTuples(pdB->GetNumberOfPolys());

    for (vtkIdType i = 0; i < pdA->GetNumberOfPolys(); i++) {
        planesA->SetValue(i, -1);
    }

    for (vtkIdType i = 0; i < pdB->GetNumberOfPolys(); i++) {
        planesB->SetValue(i, -1);
    }

    std::vector<GlobalPlane> globalPlanes;

    for (const auto &group : groups) {
        for (const Node &node : group) {
            std::cout << "(" << node.src << node.index << "), ";
        }
        std::cout << std::endl;

        std::vector<size_t> faceIdsA, faceIdsB;

        for (const Node &node : group) {
            if (node.src == 'A') {
                faceIdsA.push_back(node.index);
            } else {
                faceIdsB.push_back(node.index);
            }
        }

        std::vector<Vector3> allPoints;

        for (auto &faceId : faceIdsA) {
            for (auto &id : meshA.faces[faceId].vertIds) {
                allPoints.emplace_back(meshA.verts[id]);
            }
        }

        for (auto &faceId : faceIdsB) {
            for (auto &id : meshB.faces[faceId].vertIds) {
                allPoints.emplace_back(meshB.verts[id]);
            }
        }

        std::cout << allPoints.size() << std::endl;

        Vector3 normal = meshA.faces[faceIdsA[0]].normal;
        double d = meshA.faces[faceIdsA[0]].d;

        std::cout << normal << ", " << d;

        FitPlane(allPoints, normal, d);

        std::cout << " -> " << normal << ", " << d << std::endl;

        size_t planeId = globalPlanes.size();

        globalPlanes.emplace_back(GlobalPlane{normal, d, faceIdsA, faceIdsB});

        for (auto &faceId : faceIdsA) {
            meshA.faces[faceId].globalPlaneId = planeId;

            planesA->SetValue(faceId, planeId);
        }

        for (auto &faceId : faceIdsB) {
            meshB.faces[faceId].globalPlaneId = planeId;

            planesB->SetValue(faceId, planeId);
        }

    }

    const int maxIters = 100;

    for (int iter = 0; iter < maxIters; iter++) {
        double dispA = FindNewVerts(meshA, globalPlanes);
        double dispB = FindNewVerts(meshB, globalPlanes);

        double maxDisp = std::max(dispA, dispB);

        UpdatePlanes(meshA);
        UpdatePlanes(meshB);

        std::for_each(std::execution::par, globalPlanes.begin(), globalPlanes.end(),
            [&](GlobalPlane &plane) {
                std::vector<Vector3> allPoints;

                for (auto &faceId : plane.faceIdsA) {
                    for (auto &id : meshA.faces[faceId].vertIds) {
                        allPoints.emplace_back(meshA.verts[id]);
                    }
                }

                for (auto &faceId : plane.faceIdsB) {
                    for (auto &id : meshB.faces[faceId].vertIds) {
                        allPoints.emplace_back(meshB.verts[id]);
                    }
                }

                FitPlane(allPoints, plane.normal, plane.d);
            });

        std::cout << iter << " maxDisp " << maxDisp << std::endl;

        if (maxDisp < 1e-6) {
            break;
        }

    }

    planesA->SetName("planes");
    planesB->SetName("planes");

    pdA->GetCellData()->AddArray(planesA);
    pdB->GetCellData()->AddArray(planesB);

    UpdatePoints(pdA, meshA);
    UpdatePoints(pdB, meshB);

    double maxDisp = 0;

    for (const auto &plane : globalPlanes) {
        for (auto &faceId : plane.faceIdsA) {
            for (auto &id : meshA.faces[faceId].vertIds) {
                const auto &p = meshA.verts[id];
                double dist = std::abs(plane.normal.x * p.x + plane.normal.y * p.y + plane.normal.z * p.z - plane.d);

                maxDisp = std::max(maxDisp, dist);
            }
        }

        for (auto &faceId : plane.faceIdsB) {
            for (auto &id : meshB.faces[faceId].vertIds) {
                const auto &p = meshB.verts[id];
                double dist = std::abs(plane.normal.x * p.x + plane.normal.y * p.y + plane.normal.z * p.z - plane.d);

                maxDisp = std::max(maxDisp, dist);
            }
        }
    }

    std::cout << "maxDisp " << maxDisp << std::endl;

    struct PairResult {
        vtkIdType pairA;
        vtkIdType pairB;
        InterPts interPtsA;
        InterPts interPtsB;
    };

    std::vector<size_t> indices(nonCongruentPairs.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::vector<PairResult> results(nonCongruentPairs.size());

    std::for_each(std::execution::par, indices.begin(), indices.end(),
        [&](size_t index) {
            const auto &pair = nonCongruentPairs[index];

            const auto &faceA = meshA.faces[pair.a];
            const auto &faceB = meshB.faces[pair.b];

            /*if (pair.a == 0 && pair.b == 2247) {
                std::cout << pair.a << " [";
                for (auto id : faceA.vertIds) {
                    std::cout << id << " " << meshA.verts[id] << ", ";
                }
                std::cout << "]" << std::endl;

                std::cout << pair.b << " [";
                for (auto id : faceB.vertIds) {
                    std::cout << id << " " << meshB.verts[id] << ", ";
                }
                std::cout << "]" << std::endl;
            } else {
                return;
            }*/

            const double nA[] = {faceA.normal.x, faceA.normal.y, faceA.normal.z};
            const double nB[] = {faceB.normal.x, faceB.normal.y, faceB.normal.z};

            double dets[3];

            vtkMath::Cross(nA, nB, dets);

            double r[] = {
                dets[0],
                dets[1],
                dets[2]
            };

            double s[3];

            vtkMath::Normalize(r);

            double dA = faceA.centroid.x*nA[0]+faceA.centroid.y*nA[1]+faceA.centroid.z*nA[2];
            double dB = faceB.centroid.x*nB[0]+faceB.centroid.y*nB[1]+faceB.centroid.z*nB[2];

            vtkMath::Cross(nA, nB, r);
            vtkMath::Normalize(r);

            int k = 0;
            if (std::abs(dets[1]) > std::abs(dets[k])) {
                k = 1;
            }

            if (std::abs(dets[2]) > std::abs(dets[k])) {
                k = 2;
            }

            if (std::abs(dets[k]) < 1e-8) {
                return;
            }

            // k=0 -> i=1, j=2
            // k=1 -> i=2, j=0
            // k=2 -> i=0, j=1

            int i = (k + 1) % 3;
            int j = (k + 2) % 3;

            s[i] = (dA*nB[j]-dB*nA[j])/dets[k];
            s[j] = (dB*nA[i]-dA*nB[i])/dets[k];
            s[k] = 0;

            Base baseA(s, r, nA);
            Base baseB(s, r, nB);

            auto transA = TransformFace2d(meshA, faceA, baseA);
            auto transB = TransformFace2d(meshB, faceB, baseB);

            auto interPtsA = IntersectFaceAndLine(baseA, transA, Src::A);
            auto interPtsB = IntersectFaceAndLine(baseB, transB, Src::B);

            // for (auto &p : interPtsA) {
            //     std::cout << "A " << p << std::endl;
            // }

            // for (auto &p : interPtsB) {
            //     std::cout << "B " << p << std::endl;
            // }

            results[index] = PairResult{pair.a, pair.b, std::move(interPtsA), std::move(interPtsB)};
        });

    auto pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetDataTypeToDouble();

    auto lines = vtkSmartPointer<vtkPolyData>::New();
    lines->SetPoints(pts);
    lines->Allocate(1000);

    auto contA = vtkSmartPointer<vtkIdTypeArray>::New();
    auto contB = vtkSmartPointer<vtkIdTypeArray>::New();

#if VTK_VERSION_NUMBER >= VTK_VERSION_CHECK(9, 7, 0)
    contA->ReserveTuples(1000);
    contB->ReserveTuples(1000);
#else
    contA->Allocate(1000);
    contB->Allocate(1000);
#endif

    contA->SetName("cA");
    contB->SetName("cB");

    lines->GetCellData()->AddArray(contA);
    lines->GetCellData()->AddArray(contB);

    auto nonManifoldEdgesA = GetNonManifoldEdges(meshA);
    auto nonManifoldEdgesB = GetNonManifoldEdges(meshB);

    auto AddLines = [&](InterPts &interA, InterPts &interB, vtkIdType idA, vtkIdType idB) {
        if (interA.size() == 0 || interB.size() == 0) {
            return;
        }

        OverlapedLines overlaps;
        OverlapLines(overlaps, interA, interB);

        OverlapedLines::const_iterator itr;

        for (itr = overlaps.begin(); itr != overlaps.end(); ++itr) {
            auto &f = std::get<0>(*itr);
            auto &s = std::get<1>(*itr);

            auto edgeF = std::make_pair(f.edge.a, f.edge.b);
            auto edgeS = std::make_pair(s.edge.a, s.edge.b);

            if ((f.src == Src::A && nonManifoldEdgesA.count(edgeF)) ||
                (s.src == Src::A && nonManifoldEdgesA.count(edgeS)) ||
                (f.src == Src::B && nonManifoldEdgesB.count(edgeF)) ||
                (s.src == Src::B && nonManifoldEdgesB.count(edgeS))) {

                throw "";
            }

            auto linePts = vtkSmartPointer<vtkIdList>::New();

            linePts->InsertNextId(pts->InsertNextPoint(f.pt));
            linePts->InsertNextId(pts->InsertNextPoint(s.pt));

            lines->InsertNextCell(VTK_LINE, linePts);

            contA->InsertNextValue(idA);
            contB->InsertNextValue(idB);
        }
    };

    try {
        for (auto &result : results) {
            if ((result.interPtsA.size() & 1) == 0 && (result.interPtsB.size() & 1) == 0) {
                AddLines(result.interPtsA, result.interPtsB, result.pairA, result.pairB);
            }
        }
    } catch (...) {
        std::cerr << "Contact failed" << std::endl;

        return nullptr;
    }

    auto clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInputData(lines);
    clean->ToleranceIsAbsoluteOn();
    clean->SetAbsoluteTolerance(1e-5);
    clean->ConvertLinesToPointsOff();
    clean->Update();

    return clean->GetOutput();
}

void Contact::Save (vtkPolyData *pd, const std::string &filename) {
    auto writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(pd);
    writer->Update();

    std::cout << filename << " saved" << std::endl;
}

void Contact::SaveAll (const std::string prefix) {
    if (pdA) {
        Save(pdA, prefix + "pdA.vtk");
    }

    if (pdB) {
        Save(pdB, prefix + "pdB.vtk");
    }

    if (lines) {
        Save(lines, prefix + "lines.vtk");
    }
}
