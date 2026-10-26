#include "Optimize2.h"

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

std::ostream& operator<< (std::ostream &out, const Vector3 &vec) {
    out << "[" << vec.x << ", " << vec.y << ", " << vec.z << "]";
    return out;
}

struct Pair {
    vtkIdType a, b;
};

std::vector<Pair> intersections;

int Intersect (vtkOBBNode *nodeA, vtkOBBNode *nodeB, vtkMatrix4x4 *, void *) {
    vtkIdList *cellsA = nodeA->Cells;
    vtkIdList *cellsB = nodeB->Cells;

    vtkIdType numCellsA = cellsA->GetNumberOfIds();
    vtkIdType numCellsB = cellsB->GetNumberOfIds();

    vtkIdType i, j, cellA, cellB;

    for (i = 0; i < numCellsA; i++) {
        cellA = cellsA->GetId(i);

        for (j = 0; j < numCellsB; j++) {
            cellB = cellsB->GetId(j);

            intersections.emplace_back(Pair{cellA, cellB});
        }
    }

    return 1;
}

struct GlobalPlane {
    Vector3 normal;
    double d;
    std::vector<size_t> faceIdsA, faceIdsB;
};

std::vector<GlobalPlane> globalPlanes;

struct Face {
    std::vector<vtkIdType> vertIds;
    Vector3 normal;
    Vector3 centroid;
    double d;
    std::optional<size_t> globalPlaneId;
};

struct Mesh {
    std::vector<Vector3> verts;
    std::vector<Face> faces;
    std::vector<std::vector<vtkIdType>> pointToFaces;
};

Mesh GetMesh (vtkPolyData *pd) {
    vtkPoints *points = pd->GetPoints();

     if (!points) {
        return {};
    }

    const vtkIdType numPoints = points->GetNumberOfPoints();
    const vtkIdType numFaces = pd->GetNumberOfPolys();

    if (numPoints == 0 || numFaces == 0) {
        return {};
    }

    std::vector<Face> faces(numFaces);
    std::vector<std::vector<vtkIdType>> pointToFaces(numPoints);

    vtkCellArray *polys = pd->GetPolys();
    auto ids = vtkSmartPointer<vtkIdList>::New();

    polys->InitTraversal();
    for (vtkIdType i = 0; i < numFaces; i++) {
        polys->GetNextCell(ids);
        vtkIdType npts = ids->GetNumberOfIds();
        faces[i].vertIds.reserve(npts);

        for (vtkIdType j = 0; j < npts; j++) {
            vtkIdType id = ids->GetId(j);
            faces[i].vertIds.push_back(id);
            pointToFaces[id].push_back(i);
        }
    }

    std::vector<Vector3> verts(numPoints);

    double p[3];

    for (vtkIdType i = 0; i < numPoints; i++) {
        points->GetPoint(i, p);
        verts[i] = {p[0], p[1], p[2]};
    }

    return { verts, faces, pointToFaces };
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

            double len = std::sqrt(face.normal.x * face.normal.x + face.normal.y * face.normal.y + face.normal.z * face.normal.z);

            face.normal.x /= len; face.normal.y /= len; face.normal.z /= len;

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

        double len = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

        if (len > 1e-12) {
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

double FindNewVerts (Mesh &mesh) {
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

void SaveMesh (vtkPolyData *pd, const Mesh &mesh, const std::string &filename) {
    vtkPoints *points = pd->GetPoints();

    for (size_t i = 0; i < mesh.verts.size(); i++) {
        points->SetPoint(static_cast<vtkIdType>(i), mesh.verts[i].x, mesh.verts[i].y, mesh.verts[i].z);
    }

    points->Modified();
    pd->Modified();

    auto writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(pd);
    writer->Write();

    std::cout << filename << " saved" << std::endl;
}

enum class Src {
    A,
    B
};

enum class End {
    None,
    A,
    B
};

enum class PointSrc {
    Calculated,
    Copied
};

template<typename T>
std::ostream& operator<< (typename std::enable_if<std::is_enum<T>::value, std::ostream>::type& stream, const T& e) {
    return stream << static_cast<typename std::underlying_type<T>::type>(e);
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

class InterPt {
public:
    InterPt () = delete;

    InterPt (double x, double y, double z, double t, vtkIdType a, vtkIdType b, End end, Src src, PointSrc pointSrc = PointSrc::Calculated) : t(t), edge{a, b}, end(end), src(src), srcA(std::nullopt), srcB(std::nullopt), pointSrc(pointSrc) {
        pt[0] = x;
        pt[1] = y;
        pt[2] = z;
    }

    double pt[3], t;
    Pair edge;
    End end;
    Src src;
    std::optional<vtkIdType> srcA, srcB;
    PointSrc pointSrc;

    friend std::ostream& operator<< (std::ostream &out, const InterPt &s) {
        out << "pt [" << s.pt[0] << ", " << s.pt[1] << ", " << s.pt[2] << "]"
            << ", t " << s.t
            << ", edge (" << s.edge.a << ", " << s.edge.b << ")"
            << ", end " << s.end
            << ", src " << s.src
            << ", pointSrc " << s.pointSrc;

        return out;
    }

    void Merge (const InterPt &other) {
        assert(src != other.src);

        if (src == Src::A) {
            srcA = GetEnd();
        } else {
            srcB = GetEnd();
        }

        if (std::abs(other.t-t) < 1e-5) {
            if (other.src == Src::A) {
                srcA = other.GetEnd();
            } else {
                srcB = other.GetEnd();
            }
        }
    }

    inline std::optional<vtkIdType> GetEnd () const {
        if (end == End::A) {
            return edge.a;
        }

        if (end == End::B) {
            return edge.b;
        }

        return std::nullopt;
    }

};

typedef std::vector<InterPt> InterPts;

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

    if (std::min(pA.y, pB.y) <= 0 && std::max(pA.y, pB.y) >= 0) {
        double div = pB.y - pA.y;

        if (std::abs(div) < 1e-12) {
            return;
        }

        double t = -pA.y / div;

        // std::cout << "t " << t << std::endl;

        if (t >= 0 && t <= 1) {
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

InterPts IntersectFaceAndLine (const Base &base, const std::vector<TransformedVert> &verts, const double *r, const double *s, Src src) {
    InterPts interPts;

    size_t numVerts = verts.size();

    for (size_t i = 0; i < numVerts; i++) {
        const auto &pA = verts[i];
        const auto &pB = verts[(i + 1) % numVerts];

        IntersectEdgeAndLine(interPts, pA, pB, src);
    }

    double q[3];

    for (auto &p : interPts) {
        base.BackTransform(p.pt, q);

        p.pt[0] = q[0]; p.pt[1] = q[1]; p.pt[2] = q[2];
    }

    // for (auto &p : interPts) {
    //     std::cout << src << " -> " << p << std::endl;
    // }

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

    std::map<InterPt, InterPts, Cmp> grouped;

    for (auto &p : interPts) {
        grouped[p].push_back(p);
    }

    std::vector<InterPts> sortedGroups;

    for (auto& [k, v] : grouped) {
        sortedGroups.push_back(v);
    }

    decltype(sortedGroups)::iterator itr;

    for (itr = sortedGroups.begin(); itr != sortedGroups.end(); ++itr) {
        if (itr == sortedGroups.begin()) {
            if (itr->size() == 2) {
                itr->pop_back();
            }
        } else if (itr == sortedGroups.end()-1) {
            if (itr->size() == 2) {
                itr->pop_back();
            }
        }

        // ...

    }

    InterPts cleanedInterPts;

    for (const auto &group : sortedGroups) {
        std::copy(group.begin(), group.end(), std::back_inserter(cleanedInterPts));
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

class Contact {
public:
    Contact () {
        pts = vtkSmartPointer<vtkPoints>::New();
        pts->SetDataTypeToDouble();

        lines = vtkSmartPointer<vtkPolyData>::New();
        lines->SetPoints(pts);
        lines->Allocate(1000);

        contA = vtkSmartPointer<vtkIdTypeArray>::New();
        contB = vtkSmartPointer<vtkIdTypeArray>::New();

        contA->ReserveTuples(1000);
        contB->ReserveTuples(1000);

        contA->SetName("cA");
        contB->SetName("cB");

        lines->GetCellData()->AddArray(contA);
        lines->GetCellData()->AddArray(contB);
    }

    void AddLines (InterPts &interA, InterPts &interB, vtkIdType idA, vtkIdType idB) {
        if (interA.size() == 0 || interB.size() == 0) {
            return;
        }

        OverlapedLines overlaps;
        OverlapLines(overlaps, interA, interB);

        OverlapedLines::const_iterator itr;

        for (itr = overlaps.begin(); itr != overlaps.end(); ++itr) {
            auto &f = std::get<0>(*itr);
            auto &s = std::get<1>(*itr);

            auto linePts = vtkSmartPointer<vtkIdList>::New();

            linePts->InsertNextId(pts->InsertNextPoint(f.pt));
            linePts->InsertNextId(pts->InsertNextPoint(s.pt));

            lines->InsertNextCell(VTK_LINE, linePts);

            contA->InsertNextValue(idA);
            contB->InsertNextValue(idB);
        }

    }

    void Save (const std::string &filename) {
        auto clean = vtkSmartPointer<vtkCleanPolyData>::New();
        clean->SetInputData(lines);
        clean->ToleranceIsAbsoluteOn();
        clean->SetAbsoluteTolerance(1e-5);
        clean->ConvertLinesToPointsOff();

        auto writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetFileName(filename.c_str());
        writer->SetInputConnection(clean->GetOutputPort());
        writer->Update();

        std::cout << filename << " saved" << std::endl;
    }

    vtkSmartPointer<vtkPoints> pts;
    vtkSmartPointer<vtkPolyData> lines;
    vtkSmartPointer<vtkIdTypeArray> contA, contB;
};


void Go (vtkPolyData *pdA, vtkPolyData *pdB) {
    std::cout << ComputeMaxPlanarityError(pdA) << std::endl;
    std::cout << ComputeMaxPlanarityError(pdB) << std::endl;

    auto treeA = vtkSmartPointer<vtkOBBTree>::New();
    treeA->SetDataSet(pdA);

    auto treeB = vtkSmartPointer<vtkOBBTree>::New();
    treeB->SetDataSet(pdB);

    treeA->BuildLocator();
    treeB->BuildLocator();

    treeA->IntersectWithOBBTree(treeB, nullptr, Intersect, nullptr);

    auto meshA = GetMesh(pdA);
    auto meshB = GetMesh(pdB);

    std::cout << meshA.verts.size() << std::endl;
    std::cout << meshB.verts.size() << std::endl;

    UpdatePlanes(meshA);
    UpdatePlanes(meshB);

    std::vector<Pair> congruentPairs;
    std::vector<Pair> nonCongruentPairs;

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

    auto groupsA = vtkSmartPointer<vtkIdTypeArray>::New();
    auto groupsB = vtkSmartPointer<vtkIdTypeArray>::New();

    groupsA->SetNumberOfTuples(pdA->GetNumberOfPolys());
    groupsB->SetNumberOfTuples(pdB->GetNumberOfPolys());

    for (vtkIdType i = 0; i < pdA->GetNumberOfPolys(); i++) {
        groupsA->SetValue(i, -1);
    }

    for (vtkIdType i = 0; i < pdB->GetNumberOfPolys(); i++) {
        groupsB->SetValue(i, -1);
    }

    vtkIdType groupId = 0;

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

            groupsA->SetValue(faceId, groupId);
        }

        for (auto &faceId : faceIdsB) {
            meshB.faces[faceId].globalPlaneId = planeId;

            groupsB->SetValue(faceId, groupId);
        }

        groupId++;

    }

    const int maxIters = 100;

    for (int iter = 0; iter < maxIters; iter++) {
        double dispA = FindNewVerts(meshA);
        double dispB = FindNewVerts(meshB);

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

    groupsA->SetName("groupsA");
    groupsB->SetName("groupsB");

    pdA->GetCellData()->AddArray(groupsA);
    pdB->GetCellData()->AddArray(groupsB);

    SaveMesh(pdA, meshA, "modified_MeshA.vtk");
    SaveMesh(pdB, meshB, "modified_MeshB.vtk");

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

            // if (pair.a == 18317 && pair.b == 1308) {
            //     std::cout << pair.a << " [";
            //     for (auto id : faceA.vertIds) {
            //         std::cout << id << ", ";
            //     }
            //     std::cout << std::endl;

            //     std::cout << pair.b << " [";
            //     for (auto id : faceB.vertIds) {
            //         std::cout << id << ", ";
            //     }
            //     std::cout << std::endl;
            // } else {
            //     continue;
            // }

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

            auto interPtsA = IntersectFaceAndLine(baseA, transA, r, s, Src::A);
            auto interPtsB = IntersectFaceAndLine(baseB, transB, r, s, Src::B);

            // for (auto &p : interPtsA) {
            //     std::cout << "A " << p << std::endl;
            // }

            // for (auto &p : interPtsB) {
            //     std::cout << "B " << p << std::endl;
            // }

            results[index] = PairResult{pair.a, pair.b, std::move(interPtsA), std::move(interPtsB)};
        });

    Contact contact;

    for (auto &result : results) {
        if ((result.interPtsA.size() & 1) == 0 && (result.interPtsB.size() & 1) == 0) {
            contact.AddLines(result.interPtsA, result.interPtsB, result.pairA, result.pairB);
        }
    }

    contact.Save("lines.vtk");

}

int main () {
    auto readerA = vtkSmartPointer<vtkPolyDataReader>::New();
    readerA->SetFileName("../bug3/MeshA.vtk");

    auto readerB = vtkSmartPointer<vtkPolyDataReader>::New();
    readerB->SetFileName("../bug3/MeshB.vtk");

    readerA->Update();
    readerB->Update();

    Go(readerA->GetOutput(), readerB->GetOutput());

    /*auto sphereA = vtkSmartPointer<vtkSphereSource>::New();
    auto sphereB = vtkSmartPointer<vtkSphereSource>::New();

    std::cout << "theta " << sphereA->GetThetaResolution() << ", phi " << sphereA->GetPhiResolution() << std::endl;

    const vtkIdType resolution = 100;

    sphereA->SetThetaResolution(resolution);
    sphereA->SetPhiResolution(resolution);

    sphereB->SetThetaResolution(resolution);
    sphereB->SetPhiResolution(resolution);

    sphereA->Update();
    sphereB->Update();

    auto transform = vtkSmartPointer<vtkTransform>::New();
    transform->RotateZ(180/resolution);
    transform->RotateX(180/resolution);

    auto tf = vtkSmartPointer<vtkTransformFilter>::New();
    tf->SetInputConnection(sphereB->GetOutputPort());
    tf->SetTransform(transform);

    tf->Update();

    Go(sphereA->GetOutput(), tf->GetPolyDataOutput());*/
}
