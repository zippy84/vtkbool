#include "Optimize2.h"

#include <iostream>
#include <vector>
#include <execution>
#include <map>
#include <set>
#include <stack>
#include <tuple>

#include <vtkSmartPointer.h>
#include <vtkOBBTree.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkIdList.h>
#include <vtkMath.h>

std::ostream& operator<< (std::ostream &out, const Vector3 &vec) {
    out << "[" << vec.x << ", " << vec.y << ", " << vec.z << "]";
    return out;
}

struct Pair {
    vtkIdType a, b;
};

struct Plane2 {
    vtkIdType cellId;
    Vector3 normal;
    Vector3 centroid;
    double d;
    double R;
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

std::vector<Plane2> GetPlanes (vtkPolyData *pd) {
    const vtkIdType numPolys = pd->GetNumberOfPolys();

    std::vector<Plane2> planes(numPolys);

    vtkCellArray *polys = pd->GetPolys();

    vtkIdType npts;
    const vtkIdType *pts;

    vtkIdType cellId = 0;

    double pA[3], pB[3];

    polys->InitTraversal();
    while (polys->GetNextCell(npts, pts)) {
        Vector3 centroid{0, 0, 0};
        Vector3 normal{0, 0, 0};

        for (vtkIdType i = 0; i < npts; i++) {
            pd->GetPoint(pts[i], pA);
            pd->GetPoint(pts[(i + 1) % npts], pB);

            centroid.x += pA[0];
            centroid.y += pA[1];
            centroid.z += pA[2];

            normal.x += (pA[1] - pB[1]) * (pA[2] + pB[2]);
            normal.y += (pA[2] - pB[2]) * (pA[0] + pB[0]);
            normal.z += (pA[0] - pB[0]) * (pA[1] + pB[1]);
        }

        centroid.x /= static_cast<double>(npts);
        centroid.y /= static_cast<double>(npts);
        centroid.z /= static_cast<double>(npts);

        double len = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

        normal.x /= len; normal.y /= len; normal.z /= len;

        double R = 0;

        double p[3];

        for (vtkIdType i = 0; i < npts; i++) {
            pd->GetPoint(pts[i], p);

            const double dx = p[0] - centroid.x;
            const double dy = p[1] - centroid.y;
            const double dz = p[2] - centroid.z;

            R = std::max(R, dx * dx + dy * dy + dz * dz);
        }

        double d = normal.x * centroid.x + normal.y * centroid.y + normal.z * centroid.z;

        planes[cellId] = Plane2{cellId, normal, centroid, d, R};

        cellId++;
    }

    return planes;
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

class Plane3 {
public:
    Plane3 (const vtkIdType &cellId, const Vector3 &normal, const Vector3 &centroid, const double &d) : cellId(cellId), normal(normal), centroid(centroid), d(d) {}
    Plane3 (const Plane2 &plane) : cellId(plane.cellId), normal(plane.normal), centroid(plane.centroid), d(plane.d) {}

    vtkIdType cellId;
    Vector3 normal;
    Vector3 centroid;
    double d;
    std::vector<std::tuple<double, double, double, vtkIdType>> points;
};

void ExtendPlanes (std::vector<Plane3> &planes, vtkPolyData *pd) {
    vtkIdType npts;
    const vtkIdType *pts;

    double p[3];

    for (auto &plane : planes) {
        pd->GetCellPoints(plane.cellId, npts, pts);

        for (vtkIdType i = 0; i < npts; i++) {
            pd->GetPoint(pts[i], p);
            plane.points.emplace_back(p[0], p[1], p[2], pts[i]);
        }
    }
}

int main () {
    auto readerA = vtkSmartPointer<vtkPolyDataReader>::New();
    readerA->SetFileName("../bug3/MeshA.vtk");

    auto readerB = vtkSmartPointer<vtkPolyDataReader>::New();
    readerB->SetFileName("../bug3/MeshB.vtk");

    readerA->Update();
    readerB->Update();

    auto pdA = readerA->GetOutput();
    auto pdB = readerB->GetOutput();

    std::cout << ComputeMaxPlanarityError(pdA) << std::endl;
    std::cout << ComputeMaxPlanarityError(pdB) << std::endl;

    auto treeA = vtkSmartPointer<vtkOBBTree>::New();
    treeA->SetDataSet(pdA);

    auto treeB = vtkSmartPointer<vtkOBBTree>::New();
    treeB->SetDataSet(pdB);

    treeA->BuildLocator();
    treeB->BuildLocator();

    treeA->IntersectWithOBBTree(treeB, nullptr, Intersect, nullptr);

    std::for_each(std::execution::par, intersections.begin(), intersections.end(),
        [&](const Pair &pair) {
            // ...
        });

    auto planesA = GetPlanes(pdA);
    auto planesB = GetPlanes(pdB);

    std::vector<Pair> congruentPairs;

    for (auto pair : intersections) {
        auto planeA = planesA[pair.a];
        auto planeB = planesB[pair.b];

        // std::cout << pair.a
        //     << ", " << planeA.normal
        //     << ", " << planeA.d
        //     << ", " << pair.b
        //     << ", " << planeB.normal
        //     << ", " << planeB.d
        //     << std::endl;

        double dotprod = planeA.normal.x * planeB.normal.x + planeA.normal.y * planeB.normal.y + planeA.normal.z * planeB.normal.z;

        if (std::abs(dotprod) > .9) {
            // std::cout << dotprod << std::endl;

            vtkIdType npts;
            const vtkIdType *pts;

            double p[3];

            pdA->GetCellPoints(planeA.cellId, npts, pts);

            bool congruent{true};

            for (vtkIdType i = 0; i < npts; i++) {
                pdA->GetPoint(pts[i], p);

                double distance = planeB.normal.x * p[0] + planeB.normal.y * p[1] + planeB.normal.z * p[2] - planeB.d;

                if (std::abs(distance) > 1e-5) {
                    congruent = false;
                    break;
                }
            }

            if (congruent) {
                // std::cout << "(" << pair.a << ", " << pair.b << ")" << std::endl;

                congruentPairs.push_back(pair);
            }
        }
    }

    struct Node {
        char src;
        vtkIdType index;

        bool operator< (const Node& other) const {
            if (src != other.src) {
                return src < other.src;
            }
            return index < other.index;
        }
    };

    std::map<Node, std::set<Node>> graph;

    for (auto pair : congruentPairs) {
        graph[Node{'A', pair.a}].emplace(Node{'B', pair.b});
        graph[Node{'B', pair.b}].emplace(Node{'A', pair.a});
    }

    std::set<Node> visited;

    std::vector<std::set<Node>> groups;

    for (auto& [k, v] : graph) {
        if (visited.find(k) == visited.end()) {
            std::set<Node> group;

            std::stack<Node> stack;

            stack.push(k);
            visited.insert(k);

            while (!stack.empty()) {
                Node node = stack.top();
                stack.pop();

                group.insert(node);

                auto itr = graph.find(node);

                if (itr != graph.end()) {
                    for (const auto &neighbor : itr->second) {
                        if (visited.find(neighbor) == visited.end()) {
                            visited.insert(neighbor);
                            stack.push(neighbor);
                        }
                    }
                }
            }

            groups.push_back(group);
        }
    }

    for (const auto &group : groups) {
        for (const Node &node : group) {
            std::cout << "(" << node.src << node.index << "), ";
        }
        std::cout << std::endl;

        std::vector<Plane3> resolvedA, resolvedB;

        for (const Node &node : group) {
            if (node.src == 'A') {
                resolvedA.push_back(planesA[node.index]);
            } else {
                resolvedB.push_back(planesB[node.index]);
            }
        }

        std::vector<Vector3> allPoints;

        ExtendPlanes(resolvedA, pdA);
        ExtendPlanes(resolvedB, pdB);

        for (auto &plane : resolvedA) {
            for (auto &point : plane.points) {
                allPoints.emplace_back(Vector3{std::get<0>(point), std::get<1>(point), std::get<2>(point)});
            }
        }

        for (auto &plane : resolvedB) {
            for (auto &point : plane.points) {
                allPoints.emplace_back(Vector3{std::get<0>(point), std::get<1>(point), std::get<2>(point)});
            }
        }

        std::cout << allPoints.size() << std::endl;

        Vector3 normal = resolvedA[0].normal;
        double d = resolvedA[0].d;

        std::cout << normal << ", " << d;

        FitPlane(allPoints, normal, d);

        std::cout << " -> " << normal << ", " << d << std::endl;

        for (auto &plane : resolvedA) {
            for (auto &pt : plane.points) {
                double distance = normal.x * std::get<0>(pt) + normal.y * std::get<1>(pt) + normal.z * std::get<2>(pt) - d;

                double x = std::get<0>(pt) + normal.x * distance;
                double y = std::get<1>(pt) + normal.y * distance;
                double z = std::get<2>(pt) + normal.z * distance;

                std::cout << std::fixed << std::setprecision(6)
                    << "A " << std::get<3>(pt)
                    << ": [" << std::get<0>(pt)
                    << ", " << std::get<1>(pt)
                    << ", " << std::get<2>(pt)
                    << "] -> ["
                    << x
                    << ", " << y
                    << ", " << z
                    << "]"
                    << std::endl;
            }
        }

        for (auto &plane : resolvedB) {
            for (auto &pt : plane.points) {
                double distance = normal.x * std::get<0>(pt) + normal.y * std::get<1>(pt) + normal.z * std::get<2>(pt) - d;

                double x = std::get<0>(pt) + normal.x * distance;
                double y = std::get<1>(pt) + normal.y * distance;
                double z = std::get<2>(pt) + normal.z * distance;

                std::cout << std::fixed << std::setprecision(6)
                    << "B " << std::get<3>(pt)
                    << ": [" << std::get<0>(pt)
                    << ", " << std::get<1>(pt)
                    << ", " << std::get<2>(pt)
                    << "] -> ["
                    << x
                    << ", " << y
                    << ", " << z
                    << "]"
                    << std::endl;
            }
        }

    }
}
