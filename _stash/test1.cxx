#include <array>
#include <vector>
#include <iostream>
#include <map>
#include <optional>
#include <string>

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkPolyDataWriter.h>
#include <vtkCleanPolyData.h>
#include <vtkCellArray.h>
#include <vtkMath.h>

std::vector<std::array<double, 2>> points{
    {1.9915899,22.038311}, {6.1931579,1.8044413}, {20.124676,13.74574}, // 0
    {20.124676,13.74574}, {6.1931579,1.8044413}, {13.549487,7.6912974},
    {13.549487,7.6912974}, {32.950517,4.3474943}, {20.124676,13.74574},
    {32.950517,4.3474943}, {13.549487,7.6912974}, {6.1931579,1.8044413}, // 3
    {32.950517,4.3474943}, {35.825276,21.706608}, {20.124676,13.74574},
    {34.742531,12.779188}, {35.825276,21.706608}, {32.950517,4.3474943},
    {34.117479,8.2812084}, {34.742531,12.779188}, {32.950517,4.3474943}, // 6
    {20.124676,13.74574}, {11.249987,18.062632}, {1.9915899,22.038311},
    {20.124676,13.74574}, {35.825276,21.706608}, {11.249987,18.062632},
    {35.825276,21.706608}, {23.119505,29.616721}, {11.249987,18.062632}, // 9
    {23.119505,29.616721}, {10.464885,38.22895}, {11.249987,18.062632},
    {11.249987,18.062632}, {10.464885,38.22895}, {1.9915899,22.038311},
    {23.119505,29.616721}, {35.825276,21.706608}, {10.464885,38.22895}, // 12
    {35.825276,21.706608}, {24.564057,29.409891}, {10.464885,38.22895},
    {24.564057,29.409891}, {17.940965,33.86098}, {10.464885,38.22895},
    {17.940965,33.86098}, {32.368538,37.184793}, {10.464885,38.22895}, // 15
    {17.940965,33.86098}, {24.564057,29.409891}, {32.368538,37.184793},
    {24.564057,29.409891}, {35.825276,21.706608}, {32.368538,37.184793}
};

struct Data {
    vtkIdType a, b, c, cellId;
    std::optional<vtkIdType> neighborId, replacedBy;

    bool resolved, atBoundary;

    friend std::ostream& operator<< (std::ostream &out, const Data &d) {
        out << "(a, b)=(" << d.a << ", " << d.b << "), c=" << d.c << ", cellId=" << d.cellId << ", neighborId=" << (d.neighborId ? std::to_string(*d.neighborId) : "null") << ", resolved=" << d.resolved;
        return out;
    }
};

void Detect (vtkPolyData *pd) {
    vtkCellArray *polys = pd->GetPolys();
    auto ids = vtkSmartPointer<vtkIdList>::New();

    vtkIdType a, b, c;

    double pA[3], pB[3], pC[3], vA[3], vB[3], vC[3], lA, lB, lC, max, epsilon;

    auto neigs = vtkSmartPointer<vtkIdList>::New();

    std::map<vtkIdType, Data> triangles;

    polys->InitTraversal();
    for (vtkIdType i = 0; i < polys->GetNumberOfCells(); i++) {
        polys->GetNextCell(ids);
        vtkIdType npts = ids->GetNumberOfIds();

        if (npts == 3) {
            a = ids->GetId(0);
            b = ids->GetId(1);
            c = ids->GetId(2);

            pd->GetPoint(a, pA);
            pd->GetPoint(b, pB);
            pd->GetPoint(c, pC);

            vtkMath::Subtract(pA, pB, vA);
            vtkMath::Subtract(pB, pC, vB);
            vtkMath::Subtract(pC, pA, vC);

            lA = vtkMath::Norm(vA);
            lB = vtkMath::Norm(vB);
            lC = vtkMath::Norm(vC);

            max = std::max({ lA, lB, lC });

            epsilon = ((lA + lB + lC) - 2 * max) / max;

            std::cout << i << " (" << epsilon << ")" << std::endl;

            if (epsilon < 1e-2) { // default sollte 1e-5 sein
                if (max == lA) {
                    pd->GetCellEdgeNeighbors(i, a, b, neigs);

                    if (neigs->GetNumberOfIds() == 1) {
                        triangles.emplace(i, Data{a, b, c, i, neigs->GetId(0), std::nullopt, false, false});
                    } else {
                        // am rand
                        triangles.emplace(i, Data{a, b, c, i, std::nullopt, std::nullopt, true, true});
                    }
                } else if (max == lB) {
                    pd->GetCellEdgeNeighbors(i, b, c, neigs);

                    if (neigs->GetNumberOfIds() == 1) {
                        triangles.emplace(i, Data{b, c, a, i, neigs->GetId(0), std::nullopt, false, false});
                    } else {
                        // am rand
                        triangles.emplace(i, Data{b, c, a, i, std::nullopt, std::nullopt, true, true});
                    }
                } else {
                    pd->GetCellEdgeNeighbors(i, c, a, neigs);

                    if (neigs->GetNumberOfIds() == 1) {
                        triangles.emplace(i, Data{c, a, b, i, neigs->GetId(0), std::nullopt, false, false});
                    } else {
                        // am rand
                        triangles.emplace(i, Data{c, a, b, i, std::nullopt, std::nullopt, true, true});
                    }
                }
            }
        }
    }

    std::map<vtkIdType, std::vector<vtkIdType>> replacedCells;

    vtkIdType npts;
    const vtkIdType *cell;

    bool startOver{false};

    do {
        startOver = false;

        for (auto& [cellId, data] : triangles) {
            std::cout << cellId << ", " << data << std::endl;

            if (data.resolved) {
                continue;
            }

            if (!data.neighborId.has_value()) {
                continue; // tritt aber nicht auf
            }

            if (triangles.count(*data.neighborId) == 0) {
                if (replacedCells.count(*data.neighborId) == 0) {
                    pd->GetCellPoints(*data.neighborId, npts, cell);

                    std::vector<vtkIdType> ids(cell, cell + npts);

                    replacedCells[*data.neighborId] = ids;
                }

                auto& ids = replacedCells[*data.neighborId];

                auto itr = std::find(ids.begin(), ids.end(), data.a);

                if (itr != ids.end()) {
                    ids.insert(itr, data.c);

                    data.resolved = true;
                    data.replacedBy = data.neighborId;

                    startOver = true;
                }
            } else {
                auto &neighbor = triangles[*data.neighborId];

                if (neighbor.a == data.b && neighbor.b == data.a) {
                    auto linkedCells = vtkSmartPointer<vtkIdList>::New();

                    pd->GetPointCells(neighbor.c, linkedCells);

                    for (vtkIdType i = 0; i < linkedCells->GetNumberOfIds(); i++) {
                        vtkIdType linkedCellId = linkedCells->GetId(i);

                        if (linkedCellId != neighbor.cellId) {
                            pd->ReplaceCellPoint(linkedCellId, neighbor.c, data.c);

                            if (triangles.count(linkedCellId) == 1) {
                                auto &linkedData = triangles[linkedCellId];

                                if (linkedData.a == neighbor.c) {
                                    linkedData.a = data.c;
                                    pd->GetCellEdgeNeighbors(linkedData.cellId, linkedData.a, linkedData.b, neigs);

                                    auto filtered = vtkSmartPointer<vtkIdList>::New();

                                    for (vtkIdType j = 0; j < neigs->GetNumberOfIds(); j++) {
                                        if (neigs->GetId(j) != cellId) {
                                            filtered->InsertNextId(neigs->GetId(j));
                                        }
                                    }

                                    if (filtered->GetNumberOfIds() == 1) {
                                        linkedData.neighborId = filtered->GetId(0);
                                    } else {
                                        // lange kante befindet sich jetzt am rand - es existiert kein nachbar, mit dem das dreieck aufgelöst werden kann

                                        linkedData.neighborId = std::nullopt;
                                        linkedData.resolved = linkedData.atBoundary = true;
                                    }
                                } else if (linkedData.b == neighbor.c) {
                                    linkedData.b = data.c;
                                    pd->GetCellEdgeNeighbors(linkedData.cellId, linkedData.a, linkedData.b, neigs);

                                    auto filtered = vtkSmartPointer<vtkIdList>::New();

                                    for (vtkIdType j = 0; j < neigs->GetNumberOfIds(); j++) {
                                        if (neigs->GetId(j) != cellId) {
                                            filtered->InsertNextId(neigs->GetId(j));
                                        }
                                    }

                                    if (filtered->GetNumberOfIds() == 1) {
                                        linkedData.neighborId = filtered->GetId(0);
                                    } else {
                                        // lange kante befindet sich jetzt am rand - es existiert kein nachbar, mit dem das dreieck aufgelöst werden kann

                                        linkedData.neighborId = std::nullopt;
                                        linkedData.resolved = linkedData.atBoundary = true;
                                    }
                                } else if (linkedData.c == neighbor.c) {
                                    linkedData.c = data.c;
                                }
                            }
                        }
                    }

                    data.resolved = neighbor.resolved = true;

                    pd->DeleteCell(cellId);
                    pd->DeleteCell(neighbor.cellId);

                    startOver = true;

                } else if (neighbor.resolved) {
                    if (neighbor.replacedBy) {
                        auto& ids = replacedCells[*neighbor.replacedBy];

                        auto itr = std::find(ids.begin(), ids.end(), data.a);

                        if (itr != ids.end()) {
                            ids.insert(itr, data.c);

                            data.resolved = true;
                            data.replacedBy = neighbor.replacedBy;

                            startOver = true;
                        }
                    } else if (neighbor.atBoundary) {
                        data.resolved = neighbor.atBoundary = true;
                    }

                }
            }
        }
    } while (startOver);

    for (auto& [cellId, data] : triangles) {
        if (data.resolved) {
            pd->DeleteCell(cellId);
        }
    }

    for (auto& [cellId, ids] : replacedCells) {
        pd->DeleteCell(cellId);

        auto pointIds = vtkSmartPointer<vtkIdList>::New();

        for (vtkIdType id : ids) {
            pointIds->InsertNextId(id);
        }

        pd->InsertNextCell(VTK_POLYGON, pointIds);
    }

    pd->RemoveDeletedCells();
}

int main () {
    auto pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetDataTypeToDouble();
    pts->SetNumberOfPoints(points.size());

    auto pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(pts);
    pd->Allocate(points.size()/3);

    vtkIdType i = 0;

    for (auto p : points) {
        std::cout << p[0] << ", " << p[1] << std::endl;

        pts->SetPoint(i++, p[0], p[1], 0);
    }

    pts->Print(std::cout);

    vtkIdType cell[3];

    for (i = 0; i < pts->GetNumberOfPoints()/3; i++) {
        cell[0] = i*3; cell[1] = i*3+1; cell[2] = i*3+2;
        pd->InsertNextCell(VTK_TRIANGLE, 3, cell);
    }

    auto clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInputData(pd);
    clean->Update();

    Detect(clean->GetOutput());

    auto writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(clean->GetOutput());
    writer->SetFileName("test_.vtk");
    writer->Update();

    return 0;
}

// hier müsste man nochmal non-manifold edges betrachten
