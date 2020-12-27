/*
Copyright 2012-2020 Ronald Römer

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

#include <vtkKdTreePointLocator.h>
#include <vtkCubeSource.h>
#include <vtkCylinderSource.h>
#include <vtkSphereSource.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkMath.h>
#include <vtkTrivialProducer.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkTriangleFilter.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkPlaneSource.h>
#include <vtkAppendPolyData.h>
#include <vtkLineSource.h>
#include <vtkTubeFilter.h>
#include <vtkCommand.h>
#include <vtkPolyDataConnectivityFilter.h>

#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <array>

#include "vtkPolyDataBooleanFilter.h"

#include "Utilities.h"

typedef std::map<vtkIdType, IdsType> LinksType;

class Test {
    vtkPolyData *pd, *lines;
    vtkKdTreePointLocator *loc;
    vtkIdList *cells, *poly, *pts;

    vtkIdTypeArray *contsA, *contsB;

public:
    Test (vtkPolyData *pd, vtkPolyData *lines) : pd(pd), lines(lines) {
        loc = vtkKdTreePointLocator::New();
        loc->SetDataSet(pd);
        loc->BuildLocator();

        lines->BuildLinks();
        pd->BuildLinks();

        cells = vtkIdList::New();
        poly = vtkIdList::New();
        pts = vtkIdList::New();

        contsA = vtkIdTypeArray::SafeDownCast(lines->GetCellData()->GetArray("cA"));
        contsB = vtkIdTypeArray::SafeDownCast(lines->GetCellData()->GetArray("cB"));
    }
    ~Test () {
        loc->FreeSearchStructure();
        loc->Delete();
        cells->Delete();
        poly->Delete();
        pts->Delete();
    }

    int run () {

        vtkIdTypeArray *origCellIdsA = vtkIdTypeArray::SafeDownCast(pd->GetCellData()->GetArray("OrigCellIdsA"));
        vtkIdTypeArray *origCellIdsB = vtkIdTypeArray::SafeDownCast(pd->GetCellData()->GetArray("OrigCellIdsB"));

        int errA = checkConnectivity(0, origCellIdsA);
        int errB = checkConnectivity(1, origCellIdsB);

        return errA > 0 || errB > 0;

    }

private:
    int checkConnectivity (int input, vtkIdTypeArray *origIds) {
        std::cout << "Checking input " << input << std::endl;

        vtkIdType numCells = pd->GetNumberOfCells(),
            numLines = lines->GetNumberOfCells();

        vtkIdType i, j, k, l, m;

        int err = 0;

        for (i = 0; i < numCells && err < 1; i++) {
            if (origIds->GetValue(i) != NO_USE) {
                pd->GetCellPoints(i, poly);
                if (poly->GetNumberOfIds() < 3) {
                    std::cout << "poly " << i << " has too few points" << std::endl;

                    err++;
                }
            }
        }

        if (err == 0) {
            double ptA[3], ptB[3];

            vtkIdList *line = vtkIdList::New();

            for (i = 0; i < numLines; i++) {

                lines->GetCellPoints(i, line);

                lines->GetPoint(line->GetId(0), ptA);
                lines->GetPoint(line->GetId(1), ptB);

                FindPoints(loc, ptA, pts);

                LinksType links;

                // sammelt polygone und deren punkte am jeweiligen ende der linie

                // bei einem normalen schnitt sollte jedes polygon zwei punkte haben

                for (j = 0; j < pts->GetNumberOfIds(); j++) {
                    pd->GetPointCells(pts->GetId(j), cells);

                    for (k = 0; k < cells->GetNumberOfIds(); k++) {
                        if (origIds->GetValue(cells->GetId(k)) != NO_USE) {
                            links[cells->GetId(k)].push_back(pts->GetId(j));
                        }
                    }
                }

                FindPoints(loc, ptB, pts);

                for (j = 0; j < pts->GetNumberOfIds(); j++) {
                    pd->GetPointCells(pts->GetId(j), cells);

                    for (k = 0; k < cells->GetNumberOfIds(); k++) {
                        if (origIds->GetValue(cells->GetId(k)) != NO_USE) {
                            links[cells->GetId(k)].push_back(pts->GetId(j));
                        }
                    }
                }

                LinksType::const_iterator itr;

                IdsType polys;

                for (itr = links.begin(); itr != links.end(); ++itr) {
                    const IdsType &pts = itr->second;

                    if (pts.size() > 1) {
                        // gibt es doppelte punkte, die sich direkt nebeneinander befinden?

                        pd->GetCellPoints(itr->first, poly);
                        vtkIdType numPts = poly->GetNumberOfIds();

                        std::map<vtkIdType, Point3d> _pts;

                        vtkIdType _id, idA, idB;

                        double pt[3];

                        for (j = 0; j < numPts; j++) {
                            _id = poly->GetId(j);
                            pd->GetPoint(_id, pt);

                            _pts.emplace(_id, Point3d{pt[0], pt[1], pt[2]});
                        }

                        for (j = 0; j < numPts; j++) {
                            idA = poly->GetId(j);
                            idB = poly->GetId(j+1 == numPts ? 0 : j+1);

                            auto _a = _pts.find(idA);
                            auto _b = _pts.find(idB);

                            if (_a->second == _b->second) {
                                std::cout << "poly " << itr->first << " has duplicated points" << std::endl;

                                err++;
                                break;
                            }
                        }

                        polys.push_back(itr->first);
                    }
                }

                if (err > 0) {
                    break;
                }

                if (polys.size() == 2) {
                    const IdsType &f = links[polys[0]],
                        &s = links[polys[1]];

                    lines->GetPointCells(line->GetId(0), cells);
                    vtkIdType linkedA = cells->GetNumberOfIds();

                    lines->GetPointCells(line->GetId(1), cells);
                    vtkIdType linkedB = cells->GetNumberOfIds();

                    if (linkedA == linkedB) {
                        std::set<vtkIdType> all;

                        all.insert(f.begin(), f.end());
                        all.insert(s.begin(), s.end());

                        if (all.size() < 4) {
                            err++;
                        }
                    } else if (linkedA == 4) {
                        if (f[0] == s[0]) {
                            err++;
                        }
                    } else {
                        if (f[1] == s[1]) {
                            err++;
                        }
                    }

                    if (err > 0) {
                        std::cout << "polys " << polys[0] << ", " << polys[1]
                            << " are connected" << std::endl;
                    }

                } else if (polys.size() == 1) {
                    const IdsType &f = links[polys[0]];

                    if (f.size() < 3) {
                        std::cout << "line " << i << " has too few neighbors" << std::endl;
                        err++;
                    } else {
                        if (f[0] == f[1] || f[0] == f[2] || f[1] == f[2]) {
                            std::cout << "poly " << polys[0] << " is self connected" << std::endl;
                            err++;
                        }
                    }
                } else if (polys.size() == 0) {
                    std::cout << "line " << i << " has zero neighbors" << std::endl;
                    err++;
                }

                // sucht nach lücken

                for (j = 0; j < 2 && err < 1; j++) {
                    double pt[3];
                    lines->GetPoint(line->GetId(j), pt);

                    IdsType ids;

                    FindPoints(loc, pt, pts);
                    vtkIdType numPtsA = pts->GetNumberOfIds();

                    std::map<vtkIdType, vtkIdType> masked;

                    for (k = 0; k < numPtsA; k++) {

                        pd->GetPointCells(pts->GetId(k), cells);
                        vtkIdType _numCells = cells->GetNumberOfIds();

                        for (l = 0; l < _numCells; l++) {
                            vtkIdType cellId = cells->GetId(l);

                            if (origIds->GetValue(cellId) != NO_USE) {
                                pd->GetCellPoints(cells->GetId(l), poly);
                                vtkIdType numPtsB = poly->GetNumberOfIds();

                                for (m = 0; m < numPtsB; m++) {
                                    vtkIdType p = poly->GetId(m);

                                    /* wenn ein polygon mehrfach die gleiche p verwendet,
                                       dann erscheint die cellId auch mehrfach in cells
                                       => es sind dann auch 4 kanten zu verarbeiten
                                       => der betrachtete relative index wird dann zusammen mit der cellId in masked gespeichert */

                                    if (masked.count(cellId) == 1 && masked[cellId] == m) {
                                        continue;
                                    }

                                    if (p == pts->GetId(k)) {

                                        vtkIdType idA = m == numPtsB-1 ? 0 : m+1,
                                            idB = m == 0 ? numPtsB-1 : m-1;

                                        ids.push_back(poly->GetId(idA));
                                        ids.push_back(poly->GetId(idB));

                                        masked[cellId] = m;

                                        break;
                                    }
                                }

                            }
                        }

                    }

                    std::sort(ids.begin(), ids.end());

                    IdsType ids2 = ids,
                        ids3;

                    ids2.erase(std::unique(ids2.begin(), ids2.end()), ids2.end());

                    for (vtkIdType id : ids2) {
                        if (std::count(ids.begin(), ids.end(), id) == 1) {
                            ids3.push_back(id);
                        }
                    }

                    std::set<Point3d> uniquePts;

                    double _pt[3];

                    std::size_t count {0};

                    for (vtkIdType id : ids3) {
                        pd->GetPoint(id, _pt);

                        auto _ins = uniquePts.emplace(_pt[0], _pt[1], _pt[2]);

                        if (!_ins.second) {
                            count++;
                        }
                    }

                    if (count*2 != ids3.size()) {
                        std::cout << "line " << i << " has gaps at " << line->GetId(j) << std::endl;
                        err++;
                    }

                }

            }

            line->Delete();

        }

        return err;

    }

};

class Observer : public vtkCommand {
public:
    bool hasError;
    std::string msg;

    Observer () : hasError(false) {}

    static Observer *New () {
        return new Observer;
    }

    virtual void Execute (vtkObject *vtkNotUsed(caller), unsigned long event, void *calldata) {
        hasError = event == vtkCommand::ErrorEvent;
        msg = static_cast<char*>(calldata);
    }

    void Clear() {
        hasError = false;
        msg.clear();
    }
};

int main (int vtkNotUsed(argc), char *argv[]) {
    std::istringstream stream(argv[1]);
    int t;

    stream >> t;

    if (t == 0) {
        vtkCubeSource *cu = vtkCubeSource::New();
        cu->SetYLength(.5);

        vtkCylinderSource *cyl = vtkCylinderSource::New();
        cyl->SetResolution(32);
        cyl->SetHeight(.5);
        cyl->SetCenter(0, .5, 0);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, cu->GetOutputPort());
        bf->SetInputConnection(1, cyl->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        cyl->Delete();
        cu->Delete();

        return ok;

    } else if (t == 1) {
        vtkCubeSource *cu = vtkCubeSource::New();
        cu->SetYLength(.5);

        vtkCylinderSource *cyl = vtkCylinderSource::New();
        cyl->SetResolution(32);
        cyl->SetHeight(.5);
        cyl->SetCenter(0, .25, 0);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, cu->GetOutputPort());
        bf->SetInputConnection(1, cyl->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        cyl->Delete();
        cu->Delete();

        return ok;

    } else if (t == 2) {
        vtkCubeSource *cu = vtkCubeSource::New();
        cu->SetYLength(.5);

        vtkCylinderSource *cyl = vtkCylinderSource::New();
        cyl->SetResolution(32);
        cyl->SetHeight(.5);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, cu->GetOutputPort());
        bf->SetInputConnection(1, cyl->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        cyl->Delete();
        cu->Delete();

        return ok;

    } else if (t == 3) {
        vtkCubeSource *cubeA = vtkCubeSource::New();

        vtkCubeSource *cubeB = vtkCubeSource::New();
        cubeB->SetXLength(.5);
        cubeB->SetZLength(.5);
        cubeB->SetCenter(0, 0, .25);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, cubeA->GetOutputPort());
        bf->SetInputConnection(1, cubeB->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        cubeB->Delete();
        cubeA->Delete();

        return ok;

    } else if (t == 4) {
        vtkCubeSource *cubeA = vtkCubeSource::New();

        vtkCubeSource *cubeB = vtkCubeSource::New();
        cubeB->SetCenter(0, 0, .75);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, cubeA->GetOutputPort());
        bf->SetInputConnection(1, cubeB->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        cubeB->Delete();
        cubeA->Delete();

        return ok;

    } else if (t == 5) {
        // mehrer punkte in einem strip ohne fläche

        vtkCubeSource *cu = vtkCubeSource::New();

        vtkPolyData *cyl = vtkPolyData::New();
        cyl->Allocate(130, 1);

        vtkPoints *pts = vtkPoints::New();
        pts->SetNumberOfPoints(576);

        vtkIdList *top = vtkIdList::New();
        top->SetNumberOfIds(32);

        vtkIdList *bottom = vtkIdList::New();
        bottom->SetNumberOfIds(32);

        vtkIdList *poly = vtkIdList::New();
        poly->SetNumberOfIds(4);

        int ind = 0;

        for (int i = 0; i < 32; i++) {
            double x0 = .5*std::cos(i*2*M_PI/32);
            double z0 = .5*std::sin(i*2*M_PI/32);

            double x1 = .5*std::cos((i+1)*2*M_PI/32);
            double z1 = .5*std::sin((i+1)*2*M_PI/32);

            pts->SetPoint(ind, x0, .75, z0);
            top->SetId(i, ind++);

            pts->SetPoint(ind, x0, -.25, z0);
            bottom->SetId(i, ind++);

            for (int j = 0; j < 4; j++) {
                pts->SetPoint(ind, x0, -.25+j/4., z0);
                poly->SetId(0, ind++);
                pts->SetPoint(ind, x0, -.25+(j+1)/4., z0);
                poly->SetId(1, ind++);
                pts->SetPoint(ind, x1, -.25+(j+1)/4., z1);
                poly->SetId(2, ind++);
                pts->SetPoint(ind, x1, -.25+j/4., z1);
                poly->SetId(3, ind++);

                cyl->InsertNextCell(VTK_QUAD, poly);
            }
        }

        cyl->ReverseCell(cyl->InsertNextCell(VTK_POLYGON, top));
        cyl->InsertNextCell(VTK_POLYGON, bottom);

        poly->Delete();
        bottom->Delete();
        top->Delete();

        cyl->SetPoints(pts);

        vtkTrivialProducer *prod = vtkTrivialProducer::New();
        prod->SetOutput(cyl);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, cu->GetOutputPort());
        bf->SetInputConnection(1, prod->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        prod->Delete();

        pts->Delete();
        cyl->Delete();
        cu->Delete();

        return ok;

    } else if (t == 6) {
        vtkCubeSource *cu = vtkCubeSource::New();

        vtkSphereSource *sp = vtkSphereSource::New();
        sp->SetRadius(.5);
        sp->SetCenter(.5, .5, .5);
        sp->SetPhiResolution(100);
        sp->SetThetaResolution(100);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, cu->GetOutputPort());
        bf->SetInputConnection(1, sp->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        sp->Delete();
        cu->Delete();

        return ok;

    }  else if (t == 7) {
        // enthält schmale schnitte an den ecken der polygone

        vtkSphereSource *spA = vtkSphereSource::New();
        spA->SetRadius(50);
        spA->SetPhiResolution(6);
        spA->SetThetaResolution(6);


        vtkSphereSource *spB = vtkSphereSource::New();
        spB->SetRadius(50);
        spB->SetCenter(0, 25, 0);
        spB->SetPhiResolution(10);
        spB->SetThetaResolution(600);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, spA->GetOutputPort());
        bf->SetInputConnection(1, spB->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        spB->Delete();
        spA->Delete();

        return ok;

    } else if (t == 8) {
        vtkSphereSource *spA = vtkSphereSource::New();
        spA->SetRadius(50);
        spA->SetPhiResolution(97);
        spA->SetThetaResolution(100);

        vtkSphereSource *spB = vtkSphereSource::New();
        spB->SetRadius(50);
        spB->SetCenter(25, 0, 0);
        spB->SetPhiResolution(81);
        spB->SetThetaResolution(195);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, spA->GetOutputPort());
        bf->SetInputConnection(1, spB->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        spB->Delete();
        spA->Delete();

        return ok;

    } else if (t == 9) {
        // enthält sehr scharfwinklige schnitte mit kanten

        vtkSphereSource *spA = vtkSphereSource::New();
        spA->SetRadius(50);
        spA->SetPhiResolution(100);
        spA->SetThetaResolution(100);

        vtkSphereSource *spB = vtkSphereSource::New();
        spB->SetRadius(50);
        spB->SetCenter(25, 0, 0);
        spB->SetPhiResolution(251);
        spB->SetThetaResolution(251);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, spA->GetOutputPort());
        bf->SetInputConnection(1, spB->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        spB->Delete();
        spA->Delete();

        return ok;

    } else if (t == 10) {
        // strip liegt auf einer kante und beginnt im gleichen punkt

        vtkCubeSource *cu = vtkCubeSource::New();

        vtkCylinderSource *cyl = vtkCylinderSource::New();
        cyl->SetResolution(32);
        cyl->SetRadius(.25);
        cyl->SetCenter(.25, 0, 0);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, cu->GetOutputPort());
        bf->SetInputConnection(1, cyl->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        cyl->Delete();
        cu->Delete();

        return ok;

    } else if (t == 11) {
        vtkCubeSource *cubeA = vtkCubeSource::New();

        vtkTriangleFilter *tf = vtkTriangleFilter::New();
        tf->SetInputConnection(cubeA->GetOutputPort());

        vtkLinearSubdivisionFilter *sf = vtkLinearSubdivisionFilter::New();
        sf->SetInputConnection(tf->GetOutputPort());
        sf->SetNumberOfSubdivisions(4);

        vtkCubeSource *cubeB = vtkCubeSource::New();
        cubeB->SetXLength(.5);
        cubeB->SetZLength(.5);
        cubeB->SetCenter(0, 0, .25);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, sf->GetOutputPort());
        bf->SetInputConnection(1, cubeB->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        sf->Delete();
        tf->Delete();
        cubeB->Delete();
        cubeA->Delete();

        return ok;

    } else if (t == 12) {
        vtkSphereSource *sp = vtkSphereSource::New();
        sp->SetRadius(.5);
        sp->SetPhiResolution(8);
        sp->SetThetaResolution(8);

        vtkPlaneSource *pl = vtkPlaneSource::New();
        pl->SetResolution(5, 5);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, sp->GetOutputPort());
        bf->SetInputConnection(1, pl->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        pl->Delete();
        sp->Delete();

        return ok;

    } else if (t == 13) {
        // testet den merger

        vtkCubeSource *cuA = vtkCubeSource::New();

        vtkCubeSource *cuB = vtkCubeSource::New();
        cuB->SetXLength(.5);
        cuB->SetYLength(.5);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, cuA->GetOutputPort());
        bf->SetInputConnection(1, cuB->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        cuB->Delete();
        cuA->Delete();

        return ok;

    } else if (t == 14) {
        vtkCubeSource *cuA = vtkCubeSource::New();

        vtkCubeSource *cuB = vtkCubeSource::New();
        cuB->SetXLength(.25);
        cuB->SetYLength(.25);
        cuB->SetCenter(.25, 0, 0);

        vtkCubeSource *cuC = vtkCubeSource::New();
        cuC->SetXLength(.25);
        cuC->SetYLength(.25);
        cuC->SetCenter(-.25, 0, 0);

        vtkAppendPolyData *app = vtkAppendPolyData::New();
        app->AddInputConnection(cuB->GetOutputPort());
        app->AddInputConnection(cuC->GetOutputPort());

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, cuA->GetOutputPort());
        bf->SetInputConnection(1, app->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        app->Delete();
        cuC->Delete();
        cuB->Delete();
        cuA->Delete();

        return ok;

    } else if (t == 15) {

        vtkCubeSource *cu = vtkCubeSource::New();

        vtkCylinderSource *cylA = vtkCylinderSource::New();
        cylA->SetResolution(10);
        cylA->SetRadius(.4);

        vtkCylinderSource *cylB = vtkCylinderSource::New();
        cylB->SetResolution(10);
        cylB->SetRadius(.3);

        vtkPolyDataBooleanFilter *bfA = vtkPolyDataBooleanFilter::New();
        bfA->SetInputConnection(0, cylA->GetOutputPort());
        bfA->SetInputConnection(1, cylB->GetOutputPort());
        bfA->SetOperModeToDifference();

        vtkPolyDataBooleanFilter *bfB = vtkPolyDataBooleanFilter::New();
        bfB->SetInputConnection(0, cu->GetOutputPort());
        bfB->SetInputConnection(1, bfA->GetOutputPort());
        bfB->MergeRegsOn();
        bfB->Update();

        Test test(bfB->GetOutput(0), bfB->GetOutput(1));
        int ok = test.run();

        bfB->Delete();

        bfA->Delete();
        cylB->Delete();
        cylA->Delete();
        cu->Delete();

        return ok;

    } else if (t == 16) {
        double pA[] = {0.001, 10.122, 100.000};
        double pB[] = {-0.000, -10.128, 100.000};

        double pC[] = {10.125, -0.003, 99.995};
        double pD[] = {-10.124, -0.002, 100.005};

        vtkLineSource *lineA = vtkLineSource::New();
        lineA->SetPoint1(pA);
        lineA->SetPoint2(pB);

        vtkLineSource *lineB = vtkLineSource::New();
        lineB->SetPoint1(pC);
        lineB->SetPoint2(pD);

        vtkTubeFilter *tubeA = vtkTubeFilter::New();
        tubeA->SetRadius(3);
        tubeA->SetNumberOfSides(30);
        tubeA->SetInputConnection(lineA->GetOutputPort());

        vtkTubeFilter *tubeB = vtkTubeFilter::New();
        tubeB->SetRadius(3);
        tubeB->SetNumberOfSides(30);
        tubeB->SetInputConnection(lineB->GetOutputPort());

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, tubeA->GetOutputPort());
        bf->SetInputConnection(1, tubeB->GetOutputPort());
        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        tubeB->Delete();
        tubeA->Delete();
        lineB->Delete();
        lineA->Delete();

        return ok;

    } else if (t == 17) {
        // einfachste art sich überschneidender strips

        vtkCubeSource *cuA = vtkCubeSource::New();

        vtkCubeSource *cuB = vtkCubeSource::New();
        cuB->SetCenter(.3, .3, 0);

        vtkCubeSource *cuC = vtkCubeSource::New();
        cuC->SetCenter(.2, .9, 0);
        cuC->SetXLength(2);
        cuC->SetYLength(2);

        vtkAppendPolyData *app = vtkAppendPolyData::New();
        app->AddInputConnection(cuA->GetOutputPort());
        app->AddInputConnection(cuB->GetOutputPort());

        Observer *obs = Observer::New();

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, app->GetOutputPort());
        bf->SetInputConnection(1, cuC->GetOutputPort());
        bf->AddObserver(vtkCommand::ErrorEvent, obs);
        bf->Update();

        int ok = static_cast<int>(!obs->hasError);

        bf->Delete();
        obs->Delete();
        app->Delete();
        cuC->Delete();
        cuB->Delete();
        cuA->Delete();

        return ok;

    } else if (t == 18) {
        // hier überschneiden sich strips, die holes bilden
        // 4 tests in einem

        vtkCubeSource *cuA = vtkCubeSource::New();
        vtkCubeSource *cuB = vtkCubeSource::New();
        vtkCubeSource *cuC = vtkCubeSource::New();

        vtkAppendPolyData *app = vtkAppendPolyData::New();
        app->AddInputConnection(cuB->GetOutputPort());
        app->AddInputConnection(cuC->GetOutputPort());

        Observer *obs = Observer::New();

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, cuA->GetOutputPort());
        bf->SetInputConnection(1, app->GetOutputPort());
        bf->AddObserver(vtkCommand::ErrorEvent, obs);

        std::vector<std::array<double, 4>> data {
            {.8, .8, .2, 1},
            {.8, .5, .5, .8},
            {.5, .5, .5, 1},
            {.5, .5, .5, .8}
        };

        std::vector<bool> errors;

        for (auto &d : data) {
            cuB->SetXLength(d[0]);
            cuB->SetYLength(d[1]);

            cuC->SetXLength(d[2]);
            cuC->SetYLength(d[3]);

            bf->Update();

            errors.push_back(obs->hasError);

            obs->Clear();
        }

        int ok = static_cast<int>(std::find(errors.begin(), errors.end(), false) != errors.end());

        bf->Delete();
        obs->Delete();
        app->Delete();
        cuC->Delete();
        cuB->Delete();
        cuA->Delete();

        return ok;

    }  else if (t == 19) {
        vtkSphereSource *spA = vtkSphereSource::New();

        vtkSphereSource *spB = vtkSphereSource::New();
        spB->SetCenter(0, -1, 0);

        vtkAppendPolyData *app = vtkAppendPolyData::New();
        app->AddInputConnection(spA->GetOutputPort());
        app->AddInputConnection(spB->GetOutputPort());

        vtkSphereSource *spC = vtkSphereSource::New();
        spC->SetCenter(0, .5, 0);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, app->GetOutputPort());
        bf->SetInputConnection(1, spC->GetOutputPort());
        bf->SetOperModeToDifference();

        vtkPolyDataConnectivityFilter *cf = vtkPolyDataConnectivityFilter::New();
        cf->SetExtractionModeToAllRegions();
        cf->SetInputConnection(bf->GetOutputPort(0));

        cf->Update();

        int ok = cf->GetNumberOfExtractedRegions() != 3;

        cf->Delete();
        bf->Delete();
        spC->Delete();
        app->Delete();
        spB->Delete();
        spA->Delete();

        return ok;

    } else if (t == 20) {
        vtkPolyData *pd = vtkPolyData::New();
        pd->Allocate(1, 1);

        vtkPoints *pts = vtkPoints::New();
        pts->SetDataTypeToDouble();

        vtkIdList *cell = vtkIdList::New();

        std::vector<std::array<double, 2>> poly {{0, 0}, {1, -1}, {2, 0}, {3, -1}, {4, 0}, {5, 0}, {6, -1}, {7, 0}, {8, 1}, {9, 0}, {10, 0}, {11, 1}, {12, 0}, {13, 1}, {14, 0}, {14, 2}, {0, 2}};

        for (auto &pt : poly) {
            cell->InsertNextId(pts->InsertNextPoint(pt[0], pt[1], 0));
        }

        pd->SetPoints(pts);

        pd->InsertNextCell(VTK_POLYGON, cell);

        vtkTrivialProducer *prod = vtkTrivialProducer::New();
        prod->SetOutput(pd);

        vtkLinearExtrusionFilter *extr = vtkLinearExtrusionFilter::New();
        extr->SetInputConnection(prod->GetOutputPort());
        extr->SetVector(0, 0, 1);

        vtkPolyDataNormals *norm = vtkPolyDataNormals::New();
        norm->SetInputConnection(extr->GetOutputPort());
        norm->AutoOrientNormalsOn();

        vtkCubeSource *cube = vtkCubeSource::New();
        cube->SetBounds(0, 14, -2, 0, 0, 1);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, norm->GetOutputPort());
        bf->SetInputConnection(1, cube->GetOutputPort());

        bf->MergeRegsOn();
        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        bf->Delete();
        cube->Delete();
        norm->Delete();
        extr->Delete();
        prod->Delete();
        cell->Delete();
        pts->Delete();
        pd->Delete();

        return ok;

    }

}
