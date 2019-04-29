/*
Copyright 2012-2019 Ronald Römer

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
#include <vtkTriangleFilter.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkVectorText.h>
#include <vtkPolyDataNormals.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkPlaneSource.h>
#include <vtkAppendPolyData.h>
#include <vtkLineSource.h>
#include <vtkTubeFilter.h>
#include <vtkCommand.h>

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

// #define DD

int Point::_tag = 0;

typedef std::map<int, IdsType> LinksType;

class Test {
    vtkPolyData *lines, *pd;
    vtkKdTreePointLocator *loc;
    vtkIdList *cells, *poly, *pts;

    vtkIntArray *contsA, *contsB;

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

        contsA = vtkIntArray::SafeDownCast(lines->GetCellData()->GetArray("cA"));
        contsB = vtkIntArray::SafeDownCast(lines->GetCellData()->GetArray("cB"));
    }
    ~Test () {
        loc->FreeSearchStructure();
        loc->Delete();
        cells->Delete();
        poly->Delete();
        pts->Delete();
    }

    int run () {

#ifdef DEBUG
        WriteVTK("merged.vtk", pd);
#endif

        vtkIntArray *origCellIdsA = vtkIntArray::SafeDownCast(pd->GetCellData()->GetArray("OrigCellIdsA"));
        vtkIntArray *origCellIdsB = vtkIntArray::SafeDownCast(pd->GetCellData()->GetArray("OrigCellIdsB"));

        int errA = checkConnectivity(0, origCellIdsA);
        int errB = checkConnectivity(1, origCellIdsB);

        return errA > 0 || errB > 0;

    }

private:
    int checkConnectivity (int input, vtkIntArray *origIds) {
        std::cout << "Checking input " << input << std::endl;

        int err = 0;

        for (int i = 0; i < pd->GetNumberOfCells() && err < 1; i++) {
            if (origIds->GetValue(i) > -1) {
                pd->GetCellPoints(i, poly);
                if (poly->GetNumberOfIds() < 3) {
                    std::cout << "poly " << i << " has too few points" << std::endl;

                    err++;
                }
            }
        }

        if (err == 0) {

            int numLines = lines->GetNumberOfCells();

            double ptA[3], ptB[3];

            vtkIdList *line = vtkIdList::New();

            for (int i = 0; i < numLines; i++) {

                lines->GetCellPoints(i, line);

                lines->GetPoint(line->GetId(0), ptA);
                lines->GetPoint(line->GetId(1), ptB);

                FindPoints(loc, ptA, pts);

                LinksType links;

                // sammelt polygone und deren punkte am jeweiligen ende der linie

                // bei einem normalen schnitt sollte jedes polygon zwei punkte haben

                for (int j = 0; j < pts->GetNumberOfIds(); j++) {
                    pd->GetPointCells(pts->GetId(j), cells);

                    for (int k = 0; k < cells->GetNumberOfIds(); k++) {
                        if (origIds->GetValue(cells->GetId(k)) > -1) {
                            links[cells->GetId(k)].push_back(pts->GetId(j));
                        }
                    }
                }

                FindPoints(loc, ptB, pts);

                for (int j = 0; j < pts->GetNumberOfIds(); j++) {
                    pd->GetPointCells(pts->GetId(j), cells);

                    for (int k = 0; k < cells->GetNumberOfIds(); k++) {
                        if (origIds->GetValue(cells->GetId(k)) > -1) {
                            links[cells->GetId(k)].push_back(pts->GetId(j));
                        }
                    }
                }

                LinksType::const_iterator itr;

                IdsType polys;

                for (itr = links.begin(); itr != links.end(); ++itr) {
                    const IdsType &pts = itr->second;
                    if (pts.size() > 1) {
                        // gibt es doppelte punkte

                        pd->GetCellPoints(itr->first, poly);
                        int numPts = poly->GetNumberOfIds();

                        double a[3], b[3];
                        for (int j = 0; j < numPts; j++) {
                            pd->GetPoint(poly->GetId(j), a);
                            pd->GetPoint(poly->GetId((j+1)%numPts), b);

                            if (GetD(a, b) < 1e-6) {
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
                    IdsType &f = links[polys[0]],
                            &s = links[polys[1]];

                    lines->GetPointCells(line->GetId(0), cells);
                    int linkedA = cells->GetNumberOfIds();

                    lines->GetPointCells(line->GetId(1), cells);
                    int linkedB = cells->GetNumberOfIds();

                    if (linkedA == linkedB) {
                        std::set<int> all;
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
                    IdsType &f = links[polys[0]];

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

                IdsType::const_iterator itr2, itr3;

                for (int j = 0; j < 2 && err < 1; j++) {
                    double pt[3];
                    lines->GetPoint(line->GetId(j), pt);

                    IdsType ids;

                    FindPoints(loc, pt, pts);
                    int numPts = pts->GetNumberOfIds();

                    std::map<int, int> masked;

                    for (int k = 0; k < numPts; k++) {

                        pd->GetPointCells(pts->GetId(k), cells);
                        int numCells = cells->GetNumberOfIds();

                        for (int l = 0; l < numCells; l++) {
                            int cellId = cells->GetId(l);

                            if (origIds->GetValue(cellId) > -1) {

                                pd->GetCellPoints(cells->GetId(l), poly);
                                int numPts_ = poly->GetNumberOfIds();

                                for (int m = 0; m < numPts_; m++) {
                                    int ptId = poly->GetId(m);

                                    /* wenn ein polygon mehrfach die gleiche ptId verwendet,
                                       dann erscheint die cellId auch mehrfach in cells
                                       => es sind dann auch 4 kanten zu verarbeiten
                                       => der betrachtete relative index wird dann zusammen mit der cellId in masked gespeichert */

                                    if (masked.count(cellId) == 1 && masked[cellId] == m) {
                                        continue;
                                    }

                                    if (ptId == pts->GetId(k)) {

                                        ids.push_back(poly->GetId((m+1)%numPts_));
                                        ids.push_back(poly->GetId((m+numPts_-1)%numPts_));

                                        masked[cellId] = m;

                                        break;
                                    }
                                }

                            }
                        }

                    }

                    std::sort(ids.begin(), ids.end());

                    IdsType ids2 = ids;
                    ids2.erase(std::unique(ids2.begin(), ids2.end()), ids2.end());

                    IdsType ids3;
                    for (itr2 = ids2.begin(); itr2 != ids2.end(); ++itr2) {
                        if (std::count(ids.begin(), ids.end(), *itr2) == 1) {
                            ids3.push_back(*itr2);
                        }
                    }

                    int c = 0;

                    for (itr2 = ids3.begin(); itr2 != ids3.end(); ++itr2) {
                        for (itr3 = itr2+1; itr3 != ids3.end(); ++itr3) {
                            double a[3], b[3];
                            pd->GetPoint(*itr2, a);
                            pd->GetPoint(*itr3, b);

                            if (GetD(a, b) < 1e-6) {
                                c++;
                            }
                        }
                    }

                    if (ids3.size() != 2*c) {
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

    Observer() : hasError(false) {}

    static Observer *New() {
        return new Observer;
    }

    virtual void Execute(vtkObject *vtkNotUsed(caller), unsigned long event, void *calldata) {
        hasError = event == vtkCommand::ErrorEvent;
        msg = static_cast<char*>(calldata);
    }

    void Clear() {
        hasError = false;
        msg.clear();
    }
};

int main (int argc, char *argv[]) {
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

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test0.vtk", bf->GetOutput(0));
#endif

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

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test1.vtk", bf->GetOutput(0));
#endif

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

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test2.vtk", bf->GetOutput(0));
#endif

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

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test3.vtk", bf->GetOutput(0));
#endif

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

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test4.vtk", bf->GetOutput(0));
#endif

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
            double x0 = .5*std::cos(i*2*PI/32);
            double z0 = .5*std::sin(i*2*PI/32);

            double x1 = .5*std::cos((i+1)*2*PI/32);
            double z1 = .5*std::sin((i+1)*2*PI/32);

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

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test5.vtk", bf->GetOutput(0));
#endif

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

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test6.vtk", bf->GetOutput(0));
#endif

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

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test7.vtk", bf->GetOutput(0));
#endif

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

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test8.vtk", bf->GetOutput(0));
#endif

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

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test9.vtk", bf->GetOutput(0));
#endif

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

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test10.vtk", bf->GetOutput(0));
#endif

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

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test11.vtk", bf->GetOutput(0));
#endif

        bf->Delete();
        sf->Delete();
        tf->Delete();
        cubeB->Delete();
        cubeA->Delete();

        return ok;

    } else if (t == 12) {
        vtkVectorText *text = vtkVectorText::New();
        text->SetText("zippy84");

        vtkLinearExtrusionFilter *ef = vtkLinearExtrusionFilter::New();
        ef->SetInputConnection(text->GetOutputPort());
        ef->SetExtrusionTypeToNormalExtrusion();
        ef->SetVector(0, 0, 1);

        // TODO: vtkPolyDataBooleanFilter muss prüfen, ob polygone konsistent orientiert sind
        // -> AddAdjacentPoints und CombineRegions funktionieren sonst nicht
        vtkPolyDataNormals *nf = vtkPolyDataNormals::New();
        nf->SetInputConnection(ef->GetOutputPort());
        nf->FlipNormalsOn();
        nf->ConsistencyOn();

        nf->Update();

        double ext[6];
        nf->GetOutput()->GetPoints()->GetBounds(ext);

        double dx = ext[0]+(ext[1]-ext[0])/2,
            dy = ext[2]+(ext[3]-ext[2])/2;

        vtkCubeSource *cube = vtkCubeSource::New();
        cube->SetXLength(10);
        cube->SetYLength(5);
        cube->SetCenter(dx, dy, .5);

        vtkTriangleFilter *tf = vtkTriangleFilter::New();
        tf->SetInputConnection(cube->GetOutputPort());

        vtkLinearSubdivisionFilter *sf = vtkLinearSubdivisionFilter::New();
        sf->SetInputConnection(tf->GetOutputPort());
        sf->SetNumberOfSubdivisions(3);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, nf->GetOutputPort());

        // ginge mit dem merger
        bf->SetInputConnection(1, cube->GetOutputPort());

        // bf->SetInputConnection(1, sf->GetOutputPort());

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference2();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test12.vtk", bf->GetOutput(0));
#endif

        bf->Delete();
        sf->Delete();
        tf->Delete();
        cube->Delete();
        ef->Delete();
        nf->Delete();
        text->Delete();

        return ok;

    } else if (t == 13) {
        vtkSphereSource *sp = vtkSphereSource::New();
        sp->SetRadius(.5);
        sp->SetPhiResolution(8);
        sp->SetThetaResolution(8);

        vtkPlaneSource *pl = vtkPlaneSource::New();
        pl->SetResolution(5, 5);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, sp->GetOutputPort());
        bf->SetInputConnection(1, pl->GetOutputPort());

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToUnion();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test13.vtk", bf->GetOutput(0));
#endif

        bf->Delete();
        pl->Delete();
        sp->Delete();

        return ok;

    } else if (t == 14) {
        // testet den merger

        vtkCubeSource *cuA = vtkCubeSource::New();

        vtkCubeSource *cuB = vtkCubeSource::New();
        cuB->SetXLength(.5);
        cuB->SetYLength(.5);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, cuA->GetOutputPort());
        bf->SetInputConnection(1, cuB->GetOutputPort());

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test14.vtk", bf->GetOutput(0));
#endif

        bf->Delete();
        cuB->Delete();
        cuA->Delete();

        return ok;

    } else if (t == 15) {
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

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test15.vtk", bf->GetOutput(0));
#endif

        bf->Delete();
        app->Delete();
        cuC->Delete();
        cuB->Delete();
        cuA->Delete();

        return ok;

    } else if (t == 16) {

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

#ifndef DD
        bfB->MergeRegsOn();
#else
        bfB->SetOperModeToDifference();
#endif

        bfB->Update();

#ifndef DD
        Test test(bfB->GetOutput(0), bfB->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test16.vtk", bfB->GetOutput(0));
#endif

        bfB->Delete();

        bfA->Delete();
        cylB->Delete();
        cylA->Delete();
        cu->Delete();

        return ok;

    } else if (t == 17) {
        // testet _Test()

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

#ifndef DD
        bf->MergeRegsOn();
#else
        bf->SetOperModeToDifference();
#endif

        bf->Update();

#ifndef DD
        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();
#else
        int ok = 0;
        WriteVTK("test17.vtk", bf->GetOutput(0));
#endif

        bf->Delete();
        tubeB->Delete();
        tubeA->Delete();
        lineB->Delete();
        lineA->Delete();

        return ok;

    } else if (t == 18) {
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

    } else if (t == 19) {
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

        std::vector<std::array<double, 4>> data{
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
    }

}
