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

#include <vtkStaticCellLocator.h>
#include <vtkCellIterator.h>

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

class Test {
    vtkPolyData *pdA, *pdB, *lines;
    vtkIdTypeArray *contsA, *contsB;

    vtkIdTypeArray *origCellIdsA, *origCellIdsB;

public:
    Test (vtkPolyDataBooleanFilter *bf) {
        bf->Update();

        lines = bf->GetOutput(0);
        pdA = bf->GetOutput(1);
        pdB = bf->GetOutput(2);

        contsA = vtkIdTypeArray::SafeDownCast(lines->GetCellData()->GetArray("cA"));
        contsB = vtkIdTypeArray::SafeDownCast(lines->GetCellData()->GetArray("cB"));

        lines->BuildLinks();

        pdA->BuildLinks();
        pdB->BuildLinks();

        origCellIdsA = vtkIdTypeArray::SafeDownCast(pdA->GetCellData()->GetArray("OrigCellIds"));
        origCellIdsB = vtkIdTypeArray::SafeDownCast(pdB->GetCellData()->GetArray("OrigCellIds"));

        WriteVTK("pdA.vtk", pdA);
        WriteVTK("pdB.vtk", pdB);
        WriteVTK("lines.vtk", lines);
    }

    int run () {
        bool errA = checkConnectivity(pdA, origCellIdsA);
        bool errB = checkConnectivity(pdB, origCellIdsB);

        return errA || errB;
    }

private:
    bool checkConnectivity (vtkPolyData *pd, vtkIdTypeArray *vtkNotUsed(origCellIds)) {
        std::cout << "Checking connectivity" << std::endl;

        bool errored = false;

        vtkStaticCellLocator *loc = vtkStaticCellLocator::New();
        loc->SetDataSet(pd);
        loc->BuildLocator();

        double pA[3], pB[3];

        vtkIdList *neigs = vtkIdList::New();

        vtkCellIterator *itr = lines->NewCellIterator();
        for (itr->InitTraversal(); !itr->IsDoneWithTraversal(); itr->GoToNextCell()) {
            vtkIdList *pts = itr->GetPointIds();

            pd->GetPoint(pts->GetId(0), pA);
            pd->GetPoint(pts->GetId(1), pB);

            loc->FindCellsAlongLine(pA, pB, 1e-6, neigs);

            std::cout << itr->GetCellId() << " -> [";

            vtkIdType i, numCells = neigs->GetNumberOfIds();

            for (i = 0; i < numCells; i++) {
                std::cout << neigs->GetId(i) << ", ";
            }

            std::cout << "]" << std::endl;
        }

        loc->FreeSearchStructure();
        loc->Delete();

        neigs->Delete();
        itr->Delete();

        return errored;

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

    void Clear () {
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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bf->SetOperModeToNone();

        Test test(bf);
        int ok = test.run();

        bf->Delete();
        cyl->Delete();
        cu->Delete();

        return ok;

    } else if (t == 11) {
        // hier ist ein bug drin

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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bfB->SetOperModeToNone();

        Test test(bfB);
        int ok = test.run();

        bfB->Delete();

        bfA->Delete();
        cylB->Delete();
        cylA->Delete();
        cu->Delete();

        return ok;

    } else if (t == 16) {
        // geht auch nicht

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
        bf->SetOperModeToNone();

        Test test(bf);
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
        bf->SetOperModeToNone();

        Test test(bf);
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
