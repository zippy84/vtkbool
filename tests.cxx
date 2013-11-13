/*
   Copyright 2012, 2013 Ronald Römer

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

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>

#include <vtkSphereSource.h>
#include <vtkPolyDataWriter.h>
#include <vtkCubeSource.h>
#include <vtkCylinderSource.h>
#include <vtkSuperquadricSource.h>

#include "vtkPolyDataBooleanFilter.h"

int main(int argc, char *argv[]) {

    std::stringstream ss0, ss1;
    int example, mode;

    ss0 << "example_";
    ss1 << "example_";

    if (argc != 3) {
        std::cerr << "The two arguments example and mode are required." << std::endl;
        return 1;
    }

    example = argv[1][0]-48;

    if (std::strlen(argv[1]) > 1 || example < 0 || example > 6) {
        std:cerr << "Example must be a single-digit number between 0 and 6." << std::endl;
        return 1;
    }

    ss0 << example << "_";

    if (std::strcmp(argv[2], "union") == 0) {
        mode = 0;
        ss0 << argv[2];
    } else if (std::strcmp(argv[2], "intersection") == 0) {
        mode = 1;
        ss0 << argv[2];
    } else if (std::strcmp(argv[2], "difference") == 0) {
        mode = 2;
        ss0 << argv[2];
    } else if (std::strcmp(argv[2], "difference2") == 0) {
        mode = 3;
        ss0 << argv[2];
    } else {
        std::cerr << "Unknown mode. union, intersection, difference and difference2 are allowed." << std::endl;
        return 1;
    }

    ss0 << "_bool.vtk";
    ss1 << example << "_cont.vtk";

    std::string fn0(ss0.str());
    std::string fn1(ss1.str());

    if (example == 0) {

        vtkCubeSource *cu = vtkCubeSource::New();
        cu->SetYLength(.5);

        vtkCylinderSource *cyl = vtkCylinderSource::New();
        cyl->SetResolution(32);
        cyl->SetHeight(.5);
        cyl->SetCenter(0, .5, 0);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetOperMode(mode);
        bf->SetInputConnection(0, cu->GetOutputPort());
        bf->SetInputConnection(1, cyl->GetOutputPort());

        vtkPolyDataWriter *w0 = vtkPolyDataWriter::New();
        w0->SetInputConnection(bf->GetOutputPort());
        w0->SetFileName(fn0.c_str());
        w0->Update();

        vtkPolyDataWriter *w1 = vtkPolyDataWriter::New();
        w1->SetInputConnection(bf->GetOutputPort(1));
        w1->SetFileName(fn1.c_str());
        w1->Update();

        w1->Delete();
        w0->Delete();
        bf->Delete();
        cyl->Delete();
        cu->Delete();

    } else if (example == 1) {

        vtkCubeSource *cu = vtkCubeSource::New();
        cu->SetYLength(.5);

        vtkCylinderSource *cyl = vtkCylinderSource::New();
        cyl->SetResolution(32);
        cyl->SetHeight(.5);
        cyl->SetCenter(0, .25, 0);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetOperMode(mode);
        bf->SetInputConnection(0, cu->GetOutputPort());
        bf->SetInputConnection(1, cyl->GetOutputPort());

        vtkPolyDataWriter *w0 = vtkPolyDataWriter::New();
        w0->SetInputConnection(bf->GetOutputPort());
        w0->SetFileName(fn0.c_str());
        w0->Update();

        vtkPolyDataWriter *w1 = vtkPolyDataWriter::New();
        w1->SetInputConnection(bf->GetOutputPort(1));
        w1->SetFileName(fn1.c_str());
        w1->Update();

        w1->Delete();
        w0->Delete();
        bf->Delete();
        cyl->Delete();
        cu->Delete();

    } else if (example == 2) {

        vtkCubeSource *cu = vtkCubeSource::New();
        cu->SetYLength(.5);

        vtkCylinderSource *cyl = vtkCylinderSource::New();
        cyl->SetResolution(32);
        cyl->SetHeight(.5);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetOperMode(mode);
        bf->SetInputConnection(0, cu->GetOutputPort());
        bf->SetInputConnection(1, cyl->GetOutputPort());

        vtkPolyDataWriter *w0 = vtkPolyDataWriter::New();
        w0->SetInputConnection(bf->GetOutputPort());
        w0->SetFileName(fn0.c_str());
        w0->Update();

        vtkPolyDataWriter *w1 = vtkPolyDataWriter::New();
        w1->SetInputConnection(bf->GetOutputPort(1));
        w1->SetFileName(fn1.c_str());
        w1->Update();

        w1->Delete();
        w0->Delete();
        bf->Delete();
        cyl->Delete();
        cu->Delete();

    } else if (example == 3) {

        vtkCubeSource *cu = vtkCubeSource::New();

        vtkSuperquadricSource *sq = vtkSuperquadricSource::New();
        sq->SetSize(1.125);
        sq->ToroidalOn();

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetOperMode(mode);
        bf->SetInputConnection(0, cu->GetOutputPort());
        bf->SetInputConnection(1, sq->GetOutputPort());

        vtkPolyDataWriter *w0 = vtkPolyDataWriter::New();
        w0->SetInputConnection(bf->GetOutputPort());
        w0->SetFileName(fn0.c_str());
        w0->Update();

        vtkPolyDataWriter *w1 = vtkPolyDataWriter::New();
        w1->SetInputConnection(bf->GetOutputPort(1));
        w1->SetFileName(fn1.c_str());
        w1->Update();

        w1->Delete();
        w0->Delete();
        bf->Delete();
        sq->Delete();
        cu->Delete();

    } else if (example == 4) {

        vtkCylinderSource *cyl = vtkCylinderSource::New();
        cyl->SetResolution(30);

        vtkCubeSource *cu = vtkCubeSource::New();
        cu->SetXLength(.5);
        cu->SetYLength(.5);
        cu->SetZLength(2);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetOperMode(mode);
        bf->SetInputConnection(0, cyl->GetOutputPort());
        bf->SetInputConnection(1, cu->GetOutputPort());

        vtkPolyDataWriter *w0 = vtkPolyDataWriter::New();
        w0->SetInputConnection(bf->GetOutputPort());
        w0->SetFileName(fn0.c_str());
        w0->Update();

        vtkPolyDataWriter *w1 = vtkPolyDataWriter::New();
        w1->SetInputConnection(bf->GetOutputPort(1));
        w1->SetFileName(fn1.c_str());
        w1->Update();

        w1->Delete();
        w0->Delete();
        bf->Delete();
        cu->Delete();
        cyl->Delete();

    } else if (example == 5) {

        vtkCubeSource *cuA = vtkCubeSource::New();

        vtkCubeSource *cuB = vtkCubeSource::New();
        cuB->SetCenter(0, 1, 0);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetOperMode(mode);
        bf->SetInputConnection(0, cuA->GetOutputPort());
        bf->SetInputConnection(1, cuB->GetOutputPort());

        vtkPolyDataWriter *w0 = vtkPolyDataWriter::New();
        w0->SetInputConnection(bf->GetOutputPort());
        w0->SetFileName(fn0.c_str());
        w0->Update();

        vtkPolyDataWriter *w1 = vtkPolyDataWriter::New();
        w1->SetInputConnection(bf->GetOutputPort(1));
        w1->SetFileName(fn1.c_str());
        w1->Update();

        w1->Delete();
        w0->Delete();
        bf->Delete();
        cuB->Delete();
        cuA->Delete();

    } else {

        vtkSphereSource *spA = vtkSphereSource::New();
        spA->SetPhiResolution(30);
        spA->SetThetaResolution(30);
        //spA->LatLongTessellationOn();

        vtkSphereSource *spB = vtkSphereSource::New();
        spB->SetCenter(.025, 0, 0);
        spB->SetPhiResolution(8);
        spB->SetThetaResolution(8);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetOperMode(mode);
        bf->SetInputConnection(0, spA->GetOutputPort());
        bf->SetInputConnection(1, spB->GetOutputPort());

        vtkPolyDataWriter *w0 = vtkPolyDataWriter::New();
        w0->SetInputConnection(bf->GetOutputPort());
        w0->SetFileName(fn0.c_str());
        w0->Update();

        vtkPolyDataWriter *w1 = vtkPolyDataWriter::New();
        w1->SetInputConnection(bf->GetOutputPort(1));
        w1->SetFileName(fn1.c_str());
        w1->Update();

        w1->Delete();
        w0->Delete();
        bf->Delete();
        spB->Delete();
        spA->Delete();

    }

    return 0;

}
