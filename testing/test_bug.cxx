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

#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkStripper.h>

#include "vtkPolyDataBooleanFilter.h"

int main () {
    {
        auto readerA = vtkSmartPointer<vtkPolyDataReader>::New();
        readerA->SetFileName("../testing/data_bug/MeshA.vtk");

        auto readerB = vtkSmartPointer<vtkPolyDataReader>::New();
        readerB->SetFileName("../testing/data_bug/MeshB.vtk");

        readerA->Update();
        readerB->Update();

        // Contact(readerA->GetOutput(), readerB->GetOutput()).SaveAll("Test1_");

        auto bf = vtkSmartPointer<vtkPolyDataBooleanFilter>::New();
        bf->SetInputConnection(0, readerA->GetOutputPort());
        bf->SetInputConnection(1, readerB->GetOutputPort());
        bf->Update();

        if (bf->GetOutput()->GetNumberOfCells() == 0) {
            return EXIT_FAILURE;
        }

        auto writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetInputConnection(bf->GetOutputPort());
        writer->SetFileName("Test1_result.vtk");
        writer->Update();
    }

    {
        /*auto sphereA = vtkSmartPointer<vtkSphereSource>::New();
        auto sphereB = vtkSmartPointer<vtkSphereSource>::New();

        sphereA->SetOutputPointsPrecision(vtkAlgorithm::DOUBLE_PRECISION);
        sphereB->SetOutputPointsPrecision(vtkAlgorithm::DOUBLE_PRECISION);

        const vtkIdType resolution = 30;

        sphereA->SetThetaResolution(resolution);
        sphereA->SetPhiResolution(resolution);
        sphereA->LatLongTessellationOn();

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

        auto stripper = vtkSmartPointer<vtkStripper>::New();
        stripper->SetInputConnection(tf->GetOutputPort());
        stripper->SetMaximumLength(10);

        stripper->Update();*/

        // Contact(sphereA->GetOutput(), stripper->GetOutput()).SaveAll("Test2_");

        /*auto bf = vtkSmartPointer<vtkPolyDataBooleanFilter>::New();
        bf->SetInputConnection(0, sphereA->GetOutputPort());
        bf->SetInputConnection(1, stripper->GetOutputPort());
        bf->Update();

        if (bf->GetOutput()->GetNumberOfCells() == 0) {
            return EXIT_FAILURE;
        }

        auto writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetInputConnection(bf->GetOutputPort());
        writer->SetFileName("Test2_result.vtk");
        writer->Update();*/
    }

    {
        /*auto transform = vtkSmartPointer<vtkTransform>::New();
        transform->RotateX(5);
        transform->RotateY(5);

        auto readerA = vtkSmartPointer<vtkPolyDataReader>::New();
        readerA->SetFileName("../testing/data_bug/CubeA.vtk");

        auto readerB = vtkSmartPointer<vtkPolyDataReader>::New();
        readerB->SetFileName("../testing/data_bug/CubeB.vtk");

        auto tfA = vtkSmartPointer<vtkTransformFilter>::New();
        tfA->SetInputConnection(readerA->GetOutputPort());
        tfA->SetTransform(transform);

        auto tfB = vtkSmartPointer<vtkTransformFilter>::New();
        tfB->SetInputConnection(readerB->GetOutputPort());
        tfB->SetTransform(transform);

        tfA->Update();
        tfB->Update();

        Contact contact(tfA->GetPolyDataOutput(), tfB->GetPolyDataOutput());

        if (contact.GetLines() == nullptr) {
            return EXIT_FAILURE;
        }

        contact.SaveAll("Test3_");*/
    }

    return EXIT_SUCCESS;

}
