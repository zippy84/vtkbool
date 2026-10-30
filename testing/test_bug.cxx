#include "Contact2.h"

#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkStripper.h>

int main () {
    {
        auto readerA = vtkSmartPointer<vtkPolyDataReader>::New();
        readerA->SetFileName("../testing/data_bug/MeshA.vtk");

        auto readerB = vtkSmartPointer<vtkPolyDataReader>::New();
        readerB->SetFileName("../testing/data_bug/MeshB.vtk");

        readerA->Update();
        readerB->Update();

        Contact(readerA->GetOutput(), readerB->GetOutput()).SaveAll("Test1_");
    }

    {
        auto sphereA = vtkSmartPointer<vtkSphereSource>::New();
        auto sphereB = vtkSmartPointer<vtkSphereSource>::New();

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

        stripper->Update();

        Contact(sphereA->GetOutput(), stripper->GetOutput()).SaveAll("Test2_");
    }

    {
        auto transform = vtkSmartPointer<vtkTransform>::New();
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

        Contact(tfA->GetPolyDataOutput(), tfB->GetPolyDataOutput()).SaveAll("Test3_");
    }

}
