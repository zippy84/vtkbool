/*
Copyright 2012-2025 Ronald Römer

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

#include "Contact.h"

#include <vtkCleanPolyData.h>

int main() {
    auto pdA = CreatePolyData({ { {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0} } });

    double z = std::tan(0.0005/180*M_PI)*.5;

    auto pdB = CreatePolyData({
        { {0, 0, 0}, {0, .5, 0}, {1, .5, 0}, {1, 0, 0} },
        { {1, .5, 0}, {0, .5, 0}, {0, 1, z}, {1, 1, z} }
    });

    auto clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInputData(pdB);
    clean->Update();

    auto lines = Contact(pdA, clean->GetOutput()).GetLines();

    if (lines->GetNumberOfCells() == 0) {
        return EXIT_FAILURE;
    }

#ifdef DEBUG
    WriteVTK("lines.vtk", lines);
#endif

    return EXIT_SUCCESS;
}
