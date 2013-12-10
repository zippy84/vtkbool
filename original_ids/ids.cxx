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

#include <vtkSphereSource.h>
#include <vtkSuperquadricSource.h>
#include <vtkPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>

#include "vtkPolyDataContactFilter.h"
#include "vtkPolyDataBooleanFilter.h"

int main () {
    
    vtkSphereSource *sp = vtkSphereSource::New();
    sp->SetThetaResolution(20);
    sp->SetPhiResolution(20);
    sp->SetCenter(.5, 0, 0);
    sp->Update();
    
    vtkPolyData *spPd = sp->GetOutput();
    
    vtkFloatArray *spData = vtkFloatArray::New();
    spData->SetName("Data");
    
    for (unsigned int i = 0; i < spPd->GetNumberOfCells(); i++) {
        spData->InsertNextValue(1.78);
    }
    
    spPd->GetCellData()->AddArray(spData);

    vtkSuperquadricSource *sq = vtkSuperquadricSource::New();
    sq->ToroidalOn();
    sq->Update();
    
    vtkPolyData *sqPd = sq->GetOutput();
    
    vtkFloatArray *sqData = vtkFloatArray::New();
    sqData->SetName("Data");
    
    for (unsigned int i = 0; i < sqPd->GetNumberOfCells(); i++) {
        sqData->InsertNextValue(4.19);
    }
    
    sqPd->GetCellData()->AddArray(sqData);
    
    vtkPolyDataContactFilter *con = vtkPolyDataContactFilter::New();
#if (VTK_MAJOR_VERSION == 5)
    con->SetInput(0, spPd);
    con->SetInput(1, sqPd);
#else
    con->SetInputData(0, spPd);
    con->SetInputData(1, sqPd);
#endif
    
    vtkPolyDataWriter *w = vtkPolyDataWriter::New();
    w->SetInputConnection(con->GetOutputPort());
    w->SetFileName("lines.vtk");
    
    vtkPolyDataWriter *w1 = vtkPolyDataWriter::New();
    w1->SetInputConnection(con->GetOutputPort(1));
    w1->SetFileName("modPdA.vtk");
    
    vtkPolyDataWriter *w2 = vtkPolyDataWriter::New();
    w2->SetInputConnection(con->GetOutputPort(2));
    w2->SetFileName("modPdB.vtk");
    
    w->Update();
    w1->Update();
    w2->Update();
    
    vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
    bf->SetOperModeToDifference();
    bf->SetInputConnection(0, sp->GetOutputPort());
    bf->SetInputConnection(1, sq->GetOutputPort());
    
    vtkPolyDataWriter *w3 = vtkPolyDataWriter::New();
    w3->SetInputConnection(bf->GetOutputPort());
    w3->SetFileName("bool.vtk");
    w3->Update();
    
    w3->Delete();
    bf->Delete();
    w2->Delete();
    w1->Delete();
    w->Delete();
    con->Delete();
    sqData->Delete();
    sq->Delete();
    spData->Delete();
    sp->Delete();
    
    return 1;
}
