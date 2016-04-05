/*
   Copyright 2012-2016 Ronald Römer

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

#include "vtkPolyDataBooleanFilter.h"
#include "GeomHelper.h"

#include <vtkKdTreePointLocator.h>
#include <vtkCubeSource.h>
#include <vtkCylinderSource.h>
#include <vtkSphereSource.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkMath.h>
#include <vtkTrivialProducer.h>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>

#define _USE_MATH_DEFINES
#include <cmath>

typedef std::vector<int> IdsType;
typedef std::map<int, IdsType> LinksType;

class P {
public:
    P (int ind, double t) : ind(ind), t(t) {}
    int ind;
    double t;

    bool operator< (const P &other) const {
        return t < other.t;
    }
};

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
        int errA = checkConnectivity(0);
        int errB = checkConnectivity(1);

        return errA > 0 || errB > 0;
    }

private:
    int checkConnectivity (int input) {
        std::cout << "Checking input " << input << std::endl;

        int err = 0;

        vtkIntArray *origIds = vtkIntArray::SafeDownCast(input == 0 ? pd->GetCellData()->GetArray("OrigCellIdsA") : pd->GetCellData()->GetArray("OrigCellIdsB"));

        int numLines = lines->GetNumberOfCells();

        double ptA[3], ptB[3], endPt[3];

        vtkIdList *line = vtkIdList::New();

        for (int i = 0; i < numLines; i++) {

            lines->GetCellPoints(i, line);

            lines->GetPoint(line->GetId(0), ptA);
            lines->GetPoint(line->GetId(1), ptB);

            std::map<int, int> srcs;

            GeomHelper::FindPoints(loc, ptA, pts);

            LinksType links;

            for (int j = 0; j < pts->GetNumberOfIds(); j++) {
                pd->GetPointCells(pts->GetId(j), cells);

                for (int k = 0; k < cells->GetNumberOfIds(); k++) {
                    if (origIds->GetValue(cells->GetId(k)) > -1) {
                        links[cells->GetId(k)].push_back(pts->GetId(j));
                    }
                }

                srcs[pts->GetId(j)] = 0;
            }

            GeomHelper::FindPoints(loc, ptB, pts);

            for (int j = 0; j < pts->GetNumberOfIds(); j++) {
                pd->GetPointCells(pts->GetId(j), cells);

                for (int k = 0; k < cells->GetNumberOfIds(); k++) {
                    if (origIds->GetValue(cells->GetId(k)) > -1) {
                        links[cells->GetId(k)].push_back(pts->GetId(j));
                    }
                }

                srcs[pts->GetId(j)] = 1;
            }

            // links muss auf die reduziert werden, denen zwei oder mehr punkte zugeordnet sind

            LinksType::iterator itr = links.begin();
            IdsType::const_iterator itr2;

            LinksType links2;
            IdsType ids;

            while (itr != links.end()) {
                itr2 = itr->second.begin();

                while (itr2 != itr->second.end()) {
                    links2[*itr2++].push_back(itr->first);
                }

                if (itr->second.size() == 1) {
                    links.erase(itr++);
                } else {
                    ids.push_back(itr->first);
                    ++itr;
                }
            }

            if (links.size() == 0) {
                // ...
                std::cout << "E0: line " << i << " has no neighbors" << std::endl;
            } else if (links.size() == 1) {
                // kann vorkommen, aber nur in einem bestimmten fall
                pd->GetCellPoints(links.begin()->first, poly);
                IdsType &pts = links.begin()->second;

                double a[3], b[3], v[3], t = 0;

                std::vector<P> ps;

                int numPts = poly->GetNumberOfIds();

                for (int j = 0; j < numPts; j++) {
                    int k = (j+1)%numPts;
                    pd->GetPoint(poly->GetId(j), a);
                    pd->GetPoint(poly->GetId(k), b);

                    vtkMath::Subtract(a, b, v);
                    t += vtkMath::Norm(v);

                    if (std::count(pts.begin(), pts.end(), poly->GetId(k)) == 1) {
                        ps.push_back(P(k, t));
                    }

                }

                if (ps.back().ind == 0) {
                    ps.back().t = 0;
                }

                std::sort(ps.begin(), ps.end());

                int s = 0;

                std::vector<P>::const_iterator itr_;

                for (itr_ = ps.begin(); itr_ != ps.end()-1; ++itr_) {
                    if ((itr_+1)->t-itr_->t < 1e-7) {
                        ++s;
                    }
                }

                if (s > 0) {
                    // punkte wiederholen sich

                    if (s == pts.size()-1) {
                        // punkte befinden sich alle an der gleichen stelle
                        std::cout << "E1_1: line " << i << " has a polygon at one end with repeated points"
                                  << " (poly " << links.begin()->first
                                  << ", n " << (s+1)
                                  << ")" << std::endl;
                    } else {
                        // eine reale kante verhanden
                        // hier fehlt der zweite nachbar
                        std::cout << "E1_2: line " << i << " has only one neighbor with repeated points"
                                  << " (poly " << links.begin()->first
                                  << ", n " << (s+1)
                                  << ")" << std::endl;
                    }

                } else {
                    if (pts.size() == 2) {
                        // kein zweiter nachbar vorhanden
                        std::cout << "E2_1: line " << i << " has only one neighbor"
                                  << " (poly "  << links.begin()->first
                                  << ")" << std::endl;
                    } else {
                        // 3 punkte, diese müssen aber einmalig sein
                        std::set<int> uniqs(pts.begin(), pts.end());

                        if (uniqs.size() < pts.size()) {
                            std::cout << "E2_2: line " << i << " is adjacent by a polygon with folded edges - repeated indices"
                                      << " (poly "  << links.begin()->first
                                      << ")" << std::endl;
                        }
                    }
                }

            } else if (links.size() == 2) {
                // diese dürfen nicht verbunden sein
                IdsType &f = links[ids[0]];
                IdsType &s = links[ids[1]];

                if (f.size() > 2 || s.size() > 2) {
                    std::cout << "E3: line " << i << " has neighbors with repeated points"
                              << " (polyA " << ids[0] << " -> n " << f.size()
                              << ", polyB " << ids[1] << " -> n " << s.size()
                              << ")" << std::endl;
                    ++err;
                }

                lines->GetPointCells(line->GetId(0), cells);
                int linkedA = cells->GetNumberOfIds();

                lines->GetPointCells(line->GetId(1), cells);
                int linkedB = cells->GetNumberOfIds();

                if (linkedA == linkedB) {
                    if (f[0] == s[0] || f[1] == s[1]) {
                        std::cout << "E4_1" << std::endl;
                        ++err;
                    }
                } else if (linkedA == 4) {
                    if (f[0] == s[0]) {
                        std::cout << "E4_2" << std::endl;
                        ++err;
                    }
                } else {
                    if (f[1] == s[1]) {
                        std::cout << "E4_3" << std::endl;
                        ++err;
                    }
                }

            } else {
                // kann nicht vorkommen
            }

            // man muss prüfen ob alle an einem punkt beteiligten polygone aneinandergrenzen

            std::set<int> endsA, endsB;

            for (int j = 0; j < 2; j++) {
                std::set<int> &ends = j == 0 ? endsA : endsB;

                lines->GetPointCells(line->GetId(j), cells);

                for (int k = 0; k < cells->GetNumberOfIds(); k++) {
                    lines->GetCellPoints(cells->GetId(k), poly);
                    int end = poly->GetId(0) == line->GetId(j) ? 1 : 0;

                    lines->GetPoint(poly->GetId(end), endPt);
                    GeomHelper::FindPoints(loc, endPt, pts);

                    for (int l = 0; l < pts->GetNumberOfIds(); l++) {
                        ends.insert(pts->GetId(l));
                    }
                }
            }

            for (itr = links2.begin(); itr != links2.end(); ++itr) {

                std::map<int, int> counts;

                for (itr2 = itr->second.begin(); itr2 != itr->second.end(); ++itr2) {
                    pd->GetCellPoints(*itr2, poly);

                    int numPts = poly->GetNumberOfIds();

                    for (int j = 0; j < numPts; j++) {
                        if (poly->GetId(j) == itr->first) {
                            counts[poly->GetId((j+1)%numPts)]++;
                            counts[poly->GetId((j+numPts-1)%numPts)]++;
                            break;
                        }
                    }
                }

                std::set<int> &ends = srcs[itr->first] == 0 ? endsA : endsB;

                std::map<int, int>::const_iterator itr3;

                int filtered = 0;

                IdsType ones;

                for (itr3 = counts.begin(); itr3 != counts.end(); ++itr3) {
                    if (itr3->second == 1) {
                        if (ends.count(itr3->first) == 1) {
                            filtered++;
                        } else {
                            ones.push_back(itr3->first);
                        }
                    }
                }

                if (filtered < 2) {
                    // die region wird nicht von zwei linien begrenzt
                    std::cout << "E5_1: polygons [";

                    for (itr2 = itr->second.begin(); itr2 != itr->second.end(); ++itr2) {
                        std::cout << *itr2 << ", ";
                    }

                    std::cout << "] on point " << itr->first << " are not bounded by two lines" << std::endl;

                    ++err;

                } else if (ones.size() > 0) {
                    // ganze region existiert, es sind aber nicht alle punkte miteinander verbunden
                    std::cout << "E5_2: polygons [";

                    for (itr2 = itr->second.begin(); itr2 != itr->second.end(); ++itr2) {
                        std::cout << *itr2 << ", ";
                    }

                    std::cout << "] on point " << itr->first << " are not connected (ones [";

                    for (itr2 = ones.begin(); itr2 != ones.end(); ++itr2) {
                        std::cout << *itr2 << ", ";
                    }

                    std::cout << "])" << std::endl;

                    ++err;
                }

            }

        }

        line->Delete();

        return err;

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
        bf->MergeAllOn();

        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        //GeomHelper::WriteVTK("test1.vtk", bf->GetOutput(0));

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
        bf->MergeAllOn();

        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        //GeomHelper::WriteVTK("test2.vtk", bf->GetOutput(0));

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
        bf->MergeAllOn();

        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        //GeomHelper::WriteVTK("test3.vtk", bf->GetOutput(0));

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
        bf->MergeAllOn();

        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        //GeomHelper::WriteVTK("test4.vtk", bf->GetOutput(0));

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
        bf->MergeAllOn();

        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        //GeomHelper::WriteVTK("test5.vtk", bf->GetOutput(0));

        bf->Delete();
        cubeB->Delete();
        cubeA->Delete();

        return ok;

    } else if (t == 5) {
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
            double x0 = .5*cos(i*2*M_PI/32);
            double z0 = .5*sin(i*2*M_PI/32);

            double x1 = .5*cos((i+1)*2*M_PI/32);
            double z1 = .5*sin((i+1)*2*M_PI/32);

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
        bf->MergeAllOn();

        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        //GeomHelper::WriteVTK("test6.vtk", bf->GetOutput(0));

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
        bf->MergeAllOn();

        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        //GeomHelper::WriteVTK("test7.vtk", bf->GetOutput(0));

        bf->Delete();
        sp->Delete();
        cu->Delete();

        return ok;

    }  else if (t == 7) {
        // enthält schmale schnitte an den ecken der polygone
        // gelöst mit einführung von CAPT_*

        // außerdem ist hier mind. eine überlappung enthalten
        // siehe ResolveOverlaps

        vtkSphereSource *spA = vtkSphereSource::New();
        spA->SetRadius(50);
        spA->SetPhiResolution(100);
        spA->SetThetaResolution(100);

        vtkSphereSource *spB = vtkSphereSource::New();
        spB->SetRadius(50);
        spB->SetCenter(25, 0, 0);
        spB->SetPhiResolution(10);
        spB->SetThetaResolution(10);

        // die gleichen probleme kann man auch mit 79, 215 und 296 statt 10 erzeugen

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, spA->GetOutputPort());
        bf->SetInputConnection(1, spB->GetOutputPort());
        bf->MergeAllOn();

        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        //GeomHelper::WriteVTK("test8.vtk", bf->GetOutput(0));

        bf->Delete();
        spB->Delete();
        spA->Delete();

        return ok;

    } else if (t == 8) {
        // bei linie 50 wird auch ein polygon übersprungen

        vtkSphereSource *spA = vtkSphereSource::New();
        spA->SetRadius(50);
        spA->SetPhiResolution(100);
        spA->SetThetaResolution(100);

        vtkSphereSource *spB = vtkSphereSource::New();
        spB->SetRadius(50);
        spB->SetCenter(25, 0, 0);
        spB->SetPhiResolution(173);
        spB->SetThetaResolution(173);

        vtkPolyDataBooleanFilter *bf = vtkPolyDataBooleanFilter::New();
        bf->SetInputConnection(0, spA->GetOutputPort());
        bf->SetInputConnection(1, spB->GetOutputPort());
        bf->MergeAllOn();

        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        //GeomHelper::WriteVTK("test9.vtk", bf->GetOutput(0));

        bf->Delete();
        spB->Delete();
        spA->Delete();

        return ok;

    } else if (t == 9) {
        // enthält sehr scharfwinklige schnitte mit kanten
        // ein punkt im toleranzbereich einer ecke, der andere in der mitte einer angrenzenden kante

        // gelöst mit ts

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
        bf->MergeAllOn();

        bf->Update();

        Test test(bf->GetOutput(0), bf->GetOutput(1));
        int ok = test.run();

        //GeomHelper::WriteVTK("test10.vtk", bf->GetOutput(0));

        bf->Delete();
        spB->Delete();
        spA->Delete();

        return ok;

    }
}
