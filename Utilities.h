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

#ifndef __Utilities_h
#define __Utilities_h

#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkKdTreePointLocator.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <tuple>

#define CPY(a, b) a[0] = b[0]; a[1] = b[1]; a[2] = b[2];
#define NOUSE (-1)

const double pi = std::acos(-1);

namespace Utilities {

    void ComputeNormal (vtkPoints *pts, vtkIdList *poly, double* n);
    double GetAngle (double *vA, double *vB, double *n);
    void FindPoints (vtkKdTreePointLocator *pl, const double *pt, vtkIdList *pts, double tol = 1e-6);
    void WriteVTK (const char *name, vtkPolyData *pd);

    class Pair {
    public:
        int f, g;
        Pair () {}
        Pair (int _f, int _g) : f(_f), g(_g) {}
        bool operator< (const Pair &other) const {
            return std::tie(f, g) < std::tie(other.f, other.g);
        }
        friend std::ostream& operator<< (std::ostream &out, const Pair &p) {
            out << "(" << p.f << ", " << p.g << ")";
            return out;
        }
    };

    typedef std::vector<int> IdsType;

    double Mod (int a, int b);
    double GetD (double *ptA, double *ptB);
    double GetD3 (double *ptA, double *ptB);
    double Normalize (double *p);
    void GetNormal (double pts[][3], double *n, const int num);

}

#endif
