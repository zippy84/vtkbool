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

#ifndef __Decomposer_h
#define __Decomposer_h

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cfloat>
#include <stack>
#include <sstream>
#include <set>
#include <deque>
#include <map>
#include <cassert>

#include "Utilities.h"
using Utilities::Pair;
using Utilities::Normalize;
using Utilities::GetNormal;

namespace Decomposer {

    class Base {
    public:
        Base (double pts[][3], const int num);
        double n[3], ei[3], ej[3], d;
    };

    void Transform (double *in, double *out, Base &base);
    void Transform (double in[][3], double out[][2], const int num, Base &base);

    bool Intersect (double *o, double *r, double *pA, double *pB, double *s, double *t = NULL, double *u = NULL);
    bool Intersect2 (double *oA, double *oB, double *pA, double *pB, double *s);

    class Vert {
    public:
        Vert (double *x, double *_pt, int _succ, int _ind = NOUSE) : succ(_succ), ind(_ind) {
            pt[0] = _pt[0];
            pt[1] = _pt[1];

            r[0] = pt[0]-x[0];
            r[1] = pt[1]-x[1];
            r[2] = 0;

            Normalize(r);

            phi = std::atan2(-r[1], -r[0]);

            if (phi < 0) {
                phi += 2*pi;
            }
            
        }

        int succ;
        double pt[2], r[3], phi;
        int ind;
    };

    class Vert2 {
    public:
        Vert2 (const Vert &v) {
            pt[0] = v.pt[0];
            pt[1] = v.pt[1];
            ind = v.ind;
        }
        double pt[2];
        int ind;
    };

    class Cut {
    public:
        Cut (int _i, double *_s, double _t, double _u) : i(_i), t(_t), u(_u) {
            s[0] = _s[0];
            s[1] = _s[1];
        }
        int i;
        double s[2], t, u;
    };

    double GetAng (double *vA, double *vB);

    class Bag {
    public:
        int f, g;
        double phi;
        Bag (int _f, int _g, double _phi) : f(_f), g(_g), phi(_phi) {}
    };

    void GetVisiblePoly (double pts[][2], const int num, std::vector<Vert2> &poly, int ind);

    bool check (double pts[][2], int num);
    
    typedef std::deque<Pair> SType;

    class Deq {
        SType S;
        int oF, oB;
    public:
        Deq () : oF(0), oB(0) {}
        void popF () {
            oF++;
        }
        void popB () {
            oB--;
        }
        void push (Pair &p) {
            S.erase(S.begin()+oB, S.end());
            S.push_back(p);
            oB++;
        }
        void push (Pair p) {
            S.erase(S.begin()+oB, S.end());
            S.push_back(p);
            oB++;
        }
        bool empty () {
            return oF == oB;
        }
        int size () {
            return oB-oF;
        }
        void restore () {
            oF = 0;
            oB = S.size();
        }
        const Pair& getF () {
            return *(S.begin()+oF);
        }
        const Pair& getB () {
            return *(S.begin()+oB-1);
        }
        const Pair& getPreF () {
            return *(S.begin()+oF+1);
        }
        const Pair& getPreB () {
            return *(S.begin()+oB-2);
        }
        void clear () {
            S.clear();
            oF = 0;
            oB = 0;
        }
        friend std::ostream& operator<< (std::ostream &out, const Deq &d) {
            out << "[";
            SType::const_iterator itr;
            for (itr = d.S.begin()+d.oF; itr != d.S.begin()+d.oB; ++itr) {
                std::cout << *itr << ", ";
            }
            out << "]";
            return out;
        }
    };

    class SubD {
    public:
        int w;
        Deq S;
        SubD () : w(9999) {}

        void update (int a, int b, int _w);
    };

    typedef std::map<Pair, SubD> PType;
    typedef std::vector<int> IdsType;

    typedef std::vector<IdsType> DecType;

    class Decompose {
        std::set<Pair> visVerts;
        std::vector<bool> refls;

        PType problems;

        double (*pts)[2];
        int num;

        bool isConcave;

        bool IsReflex (int i, int j, int k);
        void TypeA (int i, int j, int k);
        void TypeB (int i, int j, int k);

        void Recover (int i, int k);
        void Draw (int i, int k);

        std::vector<Pair> diags;

        void run ();

    public:
        Decompose (double _pts[][2], int _num) : pts(_pts), num(_num) {}
        void collect (DecType &dec);

    };

}

#endif
