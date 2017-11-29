/*
   Copyright 2012-2018 Ronald Römer

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

#include <vector>
#include <iostream>
#include <set>
#include <deque>
#include <map>
#include <unordered_set>
#include <unordered_map>

#include "Tools.h"

class Vert6 : public Point {
public:
    Vert6 (Point &_p) : Point(_p), refl(false) {}
    bool refl;

    friend std::ostream& operator<< (std::ostream &out, const Vert6 &v) {
        out << "id: " << v.id
            << ", pt: [" << v.x << "," << v.y << "]"
            << ", refl: " << v.refl;
        return out;
    }
};

typedef std::vector<Vert6> VertsType6;

class SubP {
public:
    SubP () : w(99999) {}

    void AddPair (Pair _p, int _w);
    void RestoreS ();

    int w;
    std::deque<Pair> S, S_head, S_tail;
};

typedef std::vector<IdsType> DecResType;

// siehe boost
template<typename T> void hash_combine (std::size_t &seed, T const &v) {
    seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

struct p_h {
    std::size_t operator () (const Pair& p) const {
        std::size_t seed = 0;
        hash_combine(seed, p.f);
        hash_combine(seed, p.g);
        return seed;
    }
};

class Decomposer {
    PolyType poly;
    VertsType6 verts;

    int num;

    std::unordered_set<Pair, p_h> pairs;
    std::unordered_map<Pair, SubP, p_h> subs;

    bool IsRefl (int a, int b, int c);

    void Forw (int i, int j, int k);
    void Backw (int i, int j, int k);

    void Recover (int i, int k);
    void Collect (int i, int k);

    std::vector<Pair> diags;

public:
    Decomposer (PolyType &_poly);

    void GetDecomposed (DecResType &res);

};

void SimpleRmInternals (VertsType6 &verts);
void Scale (PolyType &poly);

#endif
