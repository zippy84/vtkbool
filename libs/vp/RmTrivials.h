/*
Copyright 2012-2019 Ronald Römer

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

#ifndef __RmTrivials_h
#define __RmTrivials_h

#include <vector>
#include <iostream>
#include <memory>

#include "Tools.h"
#include "VisPoly.h"

enum class Src {
    NONE = 1,
    A = 2,
    B = 3
};

class Vert3 : public Point {
public:
    Vert3 (Point &p) : Point(p), rm(false), src(Src::NONE), t(0) {}
    Vert3 (double *_s, Src _src, double _t) : Point(_s), rm(false), src(_src), t(_t) {}

    double t;
    Src src;
    bool rm;

    int i;

    friend std::ostream& operator<< (std::ostream &out, const Vert3 &v) {
        out << "id: " << v.id
            << ", pt: [" << v.x << "," << v.y << "]"
            << ", t: " << v.t
            << ", src: " << static_cast<int>(v.src)
            << ", rm: " << v.rm
            << ", tag: " << v.tag
            << ", i: " << v.i;
        return out;
    }

};

typedef std::vector<Vert3> VertsType3;

enum class Dir {
    FORWARD = 1,
    BACKWARD = 2,
    UNDEFINED = 3
};

enum class Side {
    NONE = 1,
    IN = 2,
    OUT = 3
};

class Pair2 {
public:
    Pair2 (int _i, int _j, Dir _dir = Dir::FORWARD) : i(_i), j(_j), dir(_dir), side(Side::NONE) {}

    int i, j;
    IdsType pocket;
    Dir dir;
    Side side;

    friend std::ostream& operator<< (std::ostream &out, const Pair2 &p) {
        out << "i: " << p.i
            << ", j: " << p.j
            << ", dir: " << (p.dir == Dir::FORWARD ? "FORWARD" : (p.dir == Dir::BACKWARD ? "BACKWARD" : "UNDEFINED"))
            << ", side: " << static_cast<int>(p.side)
            << ", pocket: [";

        for (int id : p.pocket) {
            out << id << ", ";
        }

        out << "]";
        return out;
    }

};

class Grp {
public:
    Grp (Dir _dir, IdsType _ids) : dir(_dir), ids(_ids) {}
    Dir dir;
    IdsType ids;
};

class TrivialRm {
    PolyType &poly;
    VertsType3 verts;

    Tracker &tr;

    int ind;

    Point x;

public:
    TrivialRm (PolyType &_poly, Tracker &_tr, int _ind, Point &_x) : poly(_poly), tr(_tr), ind(_ind), x(_x) {}

    void GetSimplified (PolyType &res);

private:
    void GetPocket (Pair2 &pair, IdsType &pocket);
    void AssignSide (Pair2 &pair, Src src);
    bool HasArea (const IdsType &pocket);

    void RemovePockets (VertsType3 &good, double *rot, double d, Src src);

    void RemoveRedundants (const VertsType3 &good);
};

#endif
