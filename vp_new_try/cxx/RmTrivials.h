#ifndef __RmTrivials_h
#define __RmTrivials_h

#include <vector>
#include <iostream>

#include "Tools.h"
#include "VisPoly.h"

class Vert2 {
public:
    Vert2 (int _i, double *_r, double _d, double _phi) : i(_i), d(_d), phi(_phi) {
        r[0] = _r[0];
        r[1] = _r[1];
    }
    int i;
    double r[2], d, phi;
};

typedef std::vector<Vert2> VertsType2;

void AlignPts (PolyType &poly, int ind);
void RemoveInterals (PolyType &poly, int skip);

enum class Src {
    NONE = 1,
    A = 2,
    B = 3
};

class Vert3 : public Point {
public:
    Vert3 (Point &_p) : Point(_p), rm(false), src(Src::NONE), t(0) {}
    Vert3 (double *_s, Src _src, double _t) : Point(_s), rm(false), src(_src), t(_t) {}

    double t;
    Src src;
    bool rm;

    int i;

    friend std::ostream& operator<< (std::ostream &out, const Vert3 &v) {
        out << "id: " << v.id
            << ", pt: [" << v.x << "," << v.y << "]"
            << ", t: " << v.t
            << ", src: " << static_cast<int>(v.src);
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
    IN = 0,
    OUT = 1
};

class Pair2 {
public:
    Pair2 (int _i, int _j, Dir _dir = Dir::FORWARD) : i(_i), j(_j), dir(_dir) {}

    int i, j;
    IdsType pocket;
    Dir dir;
    Side side;

    friend std::ostream& operator<< (std::ostream &out, const Pair2 &p) {
        out << "i: " << p.i
            << ", j: " << p.j
            << ", dir: " << static_cast<int>(p.dir)
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

void AddPocket (Pair2 &pair, VertsType3 &verts, double *rot, double d);
void RemovePockets (VertsType3 &verts, VertsType3 &good, double *rot, double d);

// nur diese funktion sollte verwendet werden
void RemoveTrivials (PolyType &poly, PolyType &res, int ind = 0);

#endif
