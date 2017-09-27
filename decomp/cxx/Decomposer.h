#include <vector>
#include <iostream>
#include <set>
#include <deque>
#include <map>

#include "Tools.h"

class Vert4 : public Point {
public:
    Vert4 (Point &_p) : Point(_p), refl(false) {}
    bool refl;

    friend std::ostream& operator<< (std::ostream &out, const Vert4 &v) {
        out << "id: " << v.id
            << ", pt: [" << v.x << "," << v.y << "]"
            << ", refl: " << v.refl;
        return out;
    }
};

typedef std::vector<Vert4> VertsType4;

class SubP {
public:
    SubP () : w(99999) {}

    void AddPair (Pair _p, int _w);
    void RestoreS ();

    int w;
    std::deque<Pair> S, S_head, S_tail;
};

typedef std::vector<IdsType> DecResType;

class Decomposer {
    PolyType poly;
    VertsType4 verts;
    
    int num;

    std::set<Pair> pairs;
    std::map<Pair, SubP> subs;

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

