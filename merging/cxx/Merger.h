#ifndef __Merger_h
#define __Merger_h

#include <vector>
#include <map>

#include "VisPoly.h"

typedef std::vector<PolyType> PolysType;

typedef std::map<Pair, double> ConsType;
typedef std::map<int, ConsType> ResType;

class G {
public:
	G (double _d, Pair _con) : d(_d), con(_con) {}
	double d;
	Pair con;
};

class Merger {
    PolysType polys;

    void Merge (PolysType &group, PolyType &merged);

public:
    Merger () {}
    void AddPoly (PolyType &poly);
    void GetMerged (PolysType &res);
};

#endif
