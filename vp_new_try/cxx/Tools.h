#ifndef __Tools_h
#define __Tools_h

#include <vector>
#include <cmath>
#include <memory>
#include <iostream>

#define NO_USE -1
#define PI std::acos(-1)

#define E 1e-5

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

double Normalize (double *v);
double GetAngle (double *vA, double *vB);
double Ld (double *a, double *b, double *c);
double Cross (double *a, double *b, double *c);

class D {
public:
    D (double *_s) {
        s[0] = _s[0];
        s[1] = _s[1];
    }
    D (double *_s, double _t1, double _t2) : D(_s) {
        t1 = _t1;
        t2 = _t2;
    }
    double s[2], t1, t2;
};

std::shared_ptr<D> Intersect (double *o, double *r, double *pA, double *pB);
std::shared_ptr<D> Intersect2 (double *oA, double *oB, double *pA, double *pB);

bool IsFrontfaced (double *r, double *a, double *b);
bool IsNear (double *a, double *b);
double GetT (double *a, double *b, double *c);

#endif
