#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cassert>

#include <json/json.h>
#include <json/reader.h>

#include "VisPoly.h"
#include "RmTrivials.h"

std::string GetAbsolutePath (const PolyType &poly) {
    std::stringstream path;

    for (const Point& p : poly) {
        path << "L" << p.x << "," << p.y << " ";
    }

    std::string svg = "M" + path.str().substr(1) + "Z";

    return svg;
}

bool TestCW (const PolyType &poly) {
    // http://mathworld.wolfram.com/PolygonArea.html

    int num = poly.size();

    double sum = 0;

    for (int i = 0; i < num; i++) {
        const Point &a = poly[i],
            &b = poly[(i+1)%num];
        sum += a.x*b.y-b.x*a.y;
    }

    return sum < 0;
}

int main (int argc, char *argv[]) {
    Json::Value doc;

    Json::Reader reader;

    std::ifstream jn("../../complex.json");

    if (reader.parse(jn, doc)) {
        const Json::Value polys = doc["polys"];

        int i = 0;

        for (const Json::Value& p : polys) {
            //if (i != 0) { continue; }

            PolyType poly;

            int j = 0;

            for (const Json::Value& pt : p) {
                poly.push_back(Point(pt[0].asDouble(), pt[1].asDouble(), j++));
            }

            int num = poly.size();

            std::cout << num << std::endl;

            for (int j = 1; j < num; j++) {
                poly[j].pt[0] += poly[j-1].pt[0];
                poly[j].pt[1] += poly[j-1].pt[1];
            }

            std::map<int, PolyType> all;

            for (int j = 0; j < num; j++) {
                PolyType poly2;
                std::copy(poly.begin(), poly.end(), std::back_inserter(poly2));

                // das polygon ist in clockwise order
                assert(TestCW(poly2));

                //if (j != 30) { continue; }

                //RemoveTrivials(poly2, all[j], j);

                PolyType poly3;
                RemoveTrivials(poly2, poly3, j);

                GetVisPoly(poly3, all[j]);

                for (auto& p : all[j]) {
                    std::cout << p << std::endl;
                }

                // das ergebnis ist in counterclockwise order
                assert(!TestCW(all[j]));
                
            }

            Json::Value data;

            for (const auto& itr : all) {
                /*
                Json::Value pts(Json::arrayValue);
                for (const auto& pt : itr.second) {

                    Json::Value _pt(Json::arrayValue);
                    _pt.append(pt.x);
                    _pt.append(pt.y);

                    pts.append(_pt);
                }

                data[std::to_string(itr.first)] = pts;
                */

                data[std::to_string(itr.first)] = GetAbsolutePath(itr.second);
            }

            Json::FastWriter writer;

            std::stringstream name;
            name << "../data_files/data_" << i++ << ".js";

            std::ofstream f(name.str());
            f << "var polys = " << writer.write(data) << ";";
            f.close();

        }

    }

}
