#include <cassert>
#include <iostream>
#include <fstream>

#include <json/json.h>
#include <json/reader.h>

#include "Decomposer.h"
#include "Tools.h"

int main (int argc, char *argv[]) {
    Json::Value doc;

    Json::Reader reader;

    std::ifstream jn("../../vp_new_try/dev/complex.json");

    if (reader.parse(jn, doc)) {
        const Json::Value polys = doc["polys"];

        int i = 0;

        for (const Json::Value& p : polys) {
            /*
            if (i++ != 1) {
                continue;
            }
            */

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

            assert(TestCW(poly));

            Decomposer d(poly);

            DecResType decs;
            d.GetDecomposed(decs);

            for (auto& dec : decs) {
                PolyType p;

                for (int id : dec) {
                    p.push_back(poly[id]);
                }

                std::cout << GetAbsolutePath(p) << std::endl;

                assert(TestCW(p));
            }

        }

    }

}

