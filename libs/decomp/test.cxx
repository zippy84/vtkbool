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

#include <cassert>
#include <iostream>
#include <fstream>

#include <json/json.h>
#include <json/reader.h>

#include "inja.hpp"
#include "nlohmann/json.hpp"

#include "Decomposer.h"
#include "Tools.h"

int Point::_tag = 0;

void ToPoly (const Json::Value& pts, PolyType &poly) {
    int i = 0;

    for (const Json::Value& pt : pts) {
        poly.push_back(Point(pt[0].asDouble(), pt[1].asDouble(), i++));
    }

    for (int j = 1; j < poly.size(); j++) {
        poly[j].pt[0] += poly[j-1].pt[0];
        poly[j].pt[1] += poly[j-1].pt[1];
    }
}

int main (int argc, char *argv[]) {
    std::istringstream stream(argv[1]);
    int t;

    stream >> t;

    Json::Value doc;

    Json::CharReaderBuilder reader;

    std::ifstream jn("../../vp/dev/special.json");

    std::string err;

    if (Json::parseFromStream(reader, jn, &doc, &err)) {
        const Json::Value polys = doc["polys"];

        inja::Environment env;
        env.set_element_notation(inja::ElementNotation::Dot);

        int i = 0;

        for (const Json::Value& p : polys) {
            if (i == t) {
                PolyType poly;

                ToPoly(p, poly);

                assert(TestCW(poly));

                nlohmann::json data;

                Ext ext;
                GetExt(poly, ext);

                data["width"] = std::abs(ext.maxX-ext.minX);
                data["height"] = std::abs(ext.maxY-ext.minY);

                data["x"] = -ext.minX;
                data["y"] = -ext.minY;

                data["poly"] = GetAbsolutePath(poly);

                Decomposer d(poly, 1.);

                DecResType decs;
                d.GetDecomposed(decs);

                for (auto& dec : decs) {
                    PolyType p;

                    for (int id : dec) {
                        p.push_back(poly[id]);
                    }

                    data["data"].push_back({{ "path", GetAbsolutePath(p) }});

                    assert(TestCW(p));
                }

                std::stringstream name;
                name << "../dev/res/special_" << i << ".svg";

                env.write("../dev/template.svg", data, name.str());
            }

            i++;

        }

    }

}
