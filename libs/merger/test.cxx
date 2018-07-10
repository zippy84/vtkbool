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

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>

#include <json/json.h>
#include <json/reader.h>

#include "inja.hpp"
#include "nlohmann/json.hpp"

#include "Tools.h"
#include "Merger.h"

void ToPoly (const Json::Value& pts, PolyType &poly) {

    for (const Json::Value& pt : pts) {
        poly.push_back(Point(pt[0].asDouble(), pt[1].asDouble()));
    }

    for (int j = 1; j < poly.size(); j++) {
        poly[j].pt[0] += poly[j-1].pt[0];
        poly[j].pt[1] += poly[j-1].pt[1];
    }
}

int main (int argc, char *argv[]) {
    Json::Value doc;

    Json::CharReaderBuilder reader;

    std::ifstream jn("../dev/holes2.json");

    std::string err;

    if (Json::parseFromStream(reader, jn, &doc, &err)) {
        const Json::Value holes = doc["holes"];

        inja::Environment env;
        env.set_element_notation(inja::ElementNotation::Dot);

        nlohmann::json data;

        Merger m;

        for (const Json::Value& hole : holes) {
            PolyType poly;
            ToPoly(hole, poly);

            m.AddPoly(poly);

            assert(!TestCW(poly));
        }

        PolysType merged;

        m.GetMerged(merged);

        std::vector<double> xs, ys;

        for (auto& p : merged) {
            assert(!TestCW(p));

            data["data"].push_back({{ "path", GetAbsolutePath(p) }});

            for (auto &pt : p) {
                xs.push_back(pt.x);
                ys.push_back(pt.y);
            }
        }

        double xMin = *std::min_element(xs.begin(), xs.end()),
            xMax = *std::max_element(xs.begin(), xs.end());

        double yMin = *std::min_element(ys.begin(), ys.end()),
            yMax = *std::max_element(ys.begin(), ys.end());

        data["width"] = std::abs(xMax-xMin);
        data["height"] = std::abs(yMax-yMin);

        data["x"] = -xMin;
        data["y"] = -yMin;

        env.write("../dev/template.svg", data, "../dev/output2.svg");

    }
}
