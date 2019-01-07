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
    Json::Value doc;

    Json::CharReaderBuilder reader;

    std::ifstream jn("../dev/holes.json");

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

        Ext ext;

        for (auto& p : merged) {
            GetExt(p, ext);
            data["data"].push_back({{ "path", GetAbsolutePath(p) }});
            assert(!TestCW(p));
        }

        data["width"] = std::abs(ext.maxX-ext.minX);
        data["height"] = std::abs(ext.maxY-ext.minY);

        data["x"] = -ext.minX;
        data["y"] = -ext.minY;

        env.write("../dev/template.svg", data, "../dev/output.svg");

    }
}
