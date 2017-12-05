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
#include <sstream>
#include <cassert>

#include <json/json.h>
#include <json/reader.h>

#include "Tools.h"
#include "VisPoly.h"
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

    Json::Reader reader;

    std::ifstream jn("../dev/holes2.json");

    if (reader.parse(jn, doc)) {
        const Json::Value holes = doc["holes"];

        Merger m;

        for (const Json::Value& hole : holes) {
            PolyType poly;
            ToPoly(hole, poly);

            m.AddPoly(poly);

            assert(!TestCW(poly));
        }

        PolysType merged;

        m.GetMerged(merged);

        for (auto& p : merged) {
            assert(!TestCW(p));

            std::cout << GetAbsolutePath(p) << std::endl;
        }

    }
}
