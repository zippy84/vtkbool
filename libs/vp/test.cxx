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

#include <fstream>
#include <iostream>
#include <cassert>

#include <json/json.h>
#include <json/reader.h>

#include "inja.hpp"
#include "nlohmann/json.hpp"

#include "VisPoly.h"

int main (int argc, char *argv[]) {
    Json::Value doc;

    Json::CharReaderBuilder reader;

    std::ifstream jn("../dev/polys.json");

    std::string err;

    if (Json::parseFromStream(reader, jn, &doc, &err)) {
        const Json::Value polys = doc["polys"];

        inja::Environment env;
        env.set_element_notation(inja::ElementNotation::Dot);

        nlohmann::json data;

        int idx = 0;

        for (const Json::Value& p : polys) {
            PolyType poly;

            int i = 0;

            for (const Json::Value& pt : p) {
                poly.push_back(Point(pt[0].asDouble(), pt[1].asDouble(), i++));
            }

            std::cout << poly.size() << std::endl;

            for (int j = 1; j < poly.size(); j++) {
                poly[j].pt[0] += poly[j-1].pt[0];
                poly[j].pt[1] += poly[j-1].pt[1];
            }

            std::reverse(poly.begin(), poly.end());

            std::rotate(poly.begin(), poly.end()-1, poly.end());

            assert(TestCW(poly));

            PolyType res;

            GetVisPoly_wrapper(poly, res, 0);

            //std::cout << "Res " << GetAbsolutePath(res) << std::endl;

            int row = static_cast<int>(idx/4),
                col = static_cast<int>(idx%4);

            data["data"].push_back({{ "path", GetAbsolutePath(res) }, { "pts", GetAbsolutePath(poly) }, { "x", col*200 }, { "y", row*300 }});

            idx++;

        }

        env.write("../dev/template.svg", data, "../dev/output.svg");

    }

}
