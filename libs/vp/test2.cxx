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

#include <fstream>
#include <iostream>
#include <fstream>
#include <iostream>
#include <cassert>
#include <string>

#include <json/json.h>
#include <json/reader.h>

#include "VisPoly.h"

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
    int s, t;

    stream >> s;

    stream.clear();
    stream.str(argv[2]);

    stream >> t;

    Json::Value doc;

    Json::CharReaderBuilder reader;

    std::ifstream jn("../dev/special.json");

    // special.json + ./test 3 0 führt zu einem assert
    // siehe kommentar zum assert(bnd != NO_USE); in RmTrivials.cxx

    std::string err;

    if (Json::parseFromStream(reader, jn, &doc, &err)) {
        const Json::Value polys = doc["polys"];

        int i = 0;

        for (const Json::Value& p : polys) {
            if (i == s) {
                PolyType poly;
                ToPoly(p, poly);

                int num = poly.size();

                std::map<int, PolyType> all;

                for (int j = 0; j < num; j++) {

                    //if (j != t) { continue; }

                    PolyType res;

                    try {
                        GetVisPoly_wrapper(poly, res, j);

                        assert(!TestCW(res));

                        all[j] = std::move(res);
                    } catch (std::exception &e) {
                        std::cerr << j << ": " << e.what() << std::endl;
                    }

                }

                Json::Value data;

                for (const auto& itr : all) {
                    data[std::to_string(itr.first)] = GetAbsolutePath(itr.second);
                }

                Json::StreamWriterBuilder writer;

                std::stringstream name;
                name << "../dev/data_files/special_" << i << ".js";

                std::ofstream f(name.str());
                f << "var pts = '" << GetAbsolutePath(poly)
                    << "'; var polys = " << Json::writeString(writer, data)
                    << ";";
                f.close();
            }

            i++;

        }

    }

}
