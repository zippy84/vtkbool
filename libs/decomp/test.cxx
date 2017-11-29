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

    std::ifstream jn("../../vp/dev/complex.json");

    if (reader.parse(jn, doc)) {
        const Json::Value polys = doc["polys"];

        int i = 0;

        for (const Json::Value& p : polys) {
            //if (i == 2) {
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
            //}

            i++;

        }

    }

}

