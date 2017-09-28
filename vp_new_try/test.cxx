#include <fstream>
#include <iostream>

#include <json/json.h>
#include <json/reader.h>

#include <VisPoly.h>

int main (int argc, char *argv[]) {
    Json::Value doc;

    Json::Reader reader;

    std::ifstream jn("../dev/polys.json");

    if (reader.parse(jn, doc)) {
        const Json::Value polys = doc["polys"];

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

            GetVisPoly(poly);

        }

    }

}
