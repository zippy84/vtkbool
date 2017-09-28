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

    }
}
