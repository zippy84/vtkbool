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

void GetConvex (double pts[][2], const int num, std::vector<int> &convex) {
    // graham scan

    int s = 0;
    for (int i = 1; i < num; i++) {
        if (pts[i][1] < pts[s][1]) {
            s = i;
        }
    }

    std::cout << "s=" << s << std::endl;

    double angs[num];
    angs[s] = 0;

    std::vector<int> inds;

    for (int i = 0; i < num; i++) {
        inds.push_back(i);

        if (i == s) {
            continue;
        }

        double x[] = {
            pts[i][0]-pts[s][0],
            pts[i][1]-pts[s][1]
        };

        double l = sqrt(x[0]*x[0]+x[1]*x[1]);
        angs[i] = acos(x[0]/l);

        std::cout << i << ": " << l << " ang=" << angs[i]*180/pi << std::endl;
    }

    std::sort(inds.begin(), inds.end(), [&angs](const int &a, const int &b) {
        return angs[a] < angs[b];
    });

    std::vector<int>::const_iterator itr;
    for (itr = inds.begin(); itr != inds.end(); ++itr) {
        std::cout << *itr << ": " << angs[*itr] << std::endl;
    }

    std::copy(inds.begin(), inds.begin()+2, std::back_inserter(convex));

    int i = 2;
    while (i < num) {
        int a = *(convex.end()-2),
            b = *(convex.end()-1),
            c = inds[i];

        double *pA = pts[a],
            *pB = pts[b],
            *pC = pts[c];

        double d = (pB[0]-pA[0])*(pC[1]-pA[1])
            -(pC[0]-pA[0])*(pB[1]-pA[1]);

        std::cout << d << std::endl;

        if (d > 0) {
            convex.push_back(c);
            i++;
        } else {
            convex.pop_back();
        }
    }

    for (itr = convex.begin(); itr != convex.end(); ++itr) {
        std::cout << (itr-convex.begin() == 0 ? "M" : "L") << pts[*itr][0] << "," << pts[*itr][1] << " ";
    }
    std::cout << "Z" << std::endl;

}
