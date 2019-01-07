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

function copy (A) {
    var B = [];

    for (var i = 0; i < A.length; i++) {
        B.push([]);

        for (var j = 0; j < A[i].length; j++) {
            B[i].push(A[i][j]);
        }
    }

    return B;
}

function gaussEli (A, n) {
    var p = new Array(n-1);

    for (var i = 0; i < n-1; i++) {
        p[i] = 0;
    }

    for (var k = 0; k < n-1; k++) {
        var maxP = 0;

        for (var i = k; i < n; i++) {
            var q = 0;

            for (var j = k; j < n; j++) {
                q += Math.abs(A[i][j]);
            }

            var r = Math.abs(A[i][k])/q;

            if (r > maxP) {
                maxP = r;
                p[k] = i;
            }
        }

        if (p[k] != k) {
            for (var j = 0; j < n; j++) {
                var tmp = A[k][j];
                A[k][j] = A[p[k]][j];
                A[p[k]][j] = tmp;
            }
        }

        for (var i = k+1; i < n; i++) {
            A[i][k] /= A[k][k];

            for (var j = k+1; j < n; j++) {
                A[i][j] -= A[i][k]*A[k][j];
            }
        }
    }

    return p; // A wurde verändert
}

function gaussSol (B, c, p, n) {
    for (var k = 0; k < n-1; k++) {
        if (p[k] != k) {
            var tmp = c[k];
            c[k] = c[p[k]];
            c[p[k]] = tmp;
        }
    }

    for (var i = 0; i < n; i++) {
        for (var j = 0; j < i; j++) {
            c[i] -= B[i][j]*c[j];
        }
    }

    var x = new Array(n);

    for (var i = 0; i < n; i++) {
        x[i] = 0;
    }

    for (var i = n-1; i > -1; i--) {
        var s = c[i];

        for (var k = i+1; k < n; k++) {
            s -= B[i][k]*x[k];
        }

        x[i] = s/B[i][i];
    }

    return x;
}

function gauss (A, b) {
    var n = A[0].length;

    var B = copy(A);

    var p = gaussEli(B, n);

    var xApprox = gaussSol(B, b.slice(), p, n);

    var res = new Array(n);

    for (var i = 0; i < n; i++) {
        res[i] = 0;

        for (var j = 0; j < n; j++) {
            res[i] += A[i][j]*xApprox[j];
        }

        res[i] -= b[i];

        res[i] *= -1;

    }

    var corr = gaussSol(B, res.slice(), p, n);

    var x = xApprox.slice();

    for (var i = 0; i < n; i++) {
        x[i] += corr[i];
    }

    return x;

}

