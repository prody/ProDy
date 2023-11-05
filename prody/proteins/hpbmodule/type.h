/*
  Copyright (c) 2023, Xin Cao, Michelle H. Hummel, Bihua Yu, Evangelos A. Coutsias

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the “Software”), to deal
  in the Software without restriction, including without limitation the rights to
  use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
  the Software, and to permit persons to whom the Software is furnished to do so,
  subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#ifndef TYPE_H
#define TYPE_H
#include <cmath>
#include <vector>
#include <algorithm>
struct tri_type{
    int tri[3];
    int adjT[2];
    int op_ver[2];
    tri_type* next_ext;//list of tri on boundary of C(0)
    double frac;//=0, .5, 1 fraction to be used in volcomp
    //0 interior to C(0), .5 incident to only 1 tetra in C(0), 1 singular
    int id[2];//set id = void that it touches(only for
    //those on boundary of C(0) !id=0 on outside fringe, id>0 inside fringe,
    //id=-1 neither of the above
    int alf;//=1 if not in C(0), -1 exterior to C(0), -2 interior
    double x_char[4];
    tri_type() {
        alf = 0;
        x_char[0] = 0;
        x_char[1] = 0;
        x_char[2] = 0;
        x_char[3] = 0;
    }
};

struct edge_type{
    int edg[2];
    int adj_tri[20];
    int alf;
    int op_ver[20];
    int num_adj_tri;
    edge_type* next_ext;
    int disc;
    int tetra_ring[20];
    int num_ring;//total number of tetra (in or out of C(0)) that are in ring
    double x_char[3];
    double h1, h2;
    double t;
    edge_type() {
        x_char[0] = 0;
        x_char[1] = 0;
        x_char[2] = 0;
    }
};

struct ver_type{
    int adj_edg[100];
    int op_ver[100];
    int alf, num_adj_edg;
    int upverT[100];
    int last_upT;
    ver_type* next_ext;
};

struct ml_type{
    int tag;//1=T, -1=tri
    int tetra;//index
    tri_type* tri;
};

struct force_type{
    double force[3];
    int disc;
    int edge_ver;
    double direction[3];
    double correction[3];
};

struct tri2list_type{
    //just to be used when checking voids
    int tri[3];
    tri2list_type* next;
};

namespace operation {
    void cross_product(double b[3], double c[3], double c_prod[3]) {
        c_prod[0] = b[1]*c[2]-b[2]*c[1];
        c_prod[1] = -b[0]*c[2]+b[2]*c[0];
        c_prod[2] = b[0]*c[1]-b[1]*c[0];
        return;
    }

    double abs_cross_product(double b[3], double c[3]) {
        double abs_c_product;
        double c_prod[3];
        cross_product(b, c, c_prod);
        abs_c_product = sqrt(c_prod[0]*c_prod[0]+c_prod[1]*c_prod[1]+c_prod[2]*c_prod[2]);
        return abs_c_product;
    }
    /*
     * solving a system of linear equations Mx = b using Guassian elimination
     * x=r_0^2 - R^2, where r_0 is the coordinate of the characteristic point of a tetrahedron,
     * which is also the center of the circumsphere and
     * R is the radius of the circumsphere.
     * (r_b-r_0)^2 = (r_4-r_0)^2 => 2r_0(r_n-r_4)=r_n^2-r_4^2
     */
    void Minvb(double(&x)[4], double M[4][4], double b[4]) {
        double temp1[4];
        int maxrow;
        double l, rhs, temp2;
        for(int i = 0; i < 4; i ++) {
            x[i] = 0;
        }

        std::vector<double> absMcol;
        for(int i = 0; i < 4; ++ i) {
            for(int j = i; j < 4; ++ j) absMcol.push_back(std::abs(M[j][i]));
            maxrow = std::distance(absMcol.begin(), std::max_element(absMcol.begin(), absMcol.end() ) );
            maxrow = maxrow + i;
            absMcol.clear();
            for(int j = i; j < 4; ++ j) {
                temp1[j] = M[i][j];
                M[i][j] = M[maxrow][j];
                M[maxrow][j] = temp1[j];
            }
            temp2 = b[i];
            b[i] = b[maxrow];
            b[maxrow] = temp2;
            for(int j = i+1; j < 4; ++ j) {
                l = -M[j][i]/M[i][i];
                M[j][i] = 0;
                for(int k = 1+i; k < 4; ++ k)
                    M[j][k] = M[j][k] + l*M[i][k];
                b[j] = b[j]+l*b[i];
            }
        }

        x[3] = b[3]/M[3][3];
        for(int i = 2; i > -1; -- i){
            rhs = b[i];
            for(int j = i+1; j < 4; ++ j)
                rhs = rhs - M[i][j]*x[j];
            x[i] = rhs/M[i][i];
        }
    }
}

#endif // TYPE_H
