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

#ifndef FIND_VOL_MOD_H
#define FIND_VOL_MOD_H

#include "type.h"
#include "DT.h"

#define Pi 3.141592653589793
#define Pi2 6.283185307179586
#define Pi4d3 4.188790204786391
#define Pi_half 1.57079632679490
#define One_third 0.333333333333333333
#define Two_thirds 0.666666666666666667
#define One_sixth 0.166666666666667

void voltetra_V(double& V_total, int& num_T2, double** coords, int** T2);

void voltetra_ext(int* idT, double* V_atom, double* S_lag_int, edge_type* edgelist, tri_type* reftri, double** coords, double** char_pt_T);

void unravel(int& num_ring_in, int* ring_in, int* adj_tri, int& num_ring_out, int* ring_out, int& bit, int& previous, int* idT);

void voltetra(int& num_int_edge, int* int_edge_list, edge_type* refedg, double** coords, double** char_pt_T, double* V_atom, double* S_lag_int, double* lift, std::vector<std::vector<double> > &surfint, int& idx);

void compare_vol(int& num_coord, double& V_total, double* V_atom);

void voltri(double& V_total, double* V_atom, double* S_atom, double* S_lag_int, double** coords, double* weights, tri_type* trilist, double *S_lag2);

void voledg(double** coords, double* weights, edge_type* edgelist, int* idT, double& V_total, double* V_atom, double* S_atom, double* S_lag_int, int **T, std::vector<std::vector<double> > &surfint, int& idx);

void volver(int& num_ver_ext, int& maxupver, int* ver_ext, ver_type* verlist, double* weigths, double** coords, double& V_total, double* V_atom, double* S_atom, double* S_lag_int, int *idT, int **T);

void classify_ver(int& num_coord, int& num_ver_ext, double* weights, ver_type* verlist, int* ver_ext, edge_type *refedg, double **coords);

void classify_edg(int& num_edg_ext, int& num_int_edge, int& num_edg, edge_type* edgelist, edge_type* refedg, double** coords, double* lift, double* weights, tri_type* reftri, int* int_edg_list);

void classify_tri(int* alfT, int& num_tri_ext, int &num_edg, int& num_tri, tri_type* reftri, tri_type* trilist, double** coords, double* lift, double* weights);

void extract_edg_tri(int** atri, int** aedg, int num_T, int num_coord, ver_type* verlist, tri_type* reftri, edge_type* refedg, int** T, int** T_nbrs, int &num_tri, int& num_edg);

void tetra_ring(int& edg_ind, int edg_num_int, int T_ind, edge_type* refedg, int** T, int** Tnbr, int** atri, int& num_tri, tri_type* reftri, int** aedg);

void mark_edg_tri(int &edg_ind, int &tag, int &o_ver, int &tetra_ind, int &edg_int, edge_type *refedg, int **T, tri_type *reftri, int **aedg, int **atri);

void classify_T(double** coords, double* lift, int& num_T, int** T, ver_type* verlist, int& maxupver, double** char_pt_T, int* alfT, int** T2, int &num_T2);

bool findele(std::string* res_sel, int num_sel, std::string name1)
{
    for (int i=0; i < num_sel; i++) {
        if (name1 == res_sel[i]) return true;
    }
    return false;
}

/*
 * special index array:
 * alfT[0:num_T]
 * idT[-1:num_T]
 *
 */

void find_vol_intersection(int num_coord, double** coords, double* weights, std::string* res, double* allvalues, double* V_atom, double* S_atom, double* S_lag_int, int& PR)
{
    int num_T, num_T2, num_int_edge, num_tri, num_edg, num_ver_ext, num_tri_ext, num_edg_ext;
    int maxupver = 100;
    int maxupedg = 20;
    int** T = NULL;
    int** T_nbrs = NULL;
    int* ver_ext = new int[num_coord];
    double* lift = new double[num_coord];
    double V_total = 0;
    int maxupver_actual;
    tri_type* trilist;
    edge_type* edgelist;
    ver_type* verlist = new ver_type[num_coord];
    double ** dt_coords = new double*[num_coord];
    double * dt_weights = new double[num_coord];
    double * S_lag1 = new double[num_coord];
    double * S_lag2 = new double[num_coord];
    for(int i = 0; i < num_coord ; ++ i) {
        dt_coords[i] = new double[3];
        dt_weights[i] = weights[i];
        ver_ext[i] = 0;
        S_lag2[i] = 0;
        for(int j = 0; j < 3; ++ j)
            dt_coords[i][j] = coords[i][j];
    }

    Delaunay_Tetrahedrize(num_coord, dt_coords, dt_weights, num_T, T, T_nbrs);
    tri_type* reftri = new tri_type[4*num_T];
    edge_type* refedg = new edge_type[12*num_T];
    int* alfT = new int[num_T+1]; // integer, dimension(0:num_T) :: alfT !just need sign of alfT;
    int* idT = new int[num_T+2]; //integer, dimension(-1:num_T) :: idT;
    int* int_edge_list = new int[6*num_T];
    double** char_pt_T = new double*[num_T];
    int** T2 = new int*[num_T];
    int** atri = new int*[num_T];
    int** aedg = new int*[num_T];
    for(int i = 0; i < num_T; ++ i) {
        char_pt_T[i] = new double[3];
        T2[i] = new int[4];
        atri[i] = new int[4];
        aedg[i] = new int[6];
    }

    maxupver_actual=sizeof(verlist[0].op_ver)/sizeof(int);

    if (maxupver_actual != maxupver) {
        std::cerr << "maxupver values are not equal, stop" << std::endl;
        std::cerr << maxupver_actual << " " << maxupver << std::endl;
        exit(0);
    }

    for(int i = 0; i < num_coord; ++ i) {
        lift[i] = 0;
        for(int j = 0; j < 3; ++ j)
            lift[i] += coords[i][j]*coords[i][j];
        lift[i] -= weights[i];
        verlist[i].last_upT = 0;
    }

    classify_T(coords, lift, num_T, T, verlist, maxupver, char_pt_T, alfT, T2, num_T2);

    tri_type tep;
    trilist = &tep;
    edge_type etep;
    edgelist = &etep;

    extract_edg_tri(atri, aedg, num_T, num_coord, verlist, reftri, refedg, T, T_nbrs, num_tri, num_edg);

    classify_tri(alfT, num_tri_ext, num_edg, num_tri, reftri, trilist, coords, lift, weights);
    classify_edg(num_edg_ext, num_int_edge, num_edg, edgelist, refedg, coords, lift, weights, reftri, int_edge_list);
    classify_ver(num_coord, num_ver_ext, weights, verlist, ver_ext, refedg, coords);

    idT[1] = 1;
    for (int i = 1; i < num_T+1; ++ i) idT[i+1] = alfT[i];

    std::vector<std::vector<double> > surfint;
    int edg_count = 0;

    volver(num_ver_ext, maxupver, ver_ext, verlist, weights, coords, V_total, V_atom, S_atom, S_lag_int, idT, T);
    voledg(coords, weights, edgelist, idT, V_total, V_atom, S_atom, S_lag_int, T, surfint, edg_count);
    voltri(V_total, V_atom, S_atom,  S_lag_int, coords, weights, trilist, S_lag2);
    voltetra_V(V_total, num_T2, coords, T2);
    voltetra(num_int_edge, int_edge_list, refedg, coords, char_pt_T, V_atom, S_lag_int, lift, surfint, edg_count);
    voltetra_ext(idT, V_atom, S_lag_int, edgelist, reftri, coords, char_pt_T);

    if (PR == 1) {
        std::string res_sel[7] = {"ALA", "ILE", "LEU", "MET", "PHE", "TRP", "VAL"};
        std::stable_sort(surfint.begin(), surfint.end());

        PR = 0;
        for (int i = 0; i < edg_count; ++ i) {
            std::string name1 = res[ int(surfint[i][0])-1 ];
            std::string name2 = res[ int(surfint[i][1])-1 ];
            std::string cres1 = name1.substr(0, 3);
            std::string cres2 = name2.substr(0, 3);
            if (cres1 != cres2) {
                if (findele(res_sel, 7, cres1) || findele(res_sel, 7, cres2)) {
                    allvalues[3*PR]   = surfint[i][0];
                    allvalues[3*PR+1] = surfint[i][1];
                    allvalues[3*PR+2] = surfint[i][2];
                    PR = PR + 1;
                }
            }
        }
    }

    for(int i = 0; i < num_coord ; ++ i) delete[] dt_coords[i];
    delete[] dt_coords;
    delete[] dt_weights;
    surfint.clear();

    for(int i = 0; i < num_T; ++ i){
      delete[] T[i];
      delete[] T_nbrs[i];
      delete[] T2[i];
      delete[] atri[i];
      delete[] aedg[i];
    }

    delete[] T;
    delete[] T_nbrs;
    delete[] T2;
    delete[] atri;
    delete[] aedg;

    T = NULL;
    T_nbrs = NULL;
    delete[] verlist;
    delete[] reftri;
    delete[] refedg;
    for(int i = 0; i < num_T; ++ i) {
        delete[] char_pt_T[i];
    }

    delete[] char_pt_T;
    delete[] ver_ext;
    delete[] lift;
    delete[] alfT;
}


void voltetra_V(double& V_total, int& num_T2, double** coords, int** T2)
{
    double vol_tetra;
    for (int i = 0; i < num_T2; ++ i) {
        double p1[3] = {coords[T2[i][0]-1][0], coords[T2[i][0]-1][1], coords[T2[i][0]-1][2]};
        double a[3] = {coords[T2[i][1]-1][0]-p1[0], coords[T2[i][1]-1][1]-p1[1], coords[T2[i][1]-1][2]-p1[2]};
        double b[3] = {coords[T2[i][2]-1][0]-p1[0], coords[T2[i][2]-1][1]-p1[1], coords[T2[i][2]-1][2]-p1[2]};
        double c[3] = {coords[T2[i][3]-1][0]-p1[0], coords[T2[i][3]-1][1]-p1[1], coords[T2[i][3]-1][2]-p1[2]};
        double axb[3] = {a[1]*b[2]-a[2]*b[1], -a[0]*b[2]+a[2]*b[0], a[0]*b[1]-a[1]*b[0]};
        double dot_pro = 0;
        for (int j = 0; j < 3; ++ j) {
            dot_pro += c[j] * axb[j];
        }
        if (dot_pro < 0) dot_pro = -dot_pro;
        vol_tetra = dot_pro * One_sixth;
        V_total += vol_tetra;
    }
}

void voltetra_ext(int* idT, double* V_atom, double* S_lag_int, edge_type* edgelist, tri_type* reftri, double** coords, double** char_pt_T)
{
    int edg[2], num_ring, e1, e2, tetra_ring[20], ring2[40], adj_tri[20], num_ring2, previous, curr_ring, bit;
    double c_edg = 0;
    double c_prod_2 = 0;
    double surf, a[3], b[3], c1[3], area_tri, h1, h2, vol, t, c_prod[3], edge_vec[3];
    edge_type* current = edgelist->next_ext;
    while (current) {
        bit = 1;
        edg[0] = current->edg[0];
        edg[1] = current->edg[1];
        e1 = edg[0];
        e2 = edg[1];
        num_ring = current->num_ring;
        for (int i = 0; i < num_ring; ++ i) {
            tetra_ring[i] = current->tetra_ring[i];
            adj_tri[i] = current->adj_tri[i];
        }
        unravel(num_ring, tetra_ring, adj_tri, num_ring2, ring2, bit, previous, idT);
        if (bit == 1) {
            edge_vec[0] = coords[e2-1][0] - coords[e1-1][0];
            edge_vec[1] = coords[e2-1][1] - coords[e1-1][1];
            edge_vec[2] = coords[e2-1][2] - coords[e1-1][2];
        } else {
            exit(0);
            /*
             *?????
            */
            edge_vec[0] = coords[e1-1][0] - coords[e2-1][0];
            edge_vec[1] = coords[e1-1][1] - coords[e2-1][1];
            edge_vec[2] = coords[e1-1][2] - coords[e2-1][2];
        }
        c1[0] = current->x_char[0];
        c1[1] = current->x_char[1];
        c1[2] = current->x_char[2];
        surf = 0;
        previous = -1;
        for (int j = 0; j < num_ring2; ++ j) {
            curr_ring = ring2[j];
            if (curr_ring < 0) {
                if (previous == 1) {
                    for (int k = 0; k < 3; ++ k) {
                        a[k] = b[k];
                        b[k] = reftri[-ring2[j]-1].x_char[k]-c1[k];
                    }
                    operation::cross_product(b, a, c_prod);
                    c_edg = 0;
                    c_prod_2 = 0;
                    for (int k = 0; k < 3; ++ k) {
                        c_edg += c_prod[k] * edge_vec[k];
                        c_prod_2 += c_prod[k] * c_prod[k];
                    }
                    if (c_edg >= 0) area_tri = 0.5 * sqrt(c_prod_2);
                    else area_tri = -0.5*sqrt(c_prod_2);
                    surf += area_tri;
                } else {
                    for (int k = 0; k < 3; ++ k) {
                        b[k] = reftri[-ring2[j]-1].x_char[k]-c1[k];
                    }
                }
                previous = -1;
            } else {
                for (int k = 0; k < 3; ++ k) {
                    a[k] = b[k];
                    b[k] = char_pt_T[ring2[j]-1][k]-c1[k];
                }
                operation::cross_product(b, a, c_prod);
                c_edg = 0;
                c_prod_2 = 0;
                for (int k = 0; k < 3; ++ k) {
                    c_edg += c_prod[k] * edge_vec[k];
                    c_prod_2 += c_prod[k] * c_prod[k];
                }
                if (c_edg >= 0) area_tri = 0.5 * sqrt(c_prod_2);
                else area_tri = -0.5*sqrt(c_prod_2);
                surf += area_tri;
                previous = 1;
            }
        }
        if (num_ring2 != 0) {
            h1 = current->h1;
            S_lag_int[e1-1] += surf;
            vol = One_third*surf*h1;
            t = current->t;
            V_atom[e1-1] += vol;

            h2 = current->h2;
            S_lag_int[e2-1] += surf;
            vol = One_third*surf*h2;
            V_atom[e2-1] += vol;
        }
        current = current->next_ext;
    }
}

void unravel(int& num_ring_in, int* ring_in, int* adj_tri, int& num_ring_out, int* ring_out, int& bit, int& previous, int* idT)
{
    int mark, tri_num, adj_tri_temp[20], ring_temp[40];
    for (int i = 0; i < 2*num_ring_in; ++ i) ring_out[i] = 0;
    if (ring_in[num_ring_in-1] == 0) {
        ring_temp[0] = 0;
        for (int k = num_ring_in-1; k >= 1; -- k) {
            adj_tri_temp[num_ring_in-k-1] = adj_tri[k-1];
            ring_temp[num_ring_in-k] = ring_in[k-1];
            if (ring_in[k-1] == 0) {
                mark = k;
                break;
            }
        }
        for (int k = 1; k <= mark-1; ++ k) {
            adj_tri_temp[num_ring_in-mark+k-1] = adj_tri[k-1];
            ring_temp[num_ring_in-mark+k-1] = ring_in[k-1];
        }
        ring_temp[num_ring_in-1] = 0;
        for (int i = 0; i < num_ring_in-1; ++ i) adj_tri[i] = adj_tri_temp[i];
        for (int i = 0; i < num_ring_in; ++ i) ring_in[i] = ring_temp[i];
    }
    num_ring_out = 0;
    previous = 1;
    //idT[-1:num_T]
    if (idT[ring_in[0]+1] != -1) {
        for (int k = 1; k < num_ring_in; ++ k) {
            if (idT[ring_in[k]+1] == -1) {
                if (previous == 1) {
                    tri_num = adj_tri[k-1];
                    num_ring_out += 1;
                    ring_out[num_ring_out-1] = -tri_num;
                    num_ring_out += 1;
                    ring_out[num_ring_out-1] = ring_in[k];
                } else {
                    num_ring_out += 1;
                    ring_out[num_ring_out-1] = ring_in[k];
                }
                previous = -1;
            } else {
                if (previous == -1) {
                    tri_num = adj_tri[k-1];
                    num_ring_out += 1;
                    ring_out[num_ring_out-1] = -tri_num;
                }
                previous = 1;
            }
        }
        if (idT[ring_in[num_ring_in-1]+1] == -1) {
            num_ring_out += 1;
            tri_num = adj_tri[num_ring_in-1];
            ring_out[num_ring_out-1] = -tri_num;
        }
    } else {
            for (int k = num_ring_in; k >=2; -- k) {
                if (idT[ring_in[k-1]+1] != -1) {
                    mark = k;
                    break;
                }
            }
            tri_num = adj_tri[mark-1];
            ring_out[0] = -tri_num;
            for (int k = mark; k < num_ring_in; ++ k) ring_out[k-mark+1] = ring_in[k];
            num_ring_out = num_ring_in - mark + 2;
            ring_out[num_ring_out-1] = ring_in[0];
            previous = -1;
            for (int k = 1; k < mark; ++ k) {
                if (idT[ring_in[k]+1] == -1) {
                    if (previous == 1) {
                        tri_num = adj_tri[k-1];
                        num_ring_out += 1;
                        ring_out[num_ring_out-1] = -tri_num;
                        num_ring_out += 1;
                        ring_out[num_ring_out-1] = ring_in[k];
                    } else {
                        num_ring_out += 1;
                        ring_out[num_ring_out-1] = ring_in[k];
                    }
                    previous = -1;
                } else {
                    if (previous == -1) {
                        tri_num = adj_tri[k-1];
                        num_ring_out += 1;
                        ring_out[num_ring_out-1] = -tri_num;
                    }
                    previous = 1;
                }
            }
    }
}

void voltetra(int& num_int_edge, int* int_edge_list, edge_type* refedg, double** coords, double** char_pt_T, double* V_atom, double* S_lag_int, double* lift, std::vector<std::vector<double> > &surfint, int& idx)
{
    double surf, area_tri, t, h1, h2, vol, c_prod2, c1[3], c2[3], c3[3], a[3], b[3], x[3];
    int e1, e2, num_ring, edge_ind, tetra_ring[20];
    for (int k =0; k < num_int_edge; ++ k) {
        surf = 0;
        edge_ind = int_edge_list[k];
        int edg[2] = {refedg[edge_ind-1].edg[0], refedg[edge_ind-1].edg[1]};
        e1 = edg[0];
        e2 = edg[1];
        num_ring = refedg[edge_ind-1].num_ring;
        for(int i = 0; i < 20; ++ i) tetra_ring[i] = refedg[edge_ind-1].tetra_ring[i];
        double p1[4] = {coords[e1-1][0], coords[e1-1][1], coords[e1-1][2], lift[e1-1]};
        double p2[4] = {coords[e2-1][0], coords[e2-1][1], coords[e2-1][2], lift[e2-1]};
        double p12[4] = {p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2], p1[3]-p2[3]};
        t = p12[3] - 2*(p2[0]*p12[0] + p2[1]*p12[1] + p2[2]*p12[2]);
        t = t/2.0/(p12[0]*p12[0]+p12[1]*p12[1]+p12[2]*p12[2]);
        x[0] = p12[0]*t + p2[0];
        x[1] = p12[1]*t + p2[1];
        x[2] = p12[2]*t + p2[2];
        if (num_ring < 3) {
            std::cerr << "Error with number fo tetra in ring: stop" << std::endl;
            exit(0);
        }
        for (int i = 0; i < 3; ++ i) {
            c1[i] = char_pt_T[tetra_ring[0]-1][i];
            c2[i] = char_pt_T[tetra_ring[1]-1][i];
            c3[i] = char_pt_T[tetra_ring[2]-1][i];
            a[i] = c1[i] - c3[i];
            b[i] = c1[i] - c2[i];
        }
        double c_prod[3] = {a[1]*b[2]-a[2]*b[1], -a[0]*b[2]+a[2]*b[0], a[0]*b[1]-a[1]*b[0]};
        c_prod2 = 0;
        for (int i = 0; i < 3; ++ i) c_prod2 += c_prod[i] * c_prod[i];
        area_tri = 0.5*sqrt(c_prod2);
        surf += area_tri;
        for (int i = 3; i < num_ring; ++ i) {
            c2[0] = c3[0]; c2[1] = c3[1]; c2[2] = c3[2];
            c3[0] = char_pt_T[tetra_ring[i]-1][0];
            c3[1] = char_pt_T[tetra_ring[i]-1][1];
            c3[2] = char_pt_T[tetra_ring[i]-1][2];
            a[0] = c1[0] - c3[0]; a[1] = c1[1] - c3[1]; a[2] = c1[2] - c3[2];
            b[0] = c1[0] - c2[0]; b[1] = c1[1] - c2[1]; b[2] = c1[2] - c2[2];
            c_prod[0] = a[1]*b[2] - a[2]*b[1];
            c_prod[1] = -a[0]*b[2] + a[2]*b[0];
            c_prod[2] = a[0]*b[1] - a[1]*b[0];
            c_prod2 = 0;
            for (int j = 0; j < 3; ++ j) c_prod2 += c_prod[j] * c_prod[j];
            area_tri = 0.5*sqrt(c_prod2);
            surf += area_tri;
        }
        double x_p12 = 0;
        for (int i = 0; i < 3; ++ i) x_p12 += (x[i]-p1[i])*(x[i]-p1[i]);
        h1 = sqrt(x_p12);
        S_lag_int[e1-1] = S_lag_int[e1-1] + surf;
        vol = One_third*surf*h1;
        if (t <= 1.0) {
            V_atom[e1-1] += vol;
        } else {
            V_atom[e1-1] -= vol;
        }

        double x_p22 = 0;
        for (int i = 0; i < 3; ++ i) x_p22 += (x[i]-p2[i])*(x[i]-p2[i]);
        h2 = sqrt(x_p22);
        S_lag_int[e2-1] = S_lag_int[e2-1] + surf;
        vol = One_third*surf*h2;
        if (t >= 0) {
            V_atom[e2-1] += vol;
        } else {
            V_atom[e2-1] -= vol;
        }

        std::vector<double> csurf;
        if (e1 > e2) {
            csurf.push_back(e2);
            csurf.push_back(e1);
        } else {
            csurf.push_back(e1);
            csurf.push_back(e2);
        }
        csurf.push_back(surf);
        surfint.push_back(csurf);
        csurf.clear();
        idx = idx+1;
    }
}

void compare_vol(int& num_coord, double& V_total, double* V_atom)
{
     double sum_LV = 0;
     for (int i = 0; i < num_coord; ++ i) sum_LV += V_atom[i];
     if ((sum_LV-V_total)/num_coord > pow(10,-10)) {
         std::cerr << "error with volume sum" << std::endl;
         std::cerr << "sum LV-V_total" << sum_LV-V_total << std::endl;
         std::cerr << "stop" << std::endl;
         exit(0);
     }
}

void voltri(double& V_total, double* V_atom, double* S_atom, double* S_lag_int,
            double** coords, double* weights, tri_type* trilist, double* S_lag2)
{
    int i1, i2, i3, tri[3];
    double w1, w2, w3, trifrac, v1_2, v2_2, v3_2, v4_2, v6_2, v1dv2, v1dv4, v1dv6, v2dv4, v2dv6, v4dv6;
    double normN1, normN2, normN3, normN4, r_normN1, r_normN2, r_normN3, r_normN4;
    double r_denom1, r_denom2, r_denom3, r_denom4, r_denom5, r_denom6, num1, num2, num3, num4, num5, num6;
    double n1dn2, n1dn3, n1dn4, n2dn3, n2dn4, n3dn4, LIS_res12, LIS_res13, LIS_res23;
    double h_third, Two_P_1, Two_P_2, Two_P_3, a_2, a12, l2, m2, a, a_1, l, m, c, c2, b, b2, b1, b12;
    double r1, r2, r3, l_v1, l_v2, l_v3, h12_1, h12_2, h13_1, h13_3, h23_2, h23_3, k1, k2, k3, temp;
    double omega1, omega2, omega3, a1, a2, a3, a4, a5, a6, r2_circle_12, r2_circle_13, r2_circle_23;
    double two_tri12, two_tri13, two_tri23, r, s, t, h, vol_tri, V12, V13, V23, temp_2;
    double S1_sphr, S2_sphr, S3_sphr, S1_lag, S2_lag, S3_lag, S12_lag, S13_lag, S23_lag;
    double x_char[4], p[3], p1[3], p2[3], p3[3], v1[3], v2[3], v3[3], v4[3], v5[3], v6[3], n1[3];
    double x12[3], x13[3], x23[3], d1[3], d2[3], d1xd2[3], u1[3], u2[3], u3[3], u4[3], temp_vec[3];
    double normal1[3], normal2[3];
    tri_type* current = trilist->next_ext;
    while(current) {
        for (int i = 0; i < 4; ++ i) x_char[i] = current->x_char[i];
        for (int i = 0; i < 3; ++ i) tri[i] = current->tri[i];
        trifrac = current->frac;
        i1 = tri[0];
        i2 = tri[1];
        i3 = tri[2];
        for (int i = 0; i < 3; ++ i) {
            p1[i] = coords[i1-1][i];
            p2[i] = coords[i2-1][i];
            p3[i] = coords[i3-1][i];
        }
        w1 = weights[i1-1];
        w2 = weights[i2-1];
        w3 = weights[i3-1];
        v1_2 = 0;
        v2_2 = 0;
        v1dv2 = 0;

        for (int i = 0; i < 3; ++ i) {
            v1[i] = p1[i] - p2[i];
            v2[i] = p3[i] - p2[i];
            v1_2 += v1[i]*v1[i];
            v2_2 += v2[i]*v2[i];
            v1dv2 += v1[i]*v2[i];
        }
        n1[0] = v1[1]*v2[2] - v1[2]*v2[1];
        n1[1] = -v1[0]*v2[2] + v1[2]*v2[0];
        n1[2] = v1[0]*v2[1] - v1[1]*v2[0];
        normN1 = sqrt(v1_2*v2_2 - v1dv2*v1dv2);
        r_normN1 = 1./normN1;
        h = sqrt(-x_char[3]);
        v4_2 = 0;
        v6_2 = 0;
        v1dv4 = 0;
        v1dv6 = 0;
        v2dv4 = 0;
        v2dv6 = 0;
        v4dv6 = 0;
        for (int i = 0;i < 3; ++ i) {
            n1[i] *= r_normN1;
            p[i] = x_char[i] - h*n1[i];
            v4[i] = p[i] - p1[i];
            v6[i] = p[i] - p3[i];
            v4_2 += v4[i]*v4[i];
            v6_2 += v6[i]*v6[i];
            v1dv4 += v1[i] * v4[i];
            v1dv6 += v1[i] * v6[i];
            v2dv4 += v2[i] * v4[i];
            v2dv6 += v2[i] * v6[i];
            v4dv6 += v4[i] * v6[i];
        }

        normN2 = sqrt(v1_2*v4_2 - v1dv4*v1dv4);
        normN3 = sqrt(v2_2*v6_2 - v2dv6*v2dv6);
        normN4 = sqrt(v4_2*v6_2 - v4dv6*v4dv6);
        r_normN2 = 1.0/normN2;  //=1/|N2| etc.
        r_normN3 = 1.0/normN3;
        r_normN4 = 1.0/normN4;

        r_denom1 = r_normN1*r_normN2; //reciprocal of denominator for cos(a1), etc.
        r_denom2 = r_normN1*r_normN3;
        r_denom3 = r_normN1*r_normN4;
        r_denom4 = r_normN2*r_normN4;
        r_denom5 = r_normN2*r_normN3;
        r_denom6 = r_normN3*r_normN4;

        num1 = v1_2*v2dv4 - v1dv2*v1dv4; //numerator for cos(a1) = N1 * N2
        num2 = v2_2*v1dv6 - v1dv2*v2dv6; //numerator for cos(a2) = N1 * N3
        num3 = v1dv4*v2dv6 - v1dv6*v2dv4; //numerator for cos(a3)= N1 * N4
        num4 = v4_2*v1dv6 - v4dv6*v1dv4; //numerator for cos(a4) = -N2 * N4
        num5 = v4dv6*v1dv2 - v2dv4*v1dv6; //numerator for cos(a5) = -N2 * N3
        num6 = v6_2*v2dv4 - v4dv6*v2dv6; //numerator for cos(a6) = -N3 * N4

        n1dn2 = num1*r_denom1; //=dot_product(n1,n2)
        if (n1dn2 >= 1.0) a1 = 0.0;
        else if (n1dn2 <= -1.0) a1 = Pi;
        else a1 = acos(n1dn2);

        n1dn3 = num2*r_denom2;
        if (n1dn3 >= 1.0) a2 = 0.0;
        else if (n1dn3 <= -1.0) a2 = Pi;
        else a2 = acos(n1dn3);

        n1dn4 = num3*r_denom3;
        if (n1dn4 >= 1.0) a3 = 0.0;
        else if (n1dn4 <= -1.0) a3 = Pi;
        else a3 = acos(n1dn4);

        n2dn4 = num4*r_denom4;
        if (n2dn4 >= 1.0) a4 = 0.0;
        else if (n2dn4 <= -1.0) a4 = Pi;
        else a4 = acos(n2dn4);

        n2dn3 = num5*r_denom5;
        if (n2dn3 >= 1.0) a5 = 0.0;
        else if (n2dn3 <= -1.0) a5 = Pi;
        else a5 = acos(n2dn3);

        n3dn4 = num6*r_denom6;
        if (n3dn4 >= 1.0) a6 = 0.0;
        else if (n3dn4 <= -1.0) a6 = Pi;
        else a6 = acos(n3dn4);

        omega1 = a1+a3+a4-Pi;
        omega2 = a1+a2+a5-Pi;
        omega3 = a2+a3+a6-Pi;
        v3_2 = 0;
        for (int i = 0; i < 3; ++ i) {
            v5[i] = p[i] - p2[i];
            v3[i] = p3[i] - p1[i];
            v3_2 += v3[i]*v3[i];
        }
        r1=sqrt(w1); r2=sqrt(w2); r3=sqrt(w3);
        l_v1=sqrt(v1_2); l_v2=sqrt(v2_2); l_v3=sqrt(v3_2);
        k1=w1-w2; k2=w1-w3; k3=w2-w3;
        temp=0.5*(k1+v1_2)/l_v1;
        h12_1=r1-temp;
        h12_2=temp-l_v1+r2;
        r2_circle_12=w1-pow(temp,2);
        temp=0.5*(k2+v3_2)/l_v3;
        h13_1=r1-temp;
        h13_3=temp-l_v3+r3;
        r2_circle_13=w1-pow(temp,2);
        temp=0.5*(k3+v2_2)/l_v2;
        h23_2=r2-temp;
        h23_3=temp-l_v2+r3;
        r2_circle_23=w2-pow(temp,2);
        // pyramid computations
        //I CAN MAKE THIS PART FASTER
        r = 0.5*(k1/v1_2 + 1);
        for (int i = 0; i < 3; ++ i) x12[i] = p1[i] + r*(p2[i] - p1[i]);
        s = 0.5*(k2/v3_2 + 1);
        for (int i = 0; i < 3; ++ i) x13[i] = p1[i] + s*(p3[i] - p1[i]);
        t = 0.5*(k3/v2_2 + 1);
        for (int i = 0; i < 3; ++ i) x23[i] = p2[i] + t*(p3[i] - p2[i]);
        h_third = h*One_third;

        for (int i = 0; i < 3; ++ i) {
            u1[i] = p1[i] - x13[i];
           	u2[i] = x12[i] - p1[i];
            u3[i] = x_char[i] - x12[i];
            u4[i] = x13[i] - x_char[i];
        }

        operation::cross_product(u4, u1, temp_vec);
        double dot_product = 0;
        temp_2 = 0;
        for (int i = 0; i < 3; ++ i) {
            dot_product += temp_vec[i] * (-n1[i]);
            temp_2 += temp_vec[i] * temp_vec[i];
        }
        if (dot_product >= 0) Two_P_1 = sqrt(temp_2);
        else Two_P_1 = -sqrt(temp_2);

        operation::cross_product(u2, u3, temp_vec);
        dot_product = 0;
        temp_2 = 0;
        for (int i = 0; i < 3; ++ i) {
            dot_product += temp_vec[i] * (-n1[i]);
            temp_2 += temp_vec[i] * temp_vec[i];
        }
        if (dot_product >= 0) Two_P_1 += sqrt(temp_2);
        else Two_P_1 -= sqrt(temp_2);
        Two_P_1 = Two_P_1 * h_third;

        for (int i = 0; i < 3; ++ i) {
            u1[i] = p2[i] - x12[i];
            u2[i] = x23[i] - p2[i];
            u3[i] = x_char[i] - x23[i];
            u4[i] = x12[i] - x_char[i];
        }

        operation::cross_product(u4, u1, temp_vec);
        dot_product = 0;
        temp_2 = 0;
        for (int i = 0; i < 3; ++ i) {
            dot_product += temp_vec[i] * (-n1[i]);
            temp_2 += temp_vec[i] * temp_vec[i];
        }
        if (dot_product >= 0) Two_P_2 = sqrt(temp_2);
        else Two_P_2 = -sqrt(temp_2);
        operation::cross_product(u2, u3, temp_vec);
        dot_product = 0;
        temp_2 = 0;
        for (int i = 0; i < 3; ++ i) {
            dot_product += temp_vec[i] * (-n1[i]);
            temp_2 += temp_vec[i] * temp_vec[i];
        }
        if (dot_product >= 0) Two_P_2 += sqrt(temp_2);
        else Two_P_2 -= sqrt(temp_2);
        Two_P_2 = Two_P_2 * h_third;

        for (int i = 0; i < 3; ++ i) {
            u1[i] = p3[i] - x23[i];
            u2[i] = x13[i] - p3[i];
            u3[i] = x_char[i] - x13[i];
            u4[i] = x23[i] - x_char[i];
        }

        operation::cross_product(u4, u1, temp_vec);
        dot_product = 0;
        temp_2 = 0;
        for (int i = 0; i < 3; ++ i) {
            dot_product += temp_vec[i] * (-n1[i]);
            temp_2 += temp_vec[i] * temp_vec[i];
        }
        if (dot_product >= 0) Two_P_3 = sqrt(temp_2);
        else Two_P_3 = -sqrt(temp_2);
        operation::cross_product(u2, u3, temp_vec);
        dot_product = 0;
        temp_2 = 0;
        for (int i = 0; i < 3; ++ i) {
            dot_product += temp_vec[i] * (-n1[i]);
            temp_2 += temp_vec[i] * temp_vec[i];
        }
        if (dot_product >= 0) Two_P_3 += sqrt(temp_2);
        else Two_P_3 -= sqrt(temp_2);
        Two_P_3 = Two_P_3 * h_third;

	    double cross_product[3];
	    operation::cross_product(v1, v2, cross_product);
	    dot_product = 0;
	    for (int i = 0; i < 3; ++ i) dot_product += v5[i] * cross_product[i];
	    if (dot_product < 0) dot_product = -dot_product;
	    vol_tri = One_third*dot_product;
        vol_tri = vol_tri - Two_thirds*(omega1*pow(r1,3) + omega2*pow(r2,3) + omega3*pow(r3,3));
	    V12 =.083333333333333333/l_v1*pow((r1+r2-l_v1),2)*(v1_2+2*l_v1*r2-3*w2+2*l_v1*r1+6*r1*r2-3*w1);   //actually V12/pi
	    V13 =.083333333333333333/l_v3*pow((r1+r3-l_v3),2)*(v3_2+2*l_v3*r3-3*w3+2*l_v3*r1+6*r1*r3-3*w1);   //actually V12/pi
	    V23 =.083333333333333333/l_v2*pow((r2+r3-l_v2),2)*(v2_2+2*l_v2*r3-3*w3+2*l_v2*r2+6*r2*r3-3*w2);   //actually V12/pi
	    vol_tri=vol_tri+(V12*a1+V13*a3+V23*a2);  //contribution from edges

	    V_total=V_total+trifrac*vol_tri;
        V_atom[i1-1] = V_atom[i1-1]+trifrac*(a1*pow(h12_1,2)*(r1-h12_1*One_third)+a3*pow(h13_1,2)*(r1-h13_1*One_third)
                                       -omega1*w1*r1*Two_thirds);
	    V_atom[i1-1] = V_atom[i1-1] + trifrac*Two_P_1;

        V_atom[i2-1] = V_atom[i2-1] + trifrac*(a1*pow(h12_2,2)*(r2-h12_2*One_third)+a2*pow(h23_2,2)*(r2-h23_2*One_third)
                                       -omega2*w2*r2*Two_thirds);

        V_atom[i2-1] = V_atom[i2-1] + trifrac*Two_P_2;

        V_atom[i3-1] = V_atom[i3-1] + trifrac*(a3*pow(h13_3,2)*(r3-h13_3*One_third)+a2*pow(h23_3,2)*(r3-h23_3*One_third)
                                       -omega3*w3*r3*Two_thirds);

        V_atom[i3-1] = V_atom[i3-1]+trifrac*Two_P_3;
        double base1, base2, base3, temp1[3], temp2[3], temp3[3];
	    base1 = 0;
	    base2 = 0;
	    base3 = 0;
	    for (int i = 0; i < 3; ++ i) {
            temp1[i] = x12[i] - x_char[i];
            temp2[i] = x13[i] - x_char[i];
            temp3[i] = x23[i] - x_char[i];
	        base1 += temp1[i] * temp1[i];
	        base2 += temp2[i] * temp2[i];
	        base3 += temp3[i] * temp3[i];
	    }
	    base1 = sqrt(base1);
	    base2 = sqrt(base2);
	    base3 = sqrt(base3);
	    two_tri12 = base1 * h;
	    two_tri13 = base2 * h;
	    two_tri23 = base3 * h;
        S1_sphr=trifrac*(2.0*(a1*r1*h12_1+a3*r1*h13_1-omega1*w1));
        S2_sphr=trifrac*(2.0*(a1*r2*h12_2+a2*r2*h23_2-omega2*w2));
        S3_sphr=trifrac*(2.0*(a3*r3*h13_3+a2*r3*h23_3-omega3*w3));
	    if (a1 <= Pi_half) S12_lag=trifrac*(a1*r2_circle_12-two_tri12);
	    else S12_lag=trifrac*(a1*r2_circle_12+two_tri12);

	    if (a3 <= Pi_half) S13_lag=trifrac*(a3*r2_circle_13-two_tri13);
	    else S13_lag=trifrac*(a3*r2_circle_13+two_tri13);

	    if (a2 <= Pi_half) S23_lag=trifrac*(a2*r2_circle_23-two_tri23);
	    else S23_lag=trifrac*(a2*r2_circle_23+two_tri23);

	    S1_lag=S12_lag+S13_lag;
	    S2_lag=S12_lag+S23_lag;
	    S3_lag=S13_lag+S23_lag;

	    S_lag_int[i1-1]=S_lag_int[i1-1]+S1_sphr-S1_lag;
	    S_lag_int[i2-1]=S_lag_int[i2-1]+S2_sphr-S2_lag;
	    S_lag_int[i3-1]=S_lag_int[i3-1]+S3_sphr-S3_lag;

	    S_atom[i1-1]=S_atom[i1-1]+S1_sphr;
	    S_atom[i2-1]=S_atom[i2-1]+S2_sphr;
	    S_atom[i3-1]=S_atom[i3-1]+S3_sphr;

	    current = current->next_ext;

	    double r12=sqrt(r2_circle_12);
	    double r13=sqrt(r2_circle_13);
            double r23=sqrt(r2_circle_23);

	    temp_2 = 0;
	    for (int i = 0; i < 3; ++ i) {
	    	    temp1[i] = x_char[i] - x12[i];
                temp_2 += temp1[i]*temp1[i];
        }
	    double d = sqrt(temp_2);
	    double phi = acos(d/r12);

	    for (int i = 0; i < 3; ++ i) {
	    	    temp2[i] = coords[i1-1][i] - coords[i2-1][i];
	    	    temp3[i] = coords[i3-1][i] - coords[i2-1][i];
	    	    temp1[i] = x_char[i] - coords[i2-1][i];
	    }
	    operation::cross_product(temp2, temp3, normal1);
	    operation::cross_product(temp2, temp1, normal2);
	    double test = 0;
	    for (int i = 0; i < 3; ++ i) test += normal1[i]*normal2[i];

        if (test >= 0) LIS_res12 = (Pi-phi)*pow(r12,2) + d*sqrt(pow(r12,2)-pow(d,2));
        else  LIS_res12 = phi*pow(r12,2) - d*sqrt(pow(r12,2)-pow(d,2));

        temp_2 = 0;
        for (int i = 0; i < 3; ++ i) {
                temp1[i] = x_char[i] - x13[i];
                temp_2 += temp1[i]*temp1[i];
        }
        d = sqrt(temp_2);
        phi = acos(d/r13);

        for (int i = 0; i < 3; ++ i) {
                temp2[i] = coords[i1-1][i] - coords[i3-1][i];
                temp3[i] = coords[i2-1][i] - coords[i3-1][i];
                temp1[i] = x_char[i] - coords[i3-1][i];
        }
        operation::cross_product(temp2, temp3, normal1);
        operation::cross_product(temp2, temp1, normal2);
        test = 0;
        for (int i = 0; i < 3; ++ i) test += normal1[i]*normal2[i];

        if (test >= 0) LIS_res13 = (Pi-phi)*pow(r13,2) + d*sqrt(pow(r13,2)-pow(d,2));
        else  LIS_res13 = phi*pow(r13,2) - d*sqrt(pow(r13,2)-pow(d,2));

        temp_2 = 0;
        for (int i = 0; i < 3; ++ i) {
                temp1[i] = x_char[i] - x23[i];
                temp_2 += temp1[i]*temp1[i];
        }
        d = sqrt(temp_2);
        phi = acos(d/r23);

        for (int i = 0; i < 3; ++ i) {
                temp2[i] = coords[i2-1][i] - coords[i3-1][i];
                temp3[i] = coords[i1-1][i] - coords[i3-1][i];
                temp1[i] = x_char[i] - coords[i3-1][i];
        }
        operation::cross_product(temp2, temp3, normal1);
        operation::cross_product(temp2, temp1, normal2);
        test = 0;
        for (int i = 0; i < 3; ++ i) test += normal1[i]*normal2[i];

        if (test >= 0) LIS_res23 = (Pi-phi)*pow(r23,2) + d*sqrt(pow(r23,2)-pow(d,2));
        else  LIS_res23 = phi*pow(r23,2) - d*sqrt(pow(r23,2)-pow(d,2));

        S_lag2[i1-1]=S_lag2[i1-1]+LIS_res12+LIS_res13;
        S_lag2[i2-1]=S_lag2[i2-1]+LIS_res12+LIS_res23;
        S_lag2[i3-1]=S_lag2[i3-1]+LIS_res23+LIS_res13;
    }
}

void voledg(double** coords, double* weights, edge_type* edgelist, int* idT, double& V_total, double* V_atom, double* S_atom, double* S_lag_int, int** T, std::vector<std::vector<double> > &surfint, int& idx)
{
    int e1, e2, e3, e4, index, num_ring, curr_T, edg[2], tetra_ring[20];
    double w1, w2, r1, r2, temp, l_v1, r_l_v1, v1_2, v2_2, v4_2, v1dv2, v1dv4, v2dv4, h1, h2, normN1, normN2, dihedral, num1,
            r_denom1, a1, n1dn2, half_r2_circle, S1_sphr, S2_sphr, S_lag;
    double p1[3], p2[3], p3[3], p4[3], v1[3], v2[3], v4[3], v1_unit[3];
    edge_type* current = edgelist->next_ext;

    while (current) {
        edg[0] = current->edg[0];
        edg[1] = current->edg[1];
        e1 = edg[0];
        e2 = edg[1];

        for (int i = 0; i < 3; ++ i) {
            p1[i] = coords[e1-1][i];
            p2[i] = coords[e2-1][i];
        }
        w1 = weights[e1-1];
        w2 = weights[e2-1];
        v1_2 = 0;
        for (int i = 0; i < 3; ++ i) {
            v1[i] = p1[i] - p2[i];
            v1_2 += v1[i]*v1[i];
        }
        l_v1 = sqrt(v1_2);
        r1 = pow(w1, 0.5);
        r2 = pow(w2, 0.5);
        temp = 0.5*(w1-w2+v1_2)/l_v1;
        h1 = r1 - temp;
        h2 = temp - l_v1 + r2;
        current->h1 = temp;
        current->h2 = l_v1 - temp;
        num_ring = current->num_ring;

        for (int i = 0; i < num_ring; ++ i) tetra_ring[i] = current->tetra_ring[i];
        dihedral = Pi2;

        for (int j = 0; j < num_ring; ++ j) {
            curr_T = tetra_ring[j];
            if (curr_T != 0) {
                if (idT[curr_T+1] < 0) {
                    index = 1;
                    for (int i = 0; i < 4; ++ i) {
                        if (T[curr_T-1][i] != e1 && T[curr_T-1][i] != e2) {
                            if (index == 1) {
                                e3 = T[curr_T-1][i];
                                index += 1;
                            } else {
                                e4 = T[curr_T-1][i];
                                break;
                            }
                        }
                    }
                    for (int k = 0; k < 3; ++ k) {
                        p3[k] = coords[e3-1][k];
                        p4[k] = coords[e4-1][k];
                        v2[k] = p3[k] - p2[k];
                        v4[k] = p4[k] - p2[k];
                    }

                    v2_2 = 0;
                    v4_2 = 0;
                    v1dv2 = 0;
                    v1dv4 = 0;
                    v2dv4 = 0;

                    for (int k = 0; k < 3; ++ k) {
                        v2_2 += v2[k]*v2[k];
                        v4_2 += v4[k]*v4[k];
                        v1dv2 += v1[k]*v2[k];
                        v1dv4 += v1[k]*v4[k];
                        v2dv4 += v2[k]*v4[k];
                    }

                    normN1 = sqrt(v1_2*v2_2 - v1dv2*v1dv2);
                    normN2 = sqrt(v1_2*v4_2 - v1dv4*v1dv4);

                    r_denom1 = 1.0/(normN1*normN2);
                    num1 = v1_2*v2dv4 - v1dv2*v1dv4;
                    n1dn2 = num1*r_denom1;

                    if (n1dn2 >= 1.0) a1 = 0;
                    else if (n1dn2 <= -1.0) a1 = Pi;
                    else a1 = acos(n1dn2);

                    dihedral = dihedral - a1;
                }
            }
        }
        V_total = V_total - 0.5 * dihedral * (pow(h1,2)*(r1-h1*One_third) + pow(h2,2)*(r2-h2*One_third) );

        V_atom[e1-1] = V_atom[e1-1] - 0.5*dihedral*pow(h1,2)*(r1 - h1*One_third);
        V_atom[e2-1] = V_atom[e2-1] - 0.5*dihedral*pow(h2,2)*(r2 - h2*One_third);
        half_r2_circle = 0.5*(w1 - pow(temp,2));
        S_lag = dihedral * half_r2_circle;
        S1_sphr = dihedral * r1 * h1;
        S2_sphr = dihedral * r2 * h2;
        S_lag_int[e1-1] = S_lag_int[e1-1] - S1_sphr + S_lag;
        S_lag_int[e2-1] = S_lag_int[e2-1] - S2_sphr + S_lag;
        S_atom[e1-1] -= S1_sphr;
        S_atom[e2-1] -= S2_sphr;

        std::vector<double> csurf;
        if (e1 > e2) {
            csurf.push_back(e2);
            csurf.push_back(e1);
        } else {
            csurf.push_back(e1);
            csurf.push_back(e2);
        }
        csurf.push_back(S_lag);
        surfint.push_back(csurf);
        idx = idx+1;
        csurf.clear();

        current = current->next_ext;
     }
}

void volver(int& num_ver_ext, int& maxupver, int* ver_ext, ver_type* verlist, double* weights, double** coords, double& V_total, double* V_atom, double* S_atom, double* S_lag_int, int* idT, int** T)
{
    int numupver, i1, i2, i3, i4, tetra_index;
    int* tetra_list = new int[maxupver];
    const double Pi4 = 12.566370614359172;
    double solid_angle, aa, ab, ac, adotb, adotc, bdotc, denom, triple_product, F, s_angle, w1, S_sphr;
    double axb[3], axc[3], bxc[3];

    for (int i = 0; i < num_ver_ext; ++ i) {
        i1 = ver_ext[i];
        solid_angle = Pi4;
        numupver = verlist[i1-1].last_upT;
        for (int j = 0; j < numupver; ++ j) tetra_list[j] = verlist[i1-1].upverT[j];
        double p1[3] = {coords[i1-1][0], coords[i1-1][1], coords[i1-1][2]};
        w1 = weights[i1-1];
        for (int j = 0; j < numupver; ++j) { // 3*numver dim of upver right now
        //loop through tetra that are incident to point p1
            tetra_index = tetra_list[j];
            if (idT[tetra_index+1] < 0) { //in C(0)
                if (i1 == T[tetra_index-1][0]) {
                    i2 = T[tetra_index-1][1];
                    i3 = T[tetra_index-1][2];
                    i4 = T[tetra_index-1][3];
                }
                else if (i1 == T[tetra_index-1][1]) {
                    i2 = T[tetra_index-1][0];
                    i3 = T[tetra_index-1][3];
                    i4 = T[tetra_index-1][2];
                }
                else if (i1 == T[tetra_index-1][2]) {
                    i2=T[tetra_index-1][3];
                    i3=T[tetra_index-1][0];
                    i4=T[tetra_index-1][1];
                }
                else {
                    i2=T[tetra_index-1][0];
                    i3=T[tetra_index-1][2];
                    i4=T[tetra_index-1][1];
                }
                double a[3] = {coords[i2-1][0]-p1[0], coords[i2-1][1]-p1[1], coords[i2-1][2]-p1[2]};
                double b[3] = {coords[i3-1][0]-p1[0], coords[i3-1][1]-p1[1], coords[i3-1][2]-p1[2]};
                double c[3] = {coords[i4-1][0]-p1[0], coords[i4-1][1]-p1[1], coords[i4-1][2]-p1[2]};
                aa = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
                ab = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
                ac = sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
                adotb = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
                adotc = a[0]*c[0]+a[1]*c[1]+a[2]*c[2];
                bdotc = b[0]*c[0]+b[1]*c[1]+b[2]*c[2];
                denom = aa*ab*ac+adotb*ac+adotc*ab+bdotc*aa;
                axb[0] = a[1]*b[2]-a[2]*b[1];
                axb[1] = -a[0]*b[2]+a[2]*b[0];
                axb[2] = a[0]*b[1]-a[1]*b[0];

                axc[0] = a[1]*c[2]-a[2]*c[1];
                axc[1] = -a[0]*c[2]+a[2]*c[0];
                axc[2] = a[0]*c[1]-a[1]*c[0];

                bxc[0] =  b[1]*c[2]-b[2]*c[1];
                bxc[1] = -b[0]*c[2]+b[2]*c[0];
                bxc[2] =  b[0]*c[1]-b[1]*c[0];


                triple_product = a[0]*bxc[0] + a[1]*bxc[1] + a[2]*bxc[2];
                //orientation is such that triple_product>+0
                F = triple_product/denom;
                s_angle = 2*atan(F);
                if (denom < 0) s_angle=s_angle + Pi2;
                solid_angle = solid_angle - s_angle;
            }
        }
        V_total = V_total + solid_angle * One_third * pow(w1,1.5);
        V_atom[i1-1] = V_atom[i1-1] + solid_angle * One_third * pow(w1,1.5);
        S_sphr = solid_angle * w1;
        S_lag_int[i1-1] = S_lag_int[i1-1] + S_sphr;
        S_atom[i1-1] = S_atom[i1-1] + S_sphr;
    }
    delete[] tetra_list;
}

void classify_ver(int& num_coord, int& num_ver_ext, double* weights, ver_type* verlist, int* ver_ext, edge_type* refedg, double** coords)
{
    int ext, testpt, num_int, num_out;
    double wp, d2, att, tp[3], pi[3];
    num_ver_ext = 0;
    for (int i = 0; i < num_coord; ++ i) {
        if (weights[i] > 0 && verlist[i].num_adj_edg > 0) {
            num_int=0;
            num_out=0;
            ext=-1 ; //default=interior
            for (int j = 0; j < verlist[i].num_adj_edg; ++ j) {
                if (refedg[verlist[i].adj_edg[j]-1].alf == -1) { //exterior
                    ext = 1;
                    break;
                }
                else if (refedg[verlist[i].adj_edg[j]-1].alf == 1) num_out = 1;
                else num_int = 1;
            }

            if (ext == 1) { //ver is exterior
                num_ver_ext += 1;
                ver_ext[num_ver_ext-1] = i+1;
                verlist[i].alf = -1;
            }
            else if (num_int == 1) { //edge is interior
                verlist[i].alf = -2;
            }
            else if (num_out == 1) { //test if attacheda
                ext = 1; //not in C(0) new default
                pi[0] = coords[i][0];
                pi[1] = coords[i][1];
                pi[2] = coords[i][2];
                for (int j = 0; j < verlist[i].num_adj_edg; ++ j) {
                    testpt = verlist[i].op_ver[j];
                    tp[0] = coords[testpt-1][0];
                    tp[1] = coords[testpt-1][1];
                    tp[2] = coords[testpt-1][2];
                    wp = weights[testpt-1];
                    wp = weights[testpt-1];
                    wp = weights[testpt-1];
                    d2 = (pi[0]-tp[0])*(pi[0]-tp[0])+(pi[1]-tp[1])*(pi[1]-tp[1])+(pi[2]-tp[2])*(pi[2]-tp[2]);
                    att = d2 + weights[i] - wp;
                    if (att < 0) { //attached
                        ext = 0;
                        break;
                    }
                }
                if (ext == 0)
                    verlist[i].alf = 1;
                else {
                    num_ver_ext += 1;
                    ver_ext[num_ver_ext-1] = i + 1;
                    verlist[i].alf = -1;
                }
            }
            else {
                std::cerr << "error in classify_ver; stop" << std::endl;
                exit(0);
            }
        }
        else {//weight(i)<=0 or no adjacent edges ==> not in C(0)
                verlist[i].alf = 1;
        }// weight(i)>0
    }
}

void classify_edg(int& num_edg_ext, int& num_int_edge, int& num_edg, edge_type* edgelist, edge_type* refedg, double** coords, double* lift, double* weights, tri_type* reftri, int* int_edg_list)
{
    edge_type* current_ext;
    double t, dist, absp12, x[4];
    int ext, ver_test, num_int, num_out, e1, e2, num_adj_tri;
    int adj_tri[20];
    num_edg_ext = 0;
    num_int_edge = 0;
    current_ext = edgelist;
    for (int i = 0; i < num_edg; ++ i) {
        refedg[i].disc = 0;
        e1 = refedg[i].edg[0];
        e2 = refedg[i].edg[1];
        double p1[4] = {coords[e1-1][0], coords[e1-1][1], coords[e1-1][2], lift[e1-1]};
        double p2[4] = {coords[e2-1][0], coords[e2-1][1], coords[e2-1][2], lift[e2-1]};
        double p12[4] = {p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2], p1[3]-p2[3]};
        absp12 = sqrt(p12[0]*p12[0] + p12[1]*p12[1] + p12[2]*p12[2]);
        t = p12[3] - 2.0*(p2[0]*p12[0] + p2[1]*p12[1] + p2[2]*p12[2]);
        t = t / 2.0 / (p12[0]*p12[0] + p12[1]*p12[1] + p12[2]*p12[2]);
        x[0] = p12[0]*t + p2[0];
        x[1] = p12[1]*t + p2[1];
        x[2] = p12[2]*t + p2[2];
        x[3] = (p1[0]-x[0])*(p1[0]-x[0]) + (p1[1]-x[1])*(p1[1]-x[1]) + (p1[2]-x[2])*(p1[2]-x[2]) - weights[refedg[i].edg[0]-1];
        refedg[i].x_char[0] = x[0];
        refedg[i].x_char[1] = x[1];
        refedg[i].x_char[2] = x[2];
        ext = 5;  //arbitrary value, later will be set to -2, -1, 0 (int, ext, out)
        if (x[3] >= 0) {
            refedg[i].alf = 1;
        } else {
            num_int = 0;
            num_out = 0;
            num_adj_tri = refedg[i].num_adj_tri;
            for (int j = 0; j < num_adj_tri; ++ j) {
                adj_tri[j] = refedg[i].adj_tri[j];
            }
            for (int j = 0; j < num_adj_tri; ++ j) {
                if (reftri[adj_tri[j]-1].alf == -1) {
                    ext = 1;
                    break;
                }
                else if (reftri[refedg[i].adj_tri[j]-1].alf == 1) {
                    num_out = 1;
                }
                else {
                    num_int = 1;
                }
            }
            if (ext == 1) {
                refedg[i].alf = -1;
                refedg[i].x_char[0] = x[0];
                refedg[i].x_char[1] = x[1];
                refedg[i].x_char[2] = x[2];
                refedg[i].t = t;
                num_edg_ext += 1;
                current_ext->next_ext = &refedg[i];
                current_ext = &refedg[i];
            }
            else if (num_int == 1) {
                refedg[i].alf = -2;
                num_int_edge += 1;
                int_edg_list[num_int_edge-1] = i+1;
            }
            else {
                ext = 1;
                for (int j = 0; j < num_adj_tri; ++ j) {
                    ver_test = refedg[i].op_ver[j];
                    double difference[3] = {x[0]-coords[ver_test-1][0], x[1]-coords[ver_test-1][1], x[2]-coords[ver_test-1][2]};
                    dist = difference[0]*difference[0] + difference[1]*difference[1] + difference[2]*difference[2] - x[3] - weights[ver_test-1];
                    if (dist <= 0) {
                        ext = 0;
                        break;
                    }
                }
                if (ext == 0) {
                    refedg[i].alf = 1;
                } else {
                    refedg[i].alf = -1;
                    refedg[i].x_char[0] = x[0];
                    refedg[i].x_char[1] = x[1];
                    refedg[i].x_char[2] = x[2];
                    refedg[i].t = t;
                    num_edg_ext += 1;
                    current_ext->next_ext = &refedg[i];
                    current_ext = &refedg[i];
                }
            }
        }
    }
    current_ext->next_ext = NULL;
}

void classify_tri(int* alfT, int& num_tri_ext, int& num_edg, int& num_tri, tri_type* reftri, tri_type* trilist, double** coords, double* lift, double *weights)
{
    int classT1, classT2, tri_c[3], attached;
    double det, id, dist;
    double x[4], A[3], B[3], C[3], f[3];
    alfT[0] = 1;
    num_tri_ext = 0;
    tri_type* current_ext = trilist;
    for (int i = 0; i < num_tri; ++ i) {
        classT1 = alfT[reftri[i].adjT[0] ];
        classT2 = alfT[reftri[i].adjT[1] ];
        if (classT1 < 0) {
            if (classT2 >= 0) {
                num_tri_ext += 1;
                current_ext->next_ext = &reftri[i];
                current_ext = &reftri[i];
                reftri[i].frac = 0.5;
                reftri[i].alf = -1;
                for (int j = 0; j < 3; ++ j) tri_c[j] = reftri[i].tri[j];
                for (int j = 0; j < 3; ++ j) {
                    A[j] = coords[tri_c[0]-1][j] - coords[tri_c[2]-1][j];
                    B[j] = coords[tri_c[1]-1][j] - coords[tri_c[2]-1][j];
                }
                C[0] = A[1]*B[2] - A[2]*B[1];
                C[1] = A[2]*B[0] - A[0]*B[2];
                C[2] = A[0]*B[1] - A[1]*B[0];
                f[0] = 0.5*(lift[tri_c[0]-1] - lift[tri_c[2]-1]);
                f[1] = 0.5*(lift[tri_c[1]-1] - lift[tri_c[2]-1]);
                f[2] = C[0]*coords[tri_c[2]-1][0] + C[1]*coords[tri_c[2]-1][1] + C[2]*coords[tri_c[2]-1][2];
                det = C[0]*C[0] + C[1]*C[1] + C[2]*C[2];
                id = 1/det;
                x[0] = id*(f[0]*(B[1]*C[2] - B[2]*C[1]) + f[1]*(A[2]*C[1] - A[1]*C[2]) + f[2]*C[0]);
                x[1] = id*(f[0]*(B[2]*C[0] - B[0]*C[2]) + f[1]*(A[0]*C[2] - A[2]*C[0]) + f[2]*C[1]);
                x[2] = id*(f[0]*(B[0]*C[1] - B[1]*C[0]) + f[1]*(A[1]*C[0] - A[0]*C[1]) + f[2]*C[2]);
                double temp[3] = {x[0] - coords[tri_c[2]-1][0], x[1] - coords[tri_c[2]-1][1], x[2] - coords[tri_c[2]-1][2]};
                x[3] = temp[0]*temp[0] + temp[1]*temp[1] + temp[2]*temp[2] - weights[tri_c[2]-1];
                reftri[i].x_char[0] = x[0];
                reftri[i].x_char[1] = x[1];
                reftri[i].x_char[2] = x[2];
                reftri[i].x_char[3] = x[3];
            } else {
                reftri[i].alf = -2;
            }
        } else {
            if (classT2 < 0) {
                num_tri_ext += 1;
                current_ext->next_ext = &(reftri[i]);
                current_ext = &reftri[i];
                reftri[i].frac = 0.5;
                reftri[i].alf = -1;
                for (int j = 0; j < 3; ++ j) tri_c[j] = reftri[i].tri[j];
                for (int j = 0; j < 3; ++ j) {
                    A[j] = coords[tri_c[0]-1][j] - coords[tri_c[2]-1][j];
                    B[j] = coords[tri_c[1]-1][j] - coords[tri_c[2]-1][j];
                }
                C[0] = A[1]*B[2] - A[2]*B[1];
                C[1] = A[2]*B[0] - A[0]*B[2];
                C[2] = A[0]*B[1] - A[1]*B[0];
                f[0] = 0.5*(lift[tri_c[0]-1] - lift[tri_c[2]-1]);
                f[1] = 0.5*(lift[tri_c[1]-1] - lift[tri_c[2]-1]);
                f[2] = C[0]*coords[tri_c[2]-1][0] + C[1]*coords[tri_c[2]-1][1] + C[2]*coords[tri_c[2]-1][2];
                det = C[0]*C[0] + C[1]*C[1] + C[2]*C[2];
                id = 1.0/det;
                x[0] = id*(f[0]*(B[1]*C[2] - B[2]*C[1]) + f[1]*(A[2]*C[1] - A[1]*C[2]) + f[2]*C[0]);
                x[1] = id*(f[0]*(B[2]*C[0] - B[0]*C[2]) + f[1]*(A[0]*C[2] - A[2]*C[0]) + f[2]*C[1]);
                x[2] = id*(f[0]*(B[0]*C[1] - B[1]*C[0]) + f[1]*(A[1]*C[0] - A[0]*C[1]) + f[2]*C[2]);
                double temp[3] = {x[0] - coords[tri_c[2]-1][0], x[1] - coords[tri_c[2]-1][1], x[2] - coords[tri_c[2]-1][2]};
                x[3] = temp[0]*temp[0] + temp[1]*temp[1] + temp[2]*temp[2] - weights[tri_c[2]-1];
                reftri[i].x_char[0] = x[0];
                reftri[i].x_char[1] = x[1];
                reftri[i].x_char[2] = x[2];
                reftri[i].x_char[3] = x[3];
            } else {
                for (int j =0 ; j < 3; ++ j) tri_c[j] = reftri[i].tri[j];
                for (int j = 0; j < 3; ++ j) {
                    A[j] = coords[tri_c[0]-1][j] - coords[tri_c[2]-1][j];
                    B[j] = coords[tri_c[1]-1][j] - coords[tri_c[2]-1][j];
                }
                C[0] = A[1]*B[2] - A[2]*B[1];
                C[1] = A[2]*B[0] - A[0]*B[2];
                C[2] = A[0]*B[1] - A[1]*B[0];
                f[0] = 0.5*(lift[tri_c[0]-1] - lift[tri_c[2]-1]);
                f[1] = 0.5*(lift[tri_c[1]-1] - lift[tri_c[2]-1]);
                f[2] = C[0]*coords[tri_c[2]-1][0] + C[1]*coords[tri_c[2]-1][1] + C[2]*coords[tri_c[2]-1][2];
                det = C[0]*C[0] + C[1]*C[1] + C[2]*C[2];
                if (det == 0) std::cerr << "three colinear points: stop" << std::endl;
                id = 1.0/det;
                x[0] = id*(f[0]*(B[1]*C[2] - B[2]*C[1]) + f[1]*(A[2]*C[1] - A[1]*C[2]) + f[2]*C[0]);
                x[1] = id*(f[0]*(B[2]*C[0] - B[0]*C[2]) + f[1]*(A[0]*C[2] - A[2]*C[0]) + f[2]*C[1]);
                x[2] = id*(f[0]*(B[0]*C[1] - B[1]*C[0]) + f[1]*(A[1]*C[0] - A[0]*C[1]) + f[2]*C[2]);
                double temp[3] = {x[0] - coords[tri_c[2]-1][0], x[1] - coords[tri_c[2]-1][1], x[2] - coords[tri_c[2]-1][2]};
                x[3] = temp[0]*temp[0] + temp[1]*temp[1] + temp[2]*temp[2] - weights[tri_c[2]-1];
                reftri[i].x_char[0] = x[0];
                reftri[i].x_char[1] = x[1];
                reftri[i].x_char[2] = x[2];
                reftri[i].x_char[3] = x[3];
                if (x[3] >= 0) {
                    reftri[i].alf = 1;
                } else {
                    attached = 0;
                    for (int j = 0; j < 2; ++ j) {
                        if (reftri[i].op_ver[j] != 0) {
                            double difference[3] = {x[0]-coords[reftri[i].op_ver[j]-1][0],
                                                    x[1]-coords[reftri[i].op_ver[j]-1][1],
                                                    x[2]-coords[reftri[i].op_ver[j]-1][2]};
                            dist = difference[0]*difference[0] + difference[1]*difference[1] + difference[2]*difference[2]
                                    - x[3] - weights[reftri[i].op_ver[j]-1];
                            if (dist <= 0) {
                                attached = 1;
                                break;
                            }
                        }
                    }
                    if (attached == 1) {
                        reftri[i].alf = 1;
                    } else {
                        num_tri_ext += 1;
                        current_ext->next_ext = &reftri[i];
                        current_ext = &reftri[i];
                        reftri[i].frac = 1.0;
                        reftri[i].alf = -1;
                    }
                }
            }
        }
    }
    current_ext->next_ext = NULL;
}

void extract_edg_tri(int** atri, int** aedg, int num_T, int num_coord, ver_type* verlist, tri_type* reftri, edge_type* refedg, int** T, int** T_nbrs, int& num_tri, int &num_edg)
{
    for (int i = 0; i < num_T; ++ i) {
        for (int j = 0; j < 4; ++ j) {
            atri[i][j] = 0;
            aedg[i][j] = 0;
        }
        aedg[i][4] = 0;
        aedg[i][5] = 0;
    }

    num_edg = 0;
    num_tri = 0;
    int num_adj_edg1, num_adj_edg2;

    for (int i = 0; i < num_coord; ++ i) verlist[i].num_adj_edg = 0;

    for (int i = 0; i < 4*num_T; ++ i) {
        reftri[i].adjT[0] = -5;
        reftri[i].adjT[1] = -5;
    }

    for (int i = 0; i < num_T; ++ i) {
        for (int j = 0; j < 6; ++ j) {
            if (aedg[i][j] == 0) {
                num_edg += 1;
                for (int k = 0; k < 20; ++ k) refedg[num_edg-1].tetra_ring[k] = -5;
                refedg[num_edg-1].num_ring = 0;
                if (j == 0) {
                    refedg[num_edg-1].edg[0] = T[i][0];
                    refedg[num_edg-1].edg[1] = T[i][1];
                    aedg[i][0] = num_edg;
                    num_adj_edg1 = verlist[T[i][0]-1].num_adj_edg + 1;
                    num_adj_edg2 = verlist[T[i][1]-1].num_adj_edg + 1;
                    verlist[T[i][0]-1].adj_edg[num_adj_edg1-1] = num_edg;
                    verlist[T[i][1]-1].adj_edg[num_adj_edg2-1] = num_edg;
                    verlist[T[i][0]-1].op_ver[num_adj_edg1-1] = T[i][1];
                    verlist[T[i][1]-1].op_ver[num_adj_edg2-1] = T[i][0];
                    verlist[T[i][0]-1].num_adj_edg = num_adj_edg1;
                    verlist[T[i][1]-1].num_adj_edg = num_adj_edg2;

                    tetra_ring(num_edg, 1, i+1, refedg, T, T_nbrs, atri, num_tri, reftri, aedg);
                }
                else if (j == 1) {
                    refedg[num_edg-1].edg[0] = T[i][2];
                    refedg[num_edg-1].edg[1] = T[i][0];
                    aedg[i][1] = num_edg;
                    num_adj_edg1 = verlist[T[i][0]-1].num_adj_edg + 1;
                    num_adj_edg2 = verlist[T[i][2]-1].num_adj_edg + 1;
                    verlist[T[i][0]-1].adj_edg[num_adj_edg1-1] = num_edg;
                    verlist[T[i][2]-1].adj_edg[num_adj_edg2-1] = num_edg;
                    verlist[T[i][0]-1].num_adj_edg = num_adj_edg1;
                    verlist[T[i][2]-1].num_adj_edg = num_adj_edg2;
                    verlist[T[i][0]-1].op_ver[num_adj_edg1-1] = T[i][2];
                    verlist[T[i][2]-1].op_ver[num_adj_edg2-1] = T[i][0];

                    tetra_ring(num_edg, 2, i+1, refedg, T, T_nbrs, atri, num_tri, reftri, aedg);
                }
                else if (j == 2) {
                    refedg[num_edg-1].edg[0] = T[i][0];
                    refedg[num_edg-1].edg[1] = T[i][3];
                    aedg[i][2] = num_edg;
                    num_adj_edg1 = verlist[T[i][0]-1].num_adj_edg + 1;
                    num_adj_edg2 = verlist[T[i][3]-1].num_adj_edg + 1;
                    verlist[T[i][0]-1].adj_edg[num_adj_edg1-1] = num_edg;
                    verlist[T[i][3]-1].adj_edg[num_adj_edg2-1] = num_edg;
                    verlist[T[i][0]-1].num_adj_edg = num_adj_edg1;
                    verlist[T[i][3]-1].num_adj_edg = num_adj_edg2;
                    verlist[T[i][0]-1].op_ver[num_adj_edg1-1] = T[i][3];
                    verlist[T[i][3]-1].op_ver[num_adj_edg2-1] = T[i][0];

                    tetra_ring(num_edg, 3, i+1, refedg, T, T_nbrs, atri, num_tri, reftri, aedg);
                }
                else if (j == 3) {
                    refedg[num_edg-1].edg[0] = T[i][1];
                    refedg[num_edg-1].edg[1] = T[i][2];
                    aedg[i][3] = num_edg;
                    num_adj_edg1 = verlist[T[i][1]-1].num_adj_edg + 1;
                    num_adj_edg2 = verlist[T[i][2]-1].num_adj_edg + 1;
                    verlist[T[i][1]-1].adj_edg[num_adj_edg1-1] = num_edg;
                    verlist[T[i][2]-1].adj_edg[num_adj_edg2-1] = num_edg;
                    verlist[T[i][1]-1].num_adj_edg = num_adj_edg1;
                    verlist[T[i][2]-1].num_adj_edg = num_adj_edg2;
                    verlist[T[i][1]-1].op_ver[num_adj_edg1-1] = T[i][2];
                    verlist[T[i][2]-1].op_ver[num_adj_edg2-1] = T[i][1];

                    tetra_ring(num_edg, 4, i+1, refedg, T, T_nbrs, atri, num_tri, reftri, aedg);
                }
                else if (j == 4) {
                    refedg[num_edg-1].edg[0] = T[i][3];
                    refedg[num_edg-1].edg[1] = T[i][1];
                    aedg[i][4] = num_edg;
                    num_adj_edg1 = verlist[T[i][1]-1].num_adj_edg + 1;
                    num_adj_edg2 = verlist[T[i][3]-1].num_adj_edg + 1;
                    verlist[T[i][1]-1].adj_edg[num_adj_edg1-1] = num_edg;
                    verlist[T[i][3]-1].adj_edg[num_adj_edg2-1] = num_edg;
                    verlist[T[i][1]-1].num_adj_edg = num_adj_edg1;
                    verlist[T[i][3]-1].num_adj_edg = num_adj_edg2;
                    verlist[T[i][1]-1].op_ver[num_adj_edg1-1] = T[i][3];
                    verlist[T[i][3]-1].op_ver[num_adj_edg2-1] = T[i][1];

                    tetra_ring(num_edg, 5, i+1, refedg, T, T_nbrs, atri, num_tri, reftri, aedg);
                }
                else if (j == 5) {
                    refedg[num_edg-1].edg[0] = T[i][2];
                    refedg[num_edg-1].edg[1] = T[i][3];
                    aedg[i][5] = num_edg;
                    num_adj_edg1 = verlist[T[i][2]-1].num_adj_edg + 1;
                    num_adj_edg2 = verlist[T[i][3]-1].num_adj_edg + 1;
                    verlist[T[i][2]-1].adj_edg[num_adj_edg1-1] = num_edg;
                    verlist[T[i][3]-1].adj_edg[num_adj_edg2-1] = num_edg;
                    verlist[T[i][2]-1].num_adj_edg = num_adj_edg1;
                    verlist[T[i][3]-1].num_adj_edg = num_adj_edg2;
                    verlist[T[i][2]-1].op_ver[num_adj_edg1-1] = T[i][3];
                    verlist[T[i][3]-1].op_ver[num_adj_edg2-1] = T[i][2];

                    tetra_ring(num_edg, 6, i+1, refedg, T, T_nbrs, atri, num_tri, reftri, aedg);
                }
            }
        }
    }
}

void tetra_ring(int& edg_ind, int edg_num_int, int T_ind, edge_type* refedg, int** T, int** Tnbr, int** atri, int& num_tri, tri_type* reftri, int **aedg)
{
    int current_tri, current_tetra, previous_tetra, num_in_ring, num_zero, current_edg_int, e1, e2,
            int_nums, tag, num_while, num_adj_tri, adj_tri_ind;
    num_adj_tri = 0;
    current_edg_int = edg_num_int;
    current_tetra = T_ind;
    refedg[edg_ind-1].tetra_ring[0] = T_ind;
    num_in_ring = 1;
    num_zero = 0;
    previous_tetra = -1;
    num_while = 0;
    e1 = refedg[edg_ind-1].edg[0];
    e2 = refedg[edg_ind-1].edg[1];
    while (true) {
        num_while += 1;
        if (num_while > 20) {
            std::cerr << "num_while > 20: exit" << std::endl;
            break;
        }
        tag = 0;
        if (current_edg_int == 1) {
            if (Tnbr[current_tetra-1][3] != previous_tetra) {
                current_tri = 1;
                num_adj_tri += 1;
                adj_tri_ind = atri[current_tetra-1][0];
                if (adj_tri_ind == 0) {
                   num_tri += 1;
                   reftri[num_tri-1].tri[0] = T[current_tetra-1][0];
                   reftri[num_tri-1].tri[1] = T[current_tetra-1][1];
                   reftri[num_tri-1].tri[2] = T[current_tetra-1][2];
                   reftri[num_tri-1].adjT[0] = current_tetra;
                   reftri[num_tri-1].op_ver[0] = T[current_tetra-1][3];
                   atri[current_tetra-1][0] = num_tri;
                   adj_tri_ind = num_tri;
                   tag = num_tri;
                }
                refedg[edg_ind-1].adj_tri[num_adj_tri-1] = adj_tri_ind;
                refedg[edg_ind-1].op_ver[num_adj_tri-1] = T[current_tetra-1][2];
                previous_tetra = current_tetra;
                current_tetra = Tnbr[current_tetra-1][5-current_tri-1];
                if (current_tetra != 0) {
                    mark_edg_tri(edg_ind, tag, T[previous_tetra-1][2], current_tetra, current_edg_int,
                            refedg, T, reftri, aedg, atri);
                } else {
                    if (tag != 0) {
                        reftri[num_tri-1].adjT[1] = 0;
                        reftri[num_tri-1].op_ver[1] = 0;
                    }
                }
            } else {
                current_tri = 2;
                num_adj_tri += 1;
                adj_tri_ind = atri[current_tetra-1][1];
                if (atri[current_tetra-1][1] == 0) {
                    num_tri += 1;
                    reftri[num_tri-1].tri[0] = T[current_tetra-1][0];
                    reftri[num_tri-1].tri[1] = T[current_tetra-1][1];
                    reftri[num_tri-1].tri[2] = T[current_tetra-1][3];
                    reftri[num_tri-1].adjT[0] = current_tetra;
                    reftri[num_tri-1].op_ver[0] = T[current_tetra-1][2];
                    atri[current_tetra-1][1] = num_tri;
                    tag = num_tri;
                    adj_tri_ind = num_tri;
                }
                refedg[edg_ind-1].adj_tri[num_adj_tri-1] = adj_tri_ind;
                refedg[edg_ind-1].op_ver[num_adj_tri-1] = T[current_tetra-1][3];
                previous_tetra = current_tetra;
                current_tetra = Tnbr[current_tetra-1][4-current_tri];
                if (current_tetra != 0) {
                    mark_edg_tri(edg_ind, tag, T[previous_tetra-1][3], current_tetra, current_edg_int,
                            refedg, T, reftri, aedg, atri);
                } else {
                    if (tag != 0) {
                        reftri[num_tri-1].adjT[1] = 0;
                        reftri[num_tri-1].op_ver[1] = 0;
                    }
                }
            }
        }
        else if (current_edg_int == 2) {
            if (Tnbr[current_tetra-1][3] != previous_tetra) {
                current_tri = 1;
                num_adj_tri += 1;
                adj_tri_ind = atri[current_tetra-1][0];
                if (adj_tri_ind == 0) {
                   num_tri += 1;
                   reftri[num_tri-1].tri[0] = T[current_tetra-1][0];
                   reftri[num_tri-1].tri[1] = T[current_tetra-1][1];
                   reftri[num_tri-1].tri[2] = T[current_tetra-1][2];
                   reftri[num_tri-1].adjT[0] = current_tetra;
                   reftri[num_tri-1].op_ver[0] = T[current_tetra-1][3];
                   atri[current_tetra-1][0] = num_tri;
                   adj_tri_ind = num_tri;
                   tag = num_tri;
                }
                refedg[edg_ind-1].adj_tri[num_adj_tri-1] = adj_tri_ind;
                refedg[edg_ind-1].op_ver[num_adj_tri-1] = T[current_tetra-1][1];
                previous_tetra = current_tetra;
                current_tetra = Tnbr[current_tetra-1][4-current_tri];
                if (current_tetra != 0) {
                    mark_edg_tri(edg_ind, tag, T[previous_tetra-1][1], current_tetra, current_edg_int,
                            refedg, T, reftri, aedg, atri);
                } else {
                    if (tag != 0) {
                        reftri[num_tri-1].adjT[1] = 0;
                        reftri[num_tri-1].op_ver[1] = 0;
                    }
                }
            } else {
                current_tri = 3;
                num_adj_tri += 1;
                adj_tri_ind = atri[current_tetra-1][2];
                if (atri[current_tetra-1][2] == 0) {
                    num_tri += 1;
                    reftri[num_tri-1].tri[0] = T[current_tetra-1][0];
                    reftri[num_tri-1].tri[1] = T[current_tetra-1][2];
                    reftri[num_tri-1].tri[2] = T[current_tetra-1][3];
                    reftri[num_tri-1].adjT[0] = current_tetra;
                    reftri[num_tri-1].op_ver[0] = T[current_tetra-1][1];
                    atri[current_tetra-1][2] = num_tri;
                    tag = num_tri;
                    adj_tri_ind = num_tri;
                }
                refedg[edg_ind-1].adj_tri[num_adj_tri-1] = adj_tri_ind;
                refedg[edg_ind-1].op_ver[num_adj_tri-1] = T[current_tetra-1][3];
                previous_tetra = current_tetra;
                current_tetra = Tnbr[current_tetra-1][4-current_tri];
                if (current_tetra != 0) {
                    mark_edg_tri(edg_ind, tag, T[previous_tetra-1][3], current_tetra, current_edg_int,
                            refedg, T, reftri, aedg, atri);
                } else {
                    if (tag != 0) {
                        reftri[num_tri-1].adjT[1] = 0;
                        reftri[num_tri-1].op_ver[1] = 0;
                    }
                }
            }
        }
        else if (current_edg_int == 3) {
            if (Tnbr[current_tetra-1][2] != previous_tetra) {
                current_tri = 2;
                num_adj_tri += 1;
                adj_tri_ind = atri[current_tetra-1][1];
                if (atri[current_tetra-1][1] == 0) {
                   num_tri += 1;
                   reftri[num_tri-1].tri[0] = T[current_tetra-1][0];
                   reftri[num_tri-1].tri[1] = T[current_tetra-1][1];
                   reftri[num_tri-1].tri[2] = T[current_tetra-1][3];
                   reftri[num_tri-1].adjT[0] = current_tetra;
                   reftri[num_tri-1].op_ver[0] = T[current_tetra-1][2];
                   atri[current_tetra-1][1] = num_tri;
                   adj_tri_ind = num_tri;
                   tag = num_tri;
                }
                refedg[edg_ind-1].adj_tri[num_adj_tri-1] = adj_tri_ind;
                refedg[edg_ind-1].op_ver[num_adj_tri-1] = T[current_tetra-1][1];
                previous_tetra = current_tetra;
                current_tetra = Tnbr[current_tetra-1][4-current_tri];
                if (current_tetra != 0) {
                    mark_edg_tri(edg_ind, tag, T[previous_tetra-1][1], current_tetra, current_edg_int,
                            refedg, T, reftri, aedg, atri);
                } else {
                    if (tag != 0) {
                        reftri[num_tri-1].adjT[1] = 0;
                        reftri[num_tri-1].op_ver[1] = 0;
                    }
                }
            } else {
                current_tri = 3;
                num_adj_tri += 1;
                adj_tri_ind = atri[current_tetra-1][2];
                if (atri[current_tetra-1][2] == 0) {
                    num_tri += 1;
                    reftri[num_tri-1].tri[0] = T[current_tetra-1][0];
                    reftri[num_tri-1].tri[1] = T[current_tetra-1][2];
                    reftri[num_tri-1].tri[2] = T[current_tetra-1][3];
                    reftri[num_tri-1].adjT[0] = current_tetra;
                    reftri[num_tri-1].op_ver[0] = T[current_tetra-1][1];
                    atri[current_tetra-1][2] = num_tri;
                    tag = num_tri;
                    adj_tri_ind = num_tri;
                }
                refedg[edg_ind-1].adj_tri[num_adj_tri-1] = adj_tri_ind;
                refedg[edg_ind-1].op_ver[num_adj_tri-1] = T[current_tetra-1][2];
                previous_tetra = current_tetra;
                current_tetra = Tnbr[current_tetra-1][4-current_tri];
                if (current_tetra != 0) {
                    mark_edg_tri(edg_ind, tag, T[previous_tetra-1][2], current_tetra, current_edg_int,
                            refedg, T, reftri, aedg, atri);
                } else {
                    if (tag != 0) {
                        reftri[num_tri-1].adjT[1] = 0;
                        reftri[num_tri-1].op_ver[1] = 0;
                    }
                }
            }
        }
        else if (current_edg_int == 4) {
            if (Tnbr[current_tetra-1][3] != previous_tetra) {
                current_tri = 1;
                num_adj_tri += 1;
                adj_tri_ind = atri[current_tetra-1][0];
                if (atri[current_tetra-1][0] == 0) {
                   num_tri += 1;
                   reftri[num_tri-1].tri[0] = T[current_tetra-1][0];
                   reftri[num_tri-1].tri[1] = T[current_tetra-1][1];
                   reftri[num_tri-1].tri[2] = T[current_tetra-1][2];
                   reftri[num_tri-1].adjT[0] = current_tetra;
                   reftri[num_tri-1].op_ver[0] = T[current_tetra-1][3];
                   atri[current_tetra-1][0] = num_tri;
                   adj_tri_ind = num_tri;
                   tag = num_tri;
                }
                refedg[edg_ind-1].adj_tri[num_adj_tri-1] = adj_tri_ind;
                refedg[edg_ind-1].op_ver[num_adj_tri-1] = T[current_tetra-1][0];
                previous_tetra = current_tetra;
                current_tetra = Tnbr[current_tetra-1][4-current_tri];
                if (current_tetra != 0) {
                    mark_edg_tri(edg_ind, tag, T[previous_tetra-1][0], current_tetra, current_edg_int,
                            refedg, T, reftri, aedg, atri);
                } else {
                    if (tag != 0) {
                        reftri[num_tri-1].adjT[1] = 0;
                        reftri[num_tri-1].op_ver[1] = 0;
                    }
                }
            } else {
                current_tri = 4;
                num_adj_tri += 1;
                adj_tri_ind = atri[current_tetra-1][3];
                if (atri[current_tetra-1][3] == 0) {
                    num_tri += 1;
                    reftri[num_tri-1].tri[0] = T[current_tetra-1][1];
                    reftri[num_tri-1].tri[1] = T[current_tetra-1][2];
                    reftri[num_tri-1].tri[2] = T[current_tetra-1][3];
                    reftri[num_tri-1].adjT[0] = current_tetra;
                    reftri[num_tri-1].op_ver[0] = T[current_tetra-1][0];
                    atri[current_tetra-1][3] = num_tri;
                    tag = num_tri;
                    adj_tri_ind = num_tri;
                }
                refedg[edg_ind-1].adj_tri[num_adj_tri-1] = adj_tri_ind;
                refedg[edg_ind-1].op_ver[num_adj_tri-1] = T[current_tetra-1][3];
                previous_tetra = current_tetra;
                current_tetra = Tnbr[current_tetra-1][4-current_tri];
                if (current_tetra != 0) {
                    mark_edg_tri(edg_ind, tag, T[previous_tetra-1][3], current_tetra, current_edg_int,
                            refedg, T, reftri, aedg, atri);
                } else {
                    if (tag != 0) {
                        reftri[num_tri-1].adjT[1] = 0;
                        reftri[num_tri-1].op_ver[1] = 0;
                    }
                }
            }
        }
        else if (current_edg_int == 5) {
            if (Tnbr[current_tetra-1][2] != previous_tetra) {
                current_tri = 2;
                num_adj_tri += 1;
                adj_tri_ind = atri[current_tetra-1][1];
                if (atri[current_tetra-1][1] == 0) {
                   num_tri += 1;
                   reftri[num_tri-1].tri[0] = T[current_tetra-1][0];
                   reftri[num_tri-1].tri[1] = T[current_tetra-1][1];
                   reftri[num_tri-1].tri[2] = T[current_tetra-1][3];
                   reftri[num_tri-1].adjT[0] = current_tetra;
                   reftri[num_tri-1].op_ver[0] = T[current_tetra-1][2];
                   atri[current_tetra-1][1] = num_tri;
                   adj_tri_ind = num_tri;
                   tag = num_tri;
                }
                refedg[edg_ind-1].adj_tri[num_adj_tri-1] = adj_tri_ind;
                refedg[edg_ind-1].op_ver[num_adj_tri-1] = T[current_tetra-1][0];
                previous_tetra = current_tetra;
                current_tetra = Tnbr[current_tetra-1][4-current_tri];
                if (current_tetra != 0) {
                    mark_edg_tri(edg_ind, tag, T[previous_tetra-1][0], current_tetra, current_edg_int,
                            refedg, T, reftri, aedg, atri);
                } else {
                    if (tag != 0) {
                        reftri[num_tri-1].adjT[1] = 0;
                        reftri[num_tri-1].op_ver[1] = 0;
                    }
                }
            } else {
                current_tri = 4;
                num_adj_tri += 1;
                adj_tri_ind = atri[current_tetra-1][3];
                if (atri[current_tetra-1][3] == 0) {
                    num_tri += 1;
                    reftri[num_tri-1].tri[0] = T[current_tetra-1][1];
                    reftri[num_tri-1].tri[1] = T[current_tetra-1][2];
                    reftri[num_tri-1].tri[2] = T[current_tetra-1][3];
                    reftri[num_tri-1].adjT[0] = current_tetra;
                    reftri[num_tri-1].op_ver[0] = T[current_tetra-1][0];
                    atri[current_tetra-1][3] = num_tri;
                    tag = num_tri;
                    adj_tri_ind = num_tri;
                }
                refedg[edg_ind-1].adj_tri[num_adj_tri-1] = adj_tri_ind;
                refedg[edg_ind-1].op_ver[num_adj_tri-1] = T[current_tetra-1][2];
                previous_tetra = current_tetra;
                current_tetra = Tnbr[current_tetra-1][4-current_tri];
                if (current_tetra != 0) {
                    mark_edg_tri(edg_ind, tag, T[previous_tetra-1][2], current_tetra, current_edg_int,
                            refedg, T, reftri, aedg, atri);
                } else {
                    if (tag != 0) {
                        reftri[num_tri-1].adjT[1] = 0;
                        reftri[num_tri-1].op_ver[1] = 0;
                    }
                }
            }
        }
        else {
            if (Tnbr[current_tetra-1][1] != previous_tetra) {
                current_tri = 3;
                num_adj_tri += 1;
                adj_tri_ind = atri[current_tetra-1][2];
                if (atri[current_tetra-1][2] == 0) {
                   num_tri += 1;
                   reftri[num_tri-1].tri[0] = T[current_tetra-1][0];
                   reftri[num_tri-1].tri[1] = T[current_tetra-1][2];
                   reftri[num_tri-1].tri[2] = T[current_tetra-1][3];
                   reftri[num_tri-1].adjT[0] = current_tetra;
                   reftri[num_tri-1].op_ver[0] = T[current_tetra-1][1];
                   atri[current_tetra-1][2] = num_tri;
                   adj_tri_ind = num_tri;
                   tag = num_tri;
                }
                refedg[edg_ind-1].adj_tri[num_adj_tri-1] = adj_tri_ind;
                refedg[edg_ind-1].op_ver[num_adj_tri-1] = T[current_tetra-1][0];
                previous_tetra = current_tetra;
                current_tetra = Tnbr[current_tetra-1][4-current_tri];
                if (current_tetra != 0) {
                    mark_edg_tri(edg_ind, tag, T[previous_tetra-1][0], current_tetra, current_edg_int,
                            refedg, T, reftri, aedg, atri);
                } else {
                    if (tag != 0) {
                        reftri[num_tri-1].adjT[1] = 0;
                        reftri[num_tri-1].op_ver[1] = 0;
                    }
                }
            } else {
                current_tri = 4;
                num_adj_tri += 1;
                adj_tri_ind = atri[current_tetra-1][3];
                if (atri[current_tetra-1][3] == 0) {
                    num_tri += 1;
                    reftri[num_tri-1].tri[0] = T[current_tetra-1][1];
                    reftri[num_tri-1].tri[1] = T[current_tetra-1][2];
                    reftri[num_tri-1].tri[2] = T[current_tetra-1][3];
                    reftri[num_tri-1].adjT[0] = current_tetra;
                    reftri[num_tri-1].op_ver[0] = T[current_tetra-1][0];
                    atri[current_tetra-1][3] = num_tri;
                    tag = num_tri;
                    adj_tri_ind = num_tri;
                }
                refedg[edg_ind-1].adj_tri[num_adj_tri-1] = adj_tri_ind;
                refedg[edg_ind-1].op_ver[num_adj_tri-1] = T[current_tetra-1][1];
                previous_tetra = current_tetra;
                current_tetra = Tnbr[current_tetra-1][4-current_tri];
                if (current_tetra != 0) {
                    mark_edg_tri(edg_ind, tag, T[previous_tetra-1][1], current_tetra, current_edg_int,
                            refedg, T, reftri, aedg, atri);
                } else {
                    if (tag != 0) {
                        reftri[num_tri-1].adjT[1] = 0;
                        reftri[num_tri-1].op_ver[1] = 0;
                    }
                }
            }
        }

        if (current_tetra != T_ind) {
            num_in_ring += 1;
            refedg[edg_ind-1].tetra_ring[num_in_ring-1] = current_tetra;
        } else {
            if (reftri[num_tri-1].adjT[1] == -5)
                reftri[num_tri-1].adjT[1] = T_ind;
            break;
        }
        if (current_tetra == 0) {
            if (num_zero == 0) {
                num_zero = 1;
                current_tetra = T_ind;
                previous_tetra = refedg[edg_ind-1].tetra_ring[1];
                current_edg_int = edg_num_int;
            }
            else if (num_zero == 1) {
                break;
            }
        }
    }
    refedg[edg_ind-1].num_ring = num_in_ring;
    refedg[edg_ind-1].num_adj_tri = num_adj_tri;
}

void mark_edg_tri(int& edg_ind, int& tag, int& o_ver, int& tetra_ind, int& edg_int, edge_type* refedg, int** T, tri_type* reftri, int** aedg, int** atri)
{
    int e1, e2, tri_int;
    e1 = refedg[edg_ind-1].edg[0];
    e2 = refedg[edg_ind-1].edg[1];
    int tetra[4] = {T[tetra_ind-1][0], T[tetra_ind-1][1], T[tetra_ind-1][2], T[tetra_ind-1][3] };
    if (tag == 0) {
        if (e1 == tetra[0] || e2 == tetra[0]) {
            if (e1 == tetra[1] || e2 == tetra[1])
                edg_int = 1;
            else if (e1 == tetra[2] || e2 == tetra[2])
                edg_int = 2;
            else if (e1 == tetra[3] || e2 == tetra[3])
                edg_int = 3;
        }
        else if (e1 == tetra[1] || e2 == tetra[1]) {
            if (e1 == tetra[2] || e2 == tetra[2])
                edg_int = 4;
            else if (e1 == tetra[3] || e2 == tetra[3])
                 edg_int = 5;
        }
        else if (e1 == tetra[2] || e2 == tetra[2]) {
            if (e1 == tetra[3] || e2 == tetra[3])
                edg_int = 6;
        }
        aedg[tetra_ind-1][edg_int-1] = edg_ind;
    }
    else {
        reftri[tag-1].adjT[1] = tetra_ind;
        if (e1 == tetra[0] || e2 == tetra[0]) {
            if (e1 == tetra[1] || e2 == tetra[1]) {
                edg_int = 1;
                if (o_ver == tetra[2]) {
                    tri_int = 1;
                    reftri[tag-1].op_ver[1] = tetra[3];
                } else {
                    tri_int = 2;
                    reftri[tag-1].op_ver[1] = tetra[2];
                }
            }
            else if (e1 == tetra[2] || e2 == tetra[2]) {
                edg_int = 2;
                if (o_ver == tetra[1]) {
                    tri_int = 1;
                    reftri[tag-1].op_ver[1] = tetra[3];
                } else {
                    tri_int = 3;
                    reftri[tag-1].op_ver[1] = tetra[1];
                }
            }
            else if (e1 == tetra[3] || e2 == tetra[3]) {
                edg_int = 3;
                if (o_ver == tetra[1]) {
                    tri_int = 2;
                    reftri[tag-1].op_ver[1] = tetra[2];
                } else {
                    tri_int = 3;
                    reftri[tag-1].op_ver[1] = tetra[1];
                }
            }
        }
        else if (e1 == tetra[1] || e2 == tetra[1]) {
            if (e1 == tetra[2] || e2 == tetra[2]) {
                edg_int = 4;
                if (o_ver == tetra[0]) {
                    tri_int = 1;
                    reftri[tag-1].op_ver[1] = tetra[3];
                } else {
                    tri_int = 4;
                    reftri[tag-1].op_ver[1] = tetra[0];
                }
            }
            else if (e1 == tetra[3] || e2 == tetra[3]) {
                edg_int = 5;
                if (o_ver == tetra[0]) {
                    tri_int = 2;
                    reftri[tag-1].op_ver[1]= tetra[2];
                } else {
                    tri_int = 4;
                    reftri[tag-1].op_ver[1] = tetra[0];
                }
            }
        }
        else if (e1 == tetra[2] || e2 == tetra[2]) {
            if (e1 == tetra[3] || e2 == tetra[3]) {
                edg_int = 6;
                if (o_ver == tetra[0]) {
                    tri_int = 3;
                    reftri[tag-1].op_ver[1] = tetra[1];
                } else {
                    tri_int = 4;
                    reftri[tag-1].op_ver[1] = tetra[0];
                }
            }
        }
        aedg[tetra_ind-1][edg_int-1] = edg_ind;
        atri[tetra_ind-1][tri_int-1] = tag;
    }
}

void classify_T(double** coords, double* lift, int& num_T, int** T, ver_type* verlist,
                int& maxupver, double** char_pt_T, int* alfT, int** T2, int& num_T2)
{
    double x[4], temp;
    num_T2 = 0;
    int index, pt;
    for(int i = 0; i < num_T; ++ i) {
        double p1[4] = {coords[T[i][0]-1 ][0], coords[T[i][0]-1 ][1], coords[T[i][0]-1 ][2], lift[T[i][0]-1 ]};
        double p2[4] = {coords[T[i][1]-1 ][0], coords[T[i][1]-1 ][1], coords[T[i][1]-1 ][2], lift[T[i][1]-1 ]};
        double p3[4] = {coords[T[i][2]-1 ][0], coords[T[i][2]-1 ][1], coords[T[i][2]-1 ][2], lift[T[i][2]-1 ]};
        double p4[4] = {coords[T[i][3]-1 ][0], coords[T[i][3]-1 ][1], coords[T[i][3]-1 ][2], lift[T[i][3]-1 ]};
        double M[4][4] = { {2.0*p1[0], 2.0*p1[1], 2.0*p1[2], -1.0},
                           {2.0*p2[0], 2.0*p2[1], 2.0*p2[2], -1.0},
                           {2.0*p3[0], 2.0*p3[1], 2.0*p3[2], -1.0},
                           {2.0*p4[0], 2.0*p4[1], 2.0*p4[2], -1.0} };
        double b[4] = {p1[3], p2[3], p3[3], p4[3]};
        operation::Minvb(x, M, b);
        char_pt_T[i][0] = x[0];
        char_pt_T[i][1] = x[1];
        char_pt_T[i][2] = x[2];
        temp = x[0]*x[0]+x[1]*x[1]+x[2]*x[2]-x[3];
        if (temp < 0) {
            alfT[i+1] = -1;
            num_T2 += 1;
            for (int j = 0; j < 4; ++ j) {
                T2[num_T2-1][j] = T[i][j];
                pt = T[i][j];
                index = verlist[pt-1].last_upT+1;
                if (index > maxupver) {
                    std::cerr << "increase maxupver. stop.";
                    exit(0);
                }
                verlist[pt-1].upverT[index-1] = i+1;
                verlist[pt-1].last_upT = index;
            }
        }
        else if (temp > 0) {
            alfT[i+1] = 1;
        }
        else alfT[i+1] = 0;
    }
}

#endif // FIND_VOL_MOD_H
