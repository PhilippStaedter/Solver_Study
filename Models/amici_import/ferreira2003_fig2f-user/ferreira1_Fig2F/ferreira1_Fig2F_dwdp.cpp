#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_ferreira1_Fig2F(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*Glucose*Lysine*v1a_p1;
            break;
        case 1:
            dwdp[0] = 1.0*Glucose*Lysine*v1a_k1a;
            break;
        case 2:
            dwdp[1] = 1.0*Schiff;
            break;
        case 3:
            dwdp[2] = 1.0*Schiff*v2a_p2;
            break;
        case 4:
            dwdp[2] = 1.0*Schiff*v2a_k2a;
            break;
        case 5:
            dwdp[3] = 1.0*Amadori*v2b_p2;
            break;
        case 6:
            dwdp[3] = 1.0*Amadori*v2b_k2b;
            break;
        case 7:
            dwdp[4] = 1.6471820345351462*pow(Glucose, 0.35999999999999999)*v3_ox*v3_p3;
            break;
        case 8:
            dwdp[4] = 1.6471820345351462*pow(Glucose, 0.35999999999999999)*v3_k3*v3_ox;
            break;
        case 9:
            dwdp[4] = 1.6471820345351462*pow(Glucose, 0.35999999999999999)*v3_k3*v3_p3;
            break;
        case 10:
            dwdp[5] = 1.0*Amadori*v4_ox*v4_p4;
            break;
        case 11:
            dwdp[5] = 1.0*Amadori*v4_k4*v4_ox;
            break;
        case 12:
            dwdp[5] = 1.0*Amadori*v4_k4*v4_p4;
            break;
        case 13:
            dwdp[6] = 1.0*Glyoxal*Lysine*v5_ox*v5_p5;
            break;
        case 14:
            dwdp[6] = 1.0*Glyoxal*Lysine*v5_k5*v5_ox;
            break;
        case 15:
            dwdp[6] = 1.0*Glyoxal*Lysine*v5_k5*v5_p5;
            break;
        case 16:
            dwdp[7] = 1.0*Glyoxal;
            break;
        case 17:
            dwdp[8] = 1.6471820345351462*pow(Schiff, 0.35999999999999999)*v6_ox*v6_p6;
            break;
        case 18:
            dwdp[8] = 1.6471820345351462*pow(Schiff, 0.35999999999999999)*v6_k3*v6_ox;
            break;
        case 19:
            dwdp[8] = 1.6471820345351462*pow(Schiff, 0.35999999999999999)*v6_k3*v6_p6;
            break;
        case 20:
            dwdp[9] = 0.082359101726757311*pow(Schiff, 0.35999999999999999)*v7a_ox*v7a_p7;
            break;
        case 21:
            dwdp[9] = 0.082359101726757311*pow(Schiff, 0.35999999999999999)*v7a_k3*v7a_ox;
            break;
        case 22:
            dwdp[9] = 0.082359101726757311*pow(Schiff, 0.35999999999999999)*v7a_k3*v7a_p7;
            break;
        case 23:
            dwdp[10] = 0.0082359101726757304*pow(Schiff, 0.35999999999999999)*v7b_ox*v7b_p7;
            break;
        case 24:
            dwdp[10] = 0.0082359101726757304*pow(Schiff, 0.35999999999999999)*v7b_k3*v7b_ox;
            break;
        case 25:
            dwdp[10] = 0.0082359101726757304*pow(Schiff, 0.35999999999999999)*v7b_k3*v7b_p7;
            break;
        case 26:
            dwdp[11] = 1.6471820345351462*pow(Schiff, 0.35999999999999999)*v7c_ox*v7c_p7;
            break;
        case 27:
            dwdp[11] = 1.6471820345351462*pow(Schiff, 0.35999999999999999)*v7c_k3*v7c_ox;
            break;
        case 28:
            dwdp[11] = 1.6471820345351462*pow(Schiff, 0.35999999999999999)*v7c_k3*v7c_p7;
            break;
    }
}