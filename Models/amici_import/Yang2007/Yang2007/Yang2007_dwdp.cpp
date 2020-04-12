#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Yang2007(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.0*K15*lin*x15*(1 + x11/KI22 + x13/KI21 + x2/KI20 + x4/KI19)/pow(k15*(1 + x1/ks) + lin, 2) + 1.0*K15*x15*(1 + x11/KI22 + x13/KI21 + x2/KI20 + x4/KI19)/(k15*(1 + x1/ks) + lin);
            break;
        case 1:
            dwdp[0] = 1.0*lin*x15*(1 + x11/KI22 + x13/KI21 + x2/KI20 + x4/KI19)/(k15*(1 + x1/ks) + lin);
            break;
        case 2:
            dwdp[0] = 1.0*K15*lin*x15*(-1 - x1/ks)*(1 + x11/KI22 + x13/KI21 + x2/KI20 + x4/KI19)/pow(k15*(1 + x1/ks) + lin, 2);
            break;
        case 3:
            dwdp[1] = 1.0*x1*x16/(k16*(1 + x2/ks) + x1);
            break;
        case 4:
            dwdp[1] = 1.0*K16*x1*x16*(-1 - x2/ks)/pow(k16*(1 + x2/ks) + x1, 2);
            break;
        case 5:
            dwdp[3] = 1.0*x1*x17/(k17*(1 + x4/ks + x4/ki18 + x3/ki16) + x1);
            break;
        case 6:
            dwdp[3] = 1.0*K17*x1*x17*(-1 - x4/ks - x4/ki18 - x3/ki16)/pow(k17*(1 + x4/ks + x4/ki18 + x3/ki16) + x1, 2);
            break;
        case 7:
            dwdp[5] = 1.0*x1*x18/(k18*(1 + x6/ks + x7/ki3) + x1);
            break;
        case 8:
            dwdp[5] = 1.0*K18*x1*x18*(-1 - x6/ks - x7/ki3)/pow(k18*(1 + x6/ks + x7/ki3) + x1, 2);
            break;
        case 9:
            dwdp[6] = 1.0*x19*x6/(k19*(1 + x7/ks + x3/ki2 + x1/ki1) + x6);
            break;
        case 10:
            dwdp[6] = 1.0*K19*x19*x6*(-1 - x7/ks - x3/ki2 - x1/ki1)/pow(k19*(1 + x7/ks + x3/ki2 + x1/ki1) + x6, 2);
            break;
        case 11:
            dwdp[7] = 1.0*x20*x6/(k20*(1 + x8/ks) + x6);
            break;
        case 12:
            dwdp[7] = 1.0*K20*x20*x6*(-1 - x8/ks)/pow(k20*(1 + x8/ks) + x6, 2);
            break;
        case 13:
            dwdp[9] = 1.0*x1*x21/(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1);
            dwdp[11] = 1.0*x10*x21/(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10);
            break;
        case 14:
            dwdp[9] = 1.0*K21*x1*x21*(-1 - x10/ks - x3/ki8 - x5/ki7 - x11/ki12 - x7/ki11)/pow(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1, 2);
            dwdp[11] = 1.0*K21*x10*x21*(-1 - x12/ks - x3/ki8 - x5/ki7 - x11/ki12 - x7/ki11)/pow(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10, 2);
            break;
        case 15:
            dwdp[12] = 1.0*x12*x22/(k22*(1 + x13/ks) + x12);
            dwdp[18] = 1.0*x12*x22/(129*k22 + 129*x12);
            break;
        case 16:
            dwdp[12] = 1.0*K22*x12*x22*(-1 - x13/ks)/pow(k22*(1 + x13/ks) + x12, 2);
            dwdp[18] = -129.0*K22*x12*x22/pow(129*k22 + 129*x12, 2);
            break;
        case 17:
            dwdp[13] = 1.0*x13*x23/(k23*(1 + x14/ks + x11/ki15 + x5/ki14) + x13);
            break;
        case 18:
            dwdp[13] = 1.0*K23*x13*x23*(-1 - x14/ks - x11/ki15 - x5/ki14)/pow(k23*(1 + x14/ks + x11/ki15 + x5/ki14) + x13, 2);
            break;
        case 19:
            dwdp[2] = 1.0*x2*x24/(k24*(1 + x3/ks) + x2);
            dwdp[4] = 1.0*x24*x4/(k24*(1 + x5/ks) + x4);
            dwdp[10] = 1.0*x10*x24/(k24*(1 + x11/ks) + x10);
            break;
        case 20:
            dwdp[2] = 1.0*K24*x2*x24*(-1 - x3/ks)/pow(k24*(1 + x3/ks) + x2, 2);
            dwdp[4] = 1.0*K24*x24*x4*(-1 - x5/ks)/pow(k24*(1 + x5/ks) + x4, 2);
            dwdp[10] = 1.0*K24*x10*x24*(-1 - x11/ks)/pow(k24*(1 + x11/ks) + x10, 2);
            break;
        case 21:
            dwdp[25] = 1.0*x2;
            break;
        case 22:
            dwdp[24] = 1.0*x3;
            break;
        case 23:
            dwdp[8] = 1.0*x8;
            dwdp[21] = 1.0*x8;
            break;
        case 24:
            dwdp[19] = 1.0*x9;
            break;
        case 25:
            dwdp[27] = 1.0*x11;
            break;
        case 26:
            dwdp[23] = 1.0*x12;
            break;
        case 27:
            dwdp[22] = 1.0*x13;
            break;
        case 28:
            dwdp[26] = 1.0*x16;
            break;
        case 29:
            dwdp[6] = 1.0*K19*k19*x1*x19*x6/(pow(ki1, 2)*pow(k19*(1 + x7/ks + x3/ki2 + x1/ki1) + x6, 2));
            break;
        case 30:
            dwdp[6] = 1.0*K19*k19*x19*x3*x6/(pow(ki2, 2)*pow(k19*(1 + x7/ks + x3/ki2 + x1/ki1) + x6, 2));
            break;
        case 31:
            dwdp[5] = 1.0*K18*k18*x1*x18*x7/(pow(ki3, 2)*pow(k18*(1 + x6/ks + x7/ki3) + x1, 2));
            break;
        case 32:
            dwdp[16] = 1.0*x2*x20;
            break;
        case 33:
            dwdp[20] = 1.0*x20*x6;
            break;
        case 34:
            dwdp[30] = 1.0*x2*x21;
            break;
        case 35:
            dwdp[9] = 1.0*K21*k21*x1*x21*x5/(pow(ki7, 2)*pow(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1, 2));
            dwdp[11] = 1.0*K21*k21*x10*x21*x5/(pow(ki7, 2)*pow(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10, 2));
            break;
        case 36:
            dwdp[9] = 1.0*K21*k21*x1*x21*x3/(pow(ki8, 2)*pow(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1, 2));
            dwdp[11] = 1.0*K21*k21*x10*x21*x3/(pow(ki8, 2)*pow(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10, 2));
            break;
        case 37:
            dwdp[28] = 1.0*x12*x21;
            break;
        case 38:
            dwdp[29] = 1.0*x10*x21;
            break;
        case 39:
            dwdp[9] = 1.0*K21*k21*x1*x21*x7/(pow(ki11, 2)*pow(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1, 2));
            dwdp[11] = 1.0*K21*k21*x10*x21*x7/(pow(ki11, 2)*pow(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10, 2));
            break;
        case 40:
            dwdp[9] = 1.0*K21*k21*x1*x11*x21/(pow(ki12, 2)*pow(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1, 2));
            dwdp[11] = 1.0*K21*k21*x10*x11*x21/(pow(ki12, 2)*pow(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10, 2));
            break;
        case 41:
            dwdp[13] = 1.0*K23*k23*x13*x23*x5/(pow(ki14, 2)*pow(k23*(1 + x14/ks + x11/ki15 + x5/ki14) + x13, 2));
            break;
        case 42:
            dwdp[13] = 1.0*K23*k23*x11*x13*x23/(pow(ki15, 2)*pow(k23*(1 + x14/ks + x11/ki15 + x5/ki14) + x13, 2));
            break;
        case 43:
            dwdp[3] = 1.0*K17*k17*x1*x17*x3/(pow(ki16, 2)*pow(k17*(1 + x4/ks + x4/ki18 + x3/ki16) + x1, 2));
            break;
        case 44:
            dwdp[15] = 1.0*x17*x2;
            break;
        case 45:
            dwdp[3] = 1.0*K17*k17*x1*x17*x4/(pow(ki18, 2)*pow(k17*(1 + x4/ks + x4/ki18 + x3/ki16) + x1, 2));
            break;
        case 46:
            dwdp[0] = -1.0*K15*lin*x15*x4/(pow(KI19, 2)*(k15*(1 + x1/ks) + lin));
            break;
        case 47:
            dwdp[0] = -1.0*K15*lin*x15*x2/(pow(KI20, 2)*(k15*(1 + x1/ks) + lin));
            break;
        case 48:
            dwdp[0] = -1.0*K15*lin*x13*x15/(pow(KI21, 2)*(k15*(1 + x1/ks) + lin));
            break;
        case 49:
            dwdp[0] = -1.0*K15*lin*x11*x15/(pow(KI22, 2)*(k15*(1 + x1/ks) + lin));
            break;
        case 50:
            dwdp[17] = 1.0*x13*x21;
            break;
        case 51:
            dwdp[14] = -2.0*KI24*a24*pow(x7, 2)/pow(pow(KI24, 2) + pow(x7, 2), 2);
            break;
        case 52:
            dwdp[14] = 1.0*pow(x7, 2)/(pow(KI24, 2) + pow(x7, 2));
            break;
        case 53:
            dwdp[0] = 1.0*K15*k15*lin*x1*x15*(1 + x11/KI22 + x13/KI21 + x2/KI20 + x4/KI19)/(pow(ks, 2)*pow(k15*(1 + x1/ks) + lin, 2));
            dwdp[1] = 1.0*K16*k16*x1*x16*x2/(pow(ks, 2)*pow(k16*(1 + x2/ks) + x1, 2));
            dwdp[2] = 1.0*K24*k24*x2*x24*x3/(pow(ks, 2)*pow(k24*(1 + x3/ks) + x2, 2));
            dwdp[3] = 1.0*K17*k17*x1*x17*x4/(pow(ks, 2)*pow(k17*(1 + x4/ks + x4/ki18 + x3/ki16) + x1, 2));
            dwdp[4] = 1.0*K24*k24*x24*x4*x5/(pow(ks, 2)*pow(k24*(1 + x5/ks) + x4, 2));
            dwdp[5] = 1.0*K18*k18*x1*x18*x6/(pow(ks, 2)*pow(k18*(1 + x6/ks + x7/ki3) + x1, 2));
            dwdp[6] = 1.0*K19*k19*x19*x6*x7/(pow(ks, 2)*pow(k19*(1 + x7/ks + x3/ki2 + x1/ki1) + x6, 2));
            dwdp[7] = 1.0*K20*k20*x20*x6*x8/(pow(ks, 2)*pow(k20*(1 + x8/ks) + x6, 2));
            dwdp[9] = 1.0*K21*k21*x1*x10*x21/(pow(ks, 2)*pow(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1, 2));
            dwdp[10] = 1.0*K24*k24*x10*x11*x24/(pow(ks, 2)*pow(k24*(1 + x11/ks) + x10, 2));
            dwdp[11] = 1.0*K21*k21*x10*x12*x21/(pow(ks, 2)*pow(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10, 2));
            dwdp[12] = 1.0*K22*k22*x12*x13*x22/(pow(ks, 2)*pow(k22*(1 + x13/ks) + x12, 2));
            dwdp[13] = 1.0*K23*k23*x13*x14*x23/(pow(ks, 2)*pow(k23*(1 + x14/ks + x11/ki15 + x5/ki14) + x13, 2));
            break;
    }
}