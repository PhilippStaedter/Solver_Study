#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_Yang2007(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -1.0*K15*k15*lin*x15*(1 + x11/KI22 + x13/KI21 + x2/KI20 + x4/KI19)/(ks*pow(k15*(1 + x1/ks) + lin, 2));
    dwdx[1] = -1.0*K16*x1*x16/pow(k16*(1 + x2/ks) + x1, 2) + 1.0*K16*x16/(k16*(1 + x2/ks) + x1);
    dwdx[2] = -1.0*K17*x1*x17/pow(k17*(1 + x4/ks + x4/ki18 + x3/ki16) + x1, 2) + 1.0*K17*x17/(k17*(1 + x4/ks + x4/ki18 + x3/ki16) + x1);
    dwdx[3] = -1.0*K18*x1*x18/pow(k18*(1 + x6/ks + x7/ki3) + x1, 2) + 1.0*K18*x18/(k18*(1 + x6/ks + x7/ki3) + x1);
    dwdx[4] = -1.0*K19*k19*x19*x6/(ki1*pow(k19*(1 + x7/ks + x3/ki2 + x1/ki1) + x6, 2));
    dwdx[5] = -1.0*K21*x1*x21/pow(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1, 2) + 1.0*K21*x21/(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1);
    dwdx[6] = 0.10000000000000001;
    dwdx[7] = -1.0*K21*k21*x1*x21/(ks*pow(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1, 2));
    dwdx[8] = -1.0*K24*x10*x24/pow(k24*(1 + x11/ks) + x10, 2) + 1.0*K24*x24/(k24*(1 + x11/ks) + x10);
    dwdx[9] = -1.0*K21*x10*x21/pow(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10, 2) + 1.0*K21*x21/(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10);
    dwdx[10] = 1.0*ki10*x21;
    dwdx[11] = 1.0*K15*lin*x15/(KI22*(k15*(1 + x1/ks) + lin));
    dwdx[12] = -1.0*K21*k21*x1*x21/(ki12*pow(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1, 2));
    dwdx[13] = -1.0*K24*k24*x10*x24/(ks*pow(k24*(1 + x11/ks) + x10, 2));
    dwdx[14] = -1.0*K21*k21*x10*x21/(ki12*pow(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10, 2));
    dwdx[15] = -1.0*K23*k23*x13*x23/(ki15*pow(k23*(1 + x14/ks + x11/ki15 + x5/ki14) + x13, 2));
    dwdx[16] = 1.0*kd11;
    dwdx[17] = -1.0*K21*k21*x10*x21/(ks*pow(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10, 2));
    dwdx[18] = -1.0*K22*x12*x22/pow(k22*(1 + x13/ks) + x12, 2) + 1.0*K22*x22/(k22*(1 + x13/ks) + x12);
    dwdx[19] = -129.0*K22*x12*x22/pow(129*k22 + 129*x12, 2) + 1.0*K22*x22/(129*k22 + 129*x12);
    dwdx[20] = 1.0*kd12;
    dwdx[21] = 1.0*ki9*x21;
    dwdx[22] = 1.0*K15*lin*x15/(KI21*(k15*(1 + x1/ks) + lin));
    dwdx[23] = -1.0*K22*k22*x12*x22/(ks*pow(k22*(1 + x13/ks) + x12, 2));
    dwdx[24] = -1.0*K23*x13*x23/pow(k23*(1 + x14/ks + x11/ki15 + x5/ki14) + x13, 2) + 1.0*K23*x23/(k23*(1 + x14/ks + x11/ki15 + x5/ki14) + x13);
    dwdx[25] = 1.0*KI23*x21;
    dwdx[26] = 1.0*kd13;
    dwdx[27] = -1.0*K23*k23*x13*x23/(ks*pow(k23*(1 + x14/ks + x11/ki15 + x5/ki14) + x13, 2));
    dwdx[28] = 1.0*K15*lin*(1 + x11/KI22 + x13/KI21 + x2/KI20 + x4/KI19)/(k15*(1 + x1/ks) + lin);
    dwdx[29] = 1.0*K16*x1/(k16*(1 + x2/ks) + x1);
    dwdx[30] = 1.0*kd16;
    dwdx[31] = 1.0*K17*x1/(k17*(1 + x4/ks + x4/ki18 + x3/ki16) + x1);
    dwdx[32] = 1.0*ki17*x2;
    dwdx[33] = 1.0*K18*x1/(k18*(1 + x6/ks + x7/ki3) + x1);
    dwdx[34] = 1.0*K19*x6/(k19*(1 + x7/ks + x3/ki2 + x1/ki1) + x6);
    dwdx[35] = 1.0*K15*lin*x15/(KI20*(k15*(1 + x1/ks) + lin));
    dwdx[36] = -1.0*K16*k16*x1*x16/(ks*pow(k16*(1 + x2/ks) + x1, 2));
    dwdx[37] = -1.0*K24*x2*x24/pow(k24*(1 + x3/ks) + x2, 2) + 1.0*K24*x24/(k24*(1 + x3/ks) + x2);
    dwdx[38] = 1.0*ki17*x17;
    dwdx[39] = 1.0*ki4*x20;
    dwdx[40] = 1.0*kd2;
    dwdx[41] = 1.0*ki6*x21;
    dwdx[42] = 1.0*K20*x6/(k20*(1 + x8/ks) + x6);
    dwdx[43] = 1.0*ki4*x2;
    dwdx[44] = 1.0*ki5*x6;
    dwdx[45] = 1.0*K21*x1/(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1);
    dwdx[46] = 1.0*K21*x10/(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10);
    dwdx[47] = 1.0*KI23*x13;
    dwdx[48] = 1.0*ki9*x12;
    dwdx[49] = 1.0*ki10*x10;
    dwdx[50] = 1.0*ki6*x2;
    dwdx[51] = 1.0*K22*x12/(k22*(1 + x13/ks) + x12);
    dwdx[52] = 1.0*K22*x12/(129*k22 + 129*x12);
    dwdx[53] = 1.0*K23*x13/(k23*(1 + x14/ks + x11/ki15 + x5/ki14) + x13);
    dwdx[54] = 1.0*K24*x2/(k24*(1 + x3/ks) + x2);
    dwdx[55] = 1.0*K24*x4/(k24*(1 + x5/ks) + x4);
    dwdx[56] = 1.0*K24*x10/(k24*(1 + x11/ks) + x10);
    dwdx[57] = -1.0*K24*k24*x2*x24/(ks*pow(k24*(1 + x3/ks) + x2, 2));
    dwdx[58] = -1.0*K17*k17*x1*x17/(ki16*pow(k17*(1 + x4/ks + x4/ki18 + x3/ki16) + x1, 2));
    dwdx[59] = -1.0*K19*k19*x19*x6/(ki2*pow(k19*(1 + x7/ks + x3/ki2 + x1/ki1) + x6, 2));
    dwdx[60] = -1.0*K21*k21*x1*x21/(ki8*pow(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1, 2));
    dwdx[61] = -1.0*K21*k21*x10*x21/(ki8*pow(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10, 2));
    dwdx[62] = 1.0*kd3;
    dwdx[63] = 1.0*K15*lin*x15/(KI19*(k15*(1 + x1/ks) + lin));
    dwdx[64] = -1.0*K17*k17*x1*x17*(1.0/ks + 1.0/ki18)/pow(k17*(1 + x4/ks + x4/ki18 + x3/ki16) + x1, 2);
    dwdx[65] = -1.0*K24*x24*x4/pow(k24*(1 + x5/ks) + x4, 2) + 1.0*K24*x24/(k24*(1 + x5/ks) + x4);
    dwdx[66] = -1.0*K24*k24*x24*x4/(ks*pow(k24*(1 + x5/ks) + x4, 2));
    dwdx[67] = -1.0*K21*k21*x1*x21/(ki7*pow(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1, 2));
    dwdx[68] = -1.0*K21*k21*x10*x21/(ki7*pow(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10, 2));
    dwdx[69] = -1.0*K23*k23*x13*x23/(ki14*pow(k23*(1 + x14/ks + x11/ki15 + x5/ki14) + x13, 2));
    dwdx[70] = -1.0*K18*k18*x1*x18/(ks*pow(k18*(1 + x6/ks + x7/ki3) + x1, 2));
    dwdx[71] = -1.0*K19*x19*x6/pow(k19*(1 + x7/ks + x3/ki2 + x1/ki1) + x6, 2) + 1.0*K19*x19/(k19*(1 + x7/ks + x3/ki2 + x1/ki1) + x6);
    dwdx[72] = -1.0*K20*x20*x6/pow(k20*(1 + x8/ks) + x6, 2) + 1.0*K20*x20/(k20*(1 + x8/ks) + x6);
    dwdx[73] = 1.0*ki5*x20;
    dwdx[74] = -1.0*K18*k18*x1*x18/(ki3*pow(k18*(1 + x6/ks + x7/ki3) + x1, 2));
    dwdx[75] = -1.0*K19*k19*x19*x6/(ks*pow(k19*(1 + x7/ks + x3/ki2 + x1/ki1) + x6, 2));
    dwdx[76] = -1.0*K21*k21*x1*x21/(ki11*pow(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1, 2));
    dwdx[77] = -1.0*K21*k21*x10*x21/(ki11*pow(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10, 2));
    dwdx[78] = -2.0*a24*pow(x7, 3)/pow(pow(KI24, 2) + pow(x7, 2), 2) + 2.0*a24*x7/(pow(KI24, 2) + pow(x7, 2));
    dwdx[79] = -1.0*K20*k20*x20*x6/(ks*pow(k20*(1 + x8/ks) + x6, 2));
    dwdx[80] = 1.0*kd8;
    dwdx[81] = 1.0*kd8;
    dwdx[82] = 1.0*kd9;
}