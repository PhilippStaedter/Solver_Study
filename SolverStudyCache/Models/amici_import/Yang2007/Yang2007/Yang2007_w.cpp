#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_Yang2007(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*K15*lin*x15*(1 + x11/KI22 + x13/KI21 + x2/KI20 + x4/KI19)/(k15*(1 + x1/ks) + lin);
    w[1] = 1.0*K16*x1*x16/(k16*(1 + x2/ks) + x1);
    w[2] = 1.0*K24*x2*x24/(k24*(1 + x3/ks) + x2);
    w[3] = 1.0*K17*x1*x17/(k17*(1 + x4/ks + x4/ki18 + x3/ki16) + x1);
    w[4] = 1.0*K24*x24*x4/(k24*(1 + x5/ks) + x4);
    w[5] = 1.0*K18*x1*x18/(k18*(1 + x6/ks + x7/ki3) + x1);
    w[6] = 1.0*K19*x19*x6/(k19*(1 + x7/ks + x3/ki2 + x1/ki1) + x6);
    w[7] = 1.0*K20*x20*x6/(k20*(1 + x8/ks) + x6);
    w[8] = 1.0*kd8*x8;
    w[9] = 1.0*K21*x1*x21/(k21*(1 + x10/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x1);
    w[10] = 1.0*K24*x10*x24/(k24*(1 + x11/ks) + x10);
    w[11] = 1.0*K21*x10*x21/(k21*(1 + x12/ks + x3/ki8 + x5/ki7 + x11/ki12 + x7/ki11) + x10);
    w[12] = 1.0*K22*x12*x22/(k22*(1 + x13/ks) + x12);
    w[13] = 1.0*K23*x13*x23/(k23*(1 + x14/ks + x11/ki15 + x5/ki14) + x13);
    w[14] = 1.0*a24*pow(x7, 2)/(pow(KI24, 2) + pow(x7, 2));
    w[15] = 1.0*ki17*x17*x2;
    w[16] = 1.0*ki4*x2*x20;
    w[17] = 1.0*KI23*x13*x21;
    w[18] = 1.0*K22*x12*x22/(129*k22 + 129*x12);
    w[19] = 1.0*kd9*x9;
    w[20] = 1.0*ki5*x20*x6;
    w[21] = 1.0*kd8*x8;
    w[22] = 1.0*kd13*x13;
    w[23] = 1.0*kd12*x12;
    w[24] = 1.0*kd3*x3;
    w[25] = 1.0*kd2*x2;
    w[26] = 1.0*kd16*x16;
    w[27] = 1.0*kd11*x11;
    w[28] = 1.0*ki9*x12*x21;
    w[29] = 1.0*ki10*x10*x21;
    w[30] = 1.0*ki6*x2*x21;
    w[31] = 0.10000000000000001*x1;
}