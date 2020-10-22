#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_tiago1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6 - 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    J[2] = 1.0*dwdx13;
    J[3] = 1.0*dwdx14;
    J[4] = 1.0*dwdx15;
    J[5] = 1.0*dwdx16;
    J[7] = 1.0*dwdx18;
    J[8] = 1.0*dwdx19;
    J[11] = 1.0*dwdx23;
    J[12] = 1.0*dwdx24;
    J[13] = 1.0*dwdx25;
    J[14] = 1.0*dwdx26;
    J[23] = 1.0*dwdx17;
    J[32] = 1.0*dwdx27;
    J[33] = 1.0*dwdx28;
    J[34] = 1.0*dwdx4;
    J[36] = -1.0*dwdx13;
    J[51] = 1.0*dwdx5;
    J[54] = -1.0*dwdx14;
    J[68] = 1.0*dwdx6;
    J[72] = -1.0*dwdx15;
    J[85] = 1.0*dwdx7;
    J[90] = -1.0*dwdx16;
    J[102] = 1.0*dwdx8;
    J[108] = -1.0*dwdx17;
    J[119] = 1.0*dwdx9;
    J[126] = -1.0*dwdx18;
    J[136] = 1.0*dwdx10;
    J[144] = -1.0*dwdx19;
    J[153] = 1.0*dwdx0;
    J[162] = -1.0*dwdx20 - 1.0*dwdx21;
    J[179] = 1.0*dwdx20;
    J[180] = -1.0*dwdx22;
    J[196] = 1.0*dwdx21;
    J[197] = 1.0*dwdx22;
    J[198] = -1.0*dwdx23;
    J[204] = 1.0*dwdx11;
    J[216] = -1.0*dwdx24;
    J[221] = 1.0*dwdx12;
    J[234] = -1.0*dwdx25;
    J[238] = 1.0*dwdx1;
    J[252] = -1.0*dwdx26;
    J[255] = 1.0*dwdx3;
    J[270] = -1.0*dwdx27;
    J[272] = 1.0*dwdx2;
    J[288] = -1.0*dwdx28;
}