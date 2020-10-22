#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_tiago2(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 + 1.0*dwdx1 + 1.0*dwdx10 + 1.0*dwdx11 + 1.0*dwdx12 + 1.0*dwdx2 + 1.0*dwdx3 + 1.0*dwdx4 + 1.0*dwdx5 + 1.0*dwdx6 + 1.0*dwdx7 + 1.0*dwdx8 + 1.0*dwdx9;
    JB[2] = -1.0*dwdx4;
    JB[3] = -1.0*dwdx5;
    JB[4] = -1.0*dwdx6;
    JB[5] = -1.0*dwdx7;
    JB[6] = -1.0*dwdx8;
    JB[7] = -1.0*dwdx9;
    JB[8] = -1.0*dwdx10;
    JB[9] = -1.0*dwdx0;
    JB[12] = -1.0*dwdx11;
    JB[13] = -1.0*dwdx12;
    JB[14] = -1.0*dwdx1;
    JB[15] = -1.0*dwdx3;
    JB[16] = -1.0*dwdx2;
    JB[34] = -1.0*dwdx13;
    JB[36] = 1.0*dwdx13;
    JB[51] = -1.0*dwdx14;
    JB[54] = 1.0*dwdx14;
    JB[68] = -1.0*dwdx15;
    JB[72] = 1.0*dwdx15;
    JB[85] = -1.0*dwdx16;
    JB[90] = 1.0*dwdx16;
    JB[103] = -1.0*dwdx17;
    JB[108] = 1.0*dwdx17;
    JB[119] = -1.0*dwdx18;
    JB[126] = 1.0*dwdx18;
    JB[136] = -1.0*dwdx19;
    JB[144] = 1.0*dwdx19;
    JB[162] = 1.0*dwdx20 + 1.0*dwdx21;
    JB[163] = -1.0*dwdx20;
    JB[164] = -1.0*dwdx21;
    JB[180] = 1.0*dwdx22;
    JB[181] = -1.0*dwdx22;
    JB[187] = -1.0*dwdx23;
    JB[198] = 1.0*dwdx23;
    JB[204] = -1.0*dwdx24;
    JB[216] = 1.0*dwdx24;
    JB[221] = -1.0*dwdx25;
    JB[234] = 1.0*dwdx25;
    JB[238] = -1.0*dwdx26;
    JB[252] = 1.0*dwdx26;
    JB[256] = -1.0*dwdx27;
    JB[270] = 1.0*dwdx27;
    JB[273] = -1.0*dwdx28;
    JB[288] = 1.0*dwdx28;
}