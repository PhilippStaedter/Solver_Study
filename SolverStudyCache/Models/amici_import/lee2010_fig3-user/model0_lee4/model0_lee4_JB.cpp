#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_lee4(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 + 1.0*dwdx1 + 1.0*dwdx2 + 1.0*dwdx3 + 1.0*dwdx4 + 1.0*dwdx5;
    JB[1] = -1.0*dwdx4;
    JB[2] = -1.0*dwdx3;
    JB[4] = -1.0*dwdx2;
    JB[5] = -1.0*dwdx1;
    JB[6] = -1.0*dwdx0;
    JB[7] = -1.0*dwdx5;
    JB[8] = 1.0*dwdx4;
    JB[9] = 1.0*dwdx3;
    JB[10] = 1.0*dwdx0 + 1.0*dwdx5;
    JB[11] = 1.0*dwdx2;
    JB[12] = 1.0*dwdx1;
    JB[14] = 1.0*dwdx6 - 1.0*dwdx7;
    JB[15] = -1.0*dwdx6 + 1.0*dwdx7;
    JB[22] = 1.0*dwdx6;
    JB[27] = -1.0*dwdx7;
    JB[28] = 1.0*dwdx8 - 1.0*dwdx9;
    JB[30] = -1.0*dwdx8 + 1.0*dwdx9;
    JB[37] = 1.0*dwdx8;
    JB[41] = -1.0*dwdx9;
    JB[56] = 1.0*dwdx10 - 1.0*dwdx11;
    JB[60] = -1.0*dwdx10 + 1.0*dwdx11;
    JB[67] = 1.0*dwdx10;
    JB[69] = -1.0*dwdx11;
    JB[70] = 1.0*dwdx12 - 1.0*dwdx13;
    JB[75] = -1.0*dwdx12 + 1.0*dwdx13;
    JB[82] = 1.0*dwdx12;
    JB[83] = -1.0*dwdx13;
    JB[84] = 1.0*dwdx14 - 1.0*dwdx16;
    JB[90] = -1.0*dwdx14 + 1.0*dwdx15 + 1.0*dwdx16;
    JB[92] = -1.0*dwdx16;
    JB[94] = 1.0*dwdx14;
    JB[97] = -1.0*dwdx15;
    JB[98] = 1.0*dwdx18 - 1.0*dwdx19;
    JB[105] = 1.0*dwdx17 - 1.0*dwdx18 + 1.0*dwdx19;
    JB[108] = 1.0*dwdx18;
    JB[109] = -1.0*dwdx19;
    JB[111] = -1.0*dwdx17;
    JB[112] = 1.0*dwdx21;
    JB[113] = -1.0*dwdx21;
    JB[120] = 1.0*dwdx20 + 1.0*dwdx21;
    JB[121] = -1.0*dwdx20;
    JB[126] = 1.0*dwdx22;
    JB[128] = -1.0*dwdx22;
    JB[135] = 1.0*dwdx22;
    JB[140] = 1.0*dwdx23 + 1.0*dwdx24;
    JB[146] = -1.0*dwdx23;
    JB[147] = -1.0*dwdx24;
    JB[150] = 1.0*dwdx23 + 1.0*dwdx24;
    JB[154] = 1.0*dwdx26;
    JB[158] = -1.0*dwdx26;
    JB[165] = 1.0*dwdx25 + 1.0*dwdx26;
    JB[166] = -1.0*dwdx25;
    JB[168] = 1.0*dwdx27;
    JB[173] = -1.0*dwdx27;
    JB[180] = 1.0*dwdx27;
}