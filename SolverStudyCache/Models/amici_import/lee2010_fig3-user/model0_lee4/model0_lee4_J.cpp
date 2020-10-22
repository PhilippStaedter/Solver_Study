#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_lee4(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5;
    J[1] = -1.0*dwdx6 + 1.0*dwdx7;
    J[2] = -1.0*dwdx8 + 1.0*dwdx9;
    J[4] = -1.0*dwdx10 + 1.0*dwdx11;
    J[5] = -1.0*dwdx12 + 1.0*dwdx13;
    J[6] = -1.0*dwdx14 + 1.0*dwdx16;
    J[7] = -1.0*dwdx18 + 1.0*dwdx19;
    J[8] = -1.0*dwdx21;
    J[9] = -1.0*dwdx22;
    J[10] = -1.0*dwdx23 - 1.0*dwdx24;
    J[11] = -1.0*dwdx26;
    J[12] = -1.0*dwdx27;
    J[14] = 1.0*dwdx4;
    J[15] = 1.0*dwdx6 - 1.0*dwdx7;
    J[22] = 1.0*dwdx21;
    J[28] = 1.0*dwdx3;
    J[30] = 1.0*dwdx8 - 1.0*dwdx9;
    J[37] = 1.0*dwdx22;
    J[56] = 1.0*dwdx2;
    J[60] = 1.0*dwdx10 - 1.0*dwdx11;
    J[67] = 1.0*dwdx26;
    J[70] = 1.0*dwdx1;
    J[75] = 1.0*dwdx12 - 1.0*dwdx13;
    J[82] = 1.0*dwdx27;
    J[84] = 1.0*dwdx0;
    J[90] = 1.0*dwdx14 - 1.0*dwdx15 - 1.0*dwdx16;
    J[94] = 1.0*dwdx23;
    J[98] = 1.0*dwdx5;
    J[105] = -1.0*dwdx17 + 1.0*dwdx18 - 1.0*dwdx19;
    J[108] = 1.0*dwdx24;
    J[112] = -1.0*dwdx4;
    J[113] = -1.0*dwdx6;
    J[118] = 1.0*dwdx16;
    J[120] = -1.0*dwdx20 - 1.0*dwdx21;
    J[126] = -1.0*dwdx3;
    J[128] = -1.0*dwdx8;
    J[134] = 1.0*dwdx20;
    J[135] = -1.0*dwdx22;
    J[140] = -1.0*dwdx0 - 1.0*dwdx5;
    J[146] = -1.0*dwdx14;
    J[147] = -1.0*dwdx18;
    J[150] = -1.0*dwdx23 - 1.0*dwdx24;
    J[154] = -1.0*dwdx2;
    J[158] = -1.0*dwdx10;
    J[161] = 1.0*dwdx19;
    J[165] = -1.0*dwdx25 - 1.0*dwdx26;
    J[168] = -1.0*dwdx1;
    J[173] = -1.0*dwdx12;
    J[179] = 1.0*dwdx25;
    J[180] = -1.0*dwdx27;
    J[183] = 1.0*dwdx7;
    J[184] = 1.0*dwdx9;
    J[186] = 1.0*dwdx11;
    J[187] = 1.0*dwdx13;
    J[188] = 1.0*dwdx15;
    J[189] = 1.0*dwdx17;
}