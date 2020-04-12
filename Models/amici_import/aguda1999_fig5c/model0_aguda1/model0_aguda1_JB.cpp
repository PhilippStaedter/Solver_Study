#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_aguda1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx1 + 1.0*dwdx2 + 1.0*dwdx3;
    JB[1] = -1.0*dwdx3;
    JB[3] = -1.0*dwdx0;
    JB[5] = 1.0*dwdx2;
    JB[6] = -1.0*dwdx0;
    JB[8] = 1.0*dwdx3;
    JB[10] = 1.0*dwdx0;
    JB[12] = -1.0*dwdx5;
    JB[13] = 1.0*dwdx5;
    JB[15] = -1.0*dwdx4;
    JB[18] = -1.0*dwdx4;
    JB[20] = -1.0*dwdx5;
    JB[22] = 1.0*dwdx4;
    JB[26] = 1.0*dwdx6;
    JB[28] = -1.0*dwdx6;
    JB[32] = -1.0*dwdx6;
    JB[39] = 1.0*dwdx10 - 1.0*dwdx7 + 1.0*dwdx8;
    JB[43] = -1.0*dwdx9;
    JB[45] = 1.0*dwdx8;
    JB[46] = -1.0*dwdx8;
    JB[50] = -1.0*dwdx13;
    JB[51] = -1.0*dwdx11;
    JB[52] = 1.0*dwdx13 + 1.0*dwdx14 - 1.0*dwdx15 + 1.0*dwdx16;
    JB[54] = -1.0*dwdx11;
    JB[55] = 1.0*dwdx15 - 1.0*dwdx16;
    JB[56] = 1.0*dwdx12 + 1.0*dwdx13;
    JB[58] = 1.0*dwdx11;
    JB[60] = 1.0*dwdx17;
    JB[65] = 1.0*dwdx17 + 1.0*dwdx18;
    JB[69] = -1.0*dwdx19;
    JB[78] = 1.0*dwdx20;
    JB[81] = -1.0*dwdx20;
    JB[88] = -1.0*dwdx21;
    JB[91] = 1.0*dwdx21 + 1.0*dwdx22;
    JB[96] = 1.0*dwdx25;
    JB[97] = -1.0*dwdx25;
    JB[98] = -1.0*dwdx24;
    JB[100] = 1.0*dwdx24;
    JB[104] = 1.0*dwdx23 + 1.0*dwdx24 + 1.0*dwdx25 + 1.0*dwdx26;
    JB[111] = 1.0*dwdx27;
    JB[113] = -1.0*dwdx28;
    JB[117] = 1.0*dwdx27 + 1.0*dwdx29;
    JB[118] = -1.0*dwdx27;
    JB[123] = -1.0*dwdx30;
    JB[126] = -1.0*dwdx30;
    JB[130] = 1.0*dwdx30;
}