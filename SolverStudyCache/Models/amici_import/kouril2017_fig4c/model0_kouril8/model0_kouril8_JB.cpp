#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_kouril8(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = -1.0*dwdx0 + 1.0*dwdx1;
    JB[1] = 1.0*dwdx0 - 1.0*dwdx1;
    JB[2] = -1.0*dwdx0;
    JB[4] = 1.0*dwdx0;
    JB[9] = 1.0*dwdx1;
    JB[11] = -1.0*dwdx1;
    JB[13] = -1.0*dwdx2;
    JB[14] = 1.0*dwdx2;
    JB[15] = -1.0*dwdx2;
    JB[17] = 1.0*dwdx2;
    JB[26] = -1.0*dwdx3;
    JB[27] = 1.0*dwdx3;
    JB[28] = -1.0*dwdx3 + 1.0*dwdx4 + 1.0*dwdx5;
    JB[30] = 1.0*dwdx3 - 1.0*dwdx4;
    JB[31] = -1.0*dwdx5;
    JB[33] = -1.0*dwdx5;
    JB[34] = 1.0*dwdx5;
    JB[52] = -1.0*dwdx6;
    JB[53] = 1.0*dwdx6;
    JB[54] = -1.0*dwdx6;
    JB[56] = 1.0*dwdx6;
    JB[67] = 1.0*dwdx7;
    JB[70] = -1.0*dwdx7 + 1.0*dwdx8;
    JB[72] = -1.0*dwdx7;
    JB[73] = 1.0*dwdx7;
    JB[84] = 1.0*dwdx9;
    JB[85] = 1.0*dwdx9;
    JB[86] = -1.0*dwdx9;
    JB[93] = 1.0*dwdx10;
    JB[96] = -1.0*dwdx10;
    JB[97] = 1.0*dwdx11;
    JB[98] = -1.0*dwdx10 + 1.0*dwdx11;
    JB[99] = 1.0*dwdx10 - 1.0*dwdx11;
    JB[106] = 1.0*dwdx12;
    JB[109] = -1.0*dwdx12;
    JB[111] = -1.0*dwdx12;
    JB[112] = 1.0*dwdx12;
    JB[117] = 1.0*dwdx13;
    JB[118] = -1.0*dwdx13;
    JB[126] = 1.0*dwdx13;
    JB[128] = -1.0*dwdx13;
    JB[132] = 1.0*dwdx14;
    JB[135] = -1.0*dwdx14;
    JB[137] = -1.0*dwdx14;
    JB[138] = 1.0*dwdx14;
}