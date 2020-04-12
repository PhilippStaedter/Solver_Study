#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_martins1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[14] = -1.0*dwdx1;
    JB[15] = 1.0*dwdx0 + 1.0*dwdx1;
    JB[19] = -1.0*dwdx1;
    JB[22] = 1.0*dwdx0;
    JB[23] = -1.0*dwdx1;
    JB[25] = -1.0*dwdx0;
    JB[29] = -1.0*dwdx4;
    JB[30] = 1.0*dwdx2 + 1.0*dwdx3 + 1.0*dwdx4;
    JB[31] = -1.0*dwdx2;
    JB[32] = -1.0*dwdx3;
    JB[36] = -1.0*dwdx4;
    JB[45] = 1.0*dwdx5 + 1.0*dwdx6 + 1.0*dwdx7;
    JB[49] = -1.0*dwdx6;
    JB[50] = -1.0*dwdx5 - 1.0*dwdx6 - 1.0*dwdx7;
    JB[52] = -1.0*dwdx5;
    JB[55] = -1.0*dwdx7;
    JB[60] = 1.0*dwdx8 + 1.0*dwdx9;
    JB[62] = -1.0*dwdx8;
    JB[64] = -1.0*dwdx8 - 1.0*dwdx9;
    JB[68] = -1.0*dwdx9;
    JB[105] = 1.0*dwdx10;
    JB[111] = -1.0*dwdx10;
    JB[113] = 1.0*dwdx11;
    JB[120] = 1.0*dwdx11;
    JB[123] = -1.0*dwdx11;
    JB[147] = -1.0*dwdx12;
    JB[150] = 1.0*dwdx12;
    JB[168] = -1.0*dwdx14;
    JB[169] = -1.0*dwdx13;
    JB[180] = 1.0*dwdx13 + 1.0*dwdx14;
    JB[183] = -1.0*dwdx15;
    JB[187] = -1.0*dwdx16;
    JB[195] = 1.0*dwdx15 + 1.0*dwdx16;
}