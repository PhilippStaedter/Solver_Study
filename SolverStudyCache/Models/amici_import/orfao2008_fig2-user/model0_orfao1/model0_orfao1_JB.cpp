#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_orfao1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 + 1.0*dwdx1;
    JB[1] = -1.0*dwdx0 - 1.0*dwdx1;
    JB[13] = 1.0*dwdx3 + 1.0*dwdx4;
    JB[14] = -1.0*dwdx4;
    JB[15] = -1.0*dwdx3;
    JB[19] = 1.0*dwdx2;
    JB[20] = -1.0*dwdx2;
    JB[52] = 1.0*dwdx5;
    JB[53] = -1.0*dwdx5;
    JB[56] = 1.0*dwdx5;
    JB[58] = 1.0*dwdx5;
    JB[60] = 1.0*dwdx7;
    JB[61] = -1.0*dwdx7;
    JB[64] = 1.0*dwdx6;
    JB[65] = -1.0*dwdx6;
    JB[68] = 1.0*dwdx6;
    JB[70] = 1.0*dwdx6;
    JB[81] = 1.0*dwdx8;
    JB[82] = -1.0*dwdx8;
    JB[91] = 1.0*dwdx9;
    JB[92] = -1.0*dwdx9;
    JB[100] = 1.0*dwdx10;
    JB[101] = -1.0*dwdx10;
    JB[104] = 1.0*dwdx10;
    JB[106] = 1.0*dwdx10;
    JB[117] = 1.0*dwdx11;
    JB[118] = -1.0*dwdx11;
    JB[120] = 1.0*dwdx14;
    JB[121] = -1.0*dwdx14;
    JB[124] = 1.0*dwdx13;
    JB[125] = -1.0*dwdx13;
    JB[128] = 1.0*dwdx13;
    JB[130] = 1.0*dwdx12 + 1.0*dwdx13;
    JB[131] = -1.0*dwdx12;
}