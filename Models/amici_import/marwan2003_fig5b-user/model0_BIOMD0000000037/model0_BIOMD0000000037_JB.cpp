#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_BIOMD0000000037(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[8] = 1.0*dwdx0;
    JB[9] = -1.0*dwdx0;
    JB[13] = 1.0*dwdx1;
    JB[15] = -1.0*dwdx1;
    JB[37] = -1.0*dwdx3;
    JB[38] = -1.0*dwdx2;
    JB[39] = 1.0*dwdx2 + 1.0*dwdx3;
    JB[42] = -1.0*dwdx4;
    JB[43] = 1.0*dwdx4;
    JB[52] = 1.0*dwdx5;
    JB[53] = -1.0*dwdx6;
    JB[64] = -1.0*dwdx7;
    JB[65] = 1.0*dwdx8;
    JB[78] = 1.0*dwdx9;
    JB[79] = -1.0*dwdx9;
    JB[82] = -1.0*dwdx10;
    JB[83] = 1.0*dwdx10;
    JB[90] = -1.0*dwdx11;
    JB[91] = 1.0*dwdx11;
    JB[100] = -1.0*dwdx13;
    JB[104] = 1.0*dwdx12;
    JB[105] = -1.0*dwdx12;
    JB[106] = 1.0*dwdx13;
    JB[124] = -1.0*dwdx14;
    JB[130] = 1.0*dwdx14;
    JB[142] = -1.0*dwdx15;
    JB[143] = 1.0*dwdx15;
}