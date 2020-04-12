#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_kolodkin3(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[1] = -1.0*dwdx0;
    JB[2] = -1.0*dwdx1;
    JB[5] = 1.0*dwdx1;
    JB[11] = -1.0*dwdx2 + 1.0*dwdx3 + 3.4444444444400002*dwdx4;
    JB[13] = -1.0*dwdx4;
    JB[14] = -1.0*dwdx3;
    JB[16] = 1.0*dwdx4;
    JB[17] = 1.0*dwdx3;
    JB[22] = 1.0*dwdx5 - 1.0*dwdx6;
    JB[23] = -1.0*dwdx5;
    JB[25] = 1.0*dwdx6;
    JB[31] = 3.4444444444400002*dwdx7;
    JB[32] = 1.0*dwdx8;
    JB[33] = -1.0*dwdx7 - 1.0*dwdx8;
    JB[36] = 1.0*dwdx7;
    JB[41] = 1.0*dwdx10;
    JB[44] = -1.0*dwdx10 + 1.0*dwdx9;
    JB[47] = 1.0*dwdx10;
    JB[48] = 1.0*dwdx9;
    JB[49] = -1.0*dwdx9;
    JB[52] = -1.0*dwdx12;
    JB[55] = 1.0*dwdx11 + 1.0*dwdx12;
    JB[56] = -1.0*dwdx11;
    JB[61] = 3.4444444444400002*dwdx13;
    JB[63] = -1.0*dwdx13;
    JB[65] = 1.0*dwdx14;
    JB[66] = 1.0*dwdx13 - 1.0*dwdx14;
    JB[71] = 1.0*dwdx15;
    JB[74] = -1.0*dwdx15;
    JB[77] = 1.0*dwdx15;
    JB[84] = 1.0*dwdx16;
    JB[88] = 1.0*dwdx16;
    JB[89] = -1.0*dwdx16;
    JB[94] = 1.0*dwdx17;
    JB[98] = 1.0*dwdx17;
    JB[99] = -1.0*dwdx17;
}