#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_Ueda2001(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[14] = 1.0*dwdx0 - 1.0*dwdx1 + 1.0*dwdx2 + 1.0*dwdx3;
    JB[15] = -1.0*dwdx0;
    JB[16] = 1.0*dwdx1;
    JB[27] = -1.0*dwdx7;
    JB[28] = 1.0*dwdx7 + 1.0*dwdx8 + 1.0*dwdx9;
    JB[30] = -1.0*dwdx6;
    JB[32] = -1.0*dwdx4;
    JB[36] = -1.0*dwdx5;
    JB[40] = -1.0*dwdx10;
    JB[42] = 1.0*dwdx10 + 1.0*dwdx11 + 1.0*dwdx12;
    JB[55] = -1.0*dwdx14;
    JB[56] = 1.0*dwdx13 + 1.0*dwdx15;
    JB[70] = 1.0*dwdx16 + 1.0*dwdx17 + 1.0*dwdx18;
    JB[72] = -1.0*dwdx16;
    JB[74] = 1.0*dwdx16;
    JB[83] = -1.0*dwdx20;
    JB[84] = 1.0*dwdx19 + 1.0*dwdx21;
    JB[96] = 1.0*dwdx23;
    JB[98] = 1.0*dwdx22 - 1.0*dwdx23 + 1.0*dwdx24 + 1.0*dwdx25;
    JB[99] = -1.0*dwdx22;
    JB[100] = 1.0*dwdx23;
    JB[108] = -1.0*dwdx28;
    JB[110] = -1.0*dwdx26;
    JB[111] = -1.0*dwdx29;
    JB[112] = 1.0*dwdx29 + 1.0*dwdx30 + 1.0*dwdx31;
    JB[114] = -1.0*dwdx27;
    JB[122] = 1.0*dwdx32;
    JB[124] = -1.0*dwdx32;
    JB[126] = 1.0*dwdx32 + 1.0*dwdx33 + 1.0*dwdx34;
    JB[139] = -1.0*dwdx36;
    JB[140] = 1.0*dwdx35 + 1.0*dwdx37;
    JB[144] = -1.0*dwdx38;
    JB[146] = 1.0*dwdx38;
    JB[161] = 1.0*dwdx39;
}