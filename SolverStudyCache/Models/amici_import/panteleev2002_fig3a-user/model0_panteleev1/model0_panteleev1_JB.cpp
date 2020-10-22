#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_panteleev1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 + 1.0*dwdx1;
    JB[3] = 1.0*dwdx1;
    JB[4] = -1.0*dwdx1;
    JB[6] = 1.0*dwdx0;
    JB[7] = -1.0*dwdx0;
    JB[10] = 1.0*dwdx2 - 1.0*dwdx3 + 1.0*dwdx4;
    JB[11] = -1.0*dwdx2;
    JB[12] = 1.0*dwdx3;
    JB[14] = 1.0*dwdx2;
    JB[15] = -1.0*dwdx3;
    JB[16] = 1.0*dwdx4;
    JB[17] = -1.0*dwdx4;
    JB[19] = 1.0*dwdx5;
    JB[20] = -1.0*dwdx5 + 1.0*dwdx6 + 1.0*dwdx7;
    JB[21] = -1.0*dwdx6;
    JB[22] = -1.0*dwdx7;
    JB[23] = 1.0*dwdx5 - 1.0*dwdx7;
    JB[25] = 1.0*dwdx7;
    JB[27] = 1.0*dwdx9;
    JB[28] = -1.0*dwdx8;
    JB[30] = 1.0*dwdx8 + 1.0*dwdx9;
    JB[31] = -1.0*dwdx9;
    JB[33] = -1.0*dwdx8;
    JB[36] = 1.0*dwdx10;
    JB[38] = 1.0*dwdx11;
    JB[39] = 1.0*dwdx10;
    JB[40] = -1.0*dwdx10 - 1.0*dwdx11 + 1.0*dwdx12;
    JB[41] = -1.0*dwdx11;
    JB[43] = 1.0*dwdx11;
    JB[44] = -1.0*dwdx12;
    JB[46] = 1.0*dwdx13;
    JB[47] = -1.0*dwdx13 + 1.0*dwdx14;
    JB[49] = -1.0*dwdx14;
    JB[50] = 1.0*dwdx13 - 1.0*dwdx14;
    JB[52] = 1.0*dwdx14;
    JB[54] = 1.0*dwdx16;
    JB[55] = -1.0*dwdx15;
    JB[57] = 1.0*dwdx15;
    JB[60] = -1.0*dwdx15 + 1.0*dwdx16;
    JB[61] = -1.0*dwdx16;
    JB[63] = 1.0*dwdx17;
    JB[64] = 1.0*dwdx18;
    JB[65] = 1.0*dwdx19;
    JB[67] = -1.0*dwdx19;
    JB[68] = -1.0*dwdx19;
    JB[69] = 1.0*dwdx17;
    JB[70] = -1.0*dwdx17 + 1.0*dwdx18 + 1.0*dwdx19;
    JB[71] = -1.0*dwdx18;
    JB[73] = 1.0*dwdx20;
    JB[76] = 1.0*dwdx21;
    JB[79] = 1.0*dwdx20;
    JB[80] = -1.0*dwdx20 - 1.0*dwdx21;
}