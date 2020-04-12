#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_lou1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0;
    JB[1] = -1.0*dwdx1;
    JB[2] = -1.0*dwdx2;
    JB[4] = 1.0*dwdx1;
    JB[5] = 1.0*dwdx2;
    JB[6] = -1.0*dwdx3;
    JB[7] = 1.0*dwdx4;
    JB[8] = -1.0*dwdx5;
    JB[9] = 1.0*dwdx3;
    JB[11] = 1.0*dwdx5;
    JB[12] = -1.0*dwdx6;
    JB[13] = -1.0*dwdx7;
    JB[14] = -1.0*dwdx8 + 1.0*dwdx9;
    JB[15] = 1.0*dwdx6;
    JB[16] = 1.0*dwdx7;
    JB[17] = 1.0*dwdx8;
    JB[18] = -1.0*dwdx10;
    JB[19] = -1.0*dwdx12;
    JB[20] = -1.0*dwdx13;
    JB[21] = 1.0*dwdx10 + 1.0*dwdx11;
    JB[22] = 1.0*dwdx12;
    JB[23] = 1.0*dwdx13;
    JB[24] = -1.0*dwdx14;
    JB[25] = -1.0*dwdx15;
    JB[26] = -1.0*dwdx17;
    JB[27] = 1.0*dwdx14;
    JB[28] = 1.0*dwdx15 + 1.0*dwdx16;
    JB[29] = 1.0*dwdx17;
    JB[30] = -1.0*dwdx18;
    JB[31] = -1.0*dwdx19;
    JB[32] = -1.0*dwdx20;
    JB[33] = 1.0*dwdx18;
    JB[34] = 1.0*dwdx19;
    JB[35] = 1.0*dwdx20 + 1.0*dwdx21;
}