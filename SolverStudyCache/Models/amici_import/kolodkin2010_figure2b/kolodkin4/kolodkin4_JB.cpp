#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_kolodkin4(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[1] = -1.0*dwdx0;
    JB[3] = -1.0*dwdx1;
    JB[5] = 1.0*dwdx1;
    JB[9] = -1.0*dwdx2 + 1.0*dwdx3;
    JB[10] = 1.0*dwdx3;
    JB[12] = -1.0*dwdx3;
    JB[17] = 1.0*dwdx4;
    JB[18] = 1.0*dwdx4 - 3.4444444444444402*dwdx5;
    JB[20] = -1.0*dwdx4;
    JB[21] = 1.0*dwdx5;
    JB[27] = -1.0*dwdx6 + 1.0*dwdx7;
    JB[28] = -3.4444444444444402*dwdx7;
    JB[29] = 1.0*dwdx6;
    JB[33] = 1.0*dwdx9;
    JB[34] = 1.0*dwdx9;
    JB[35] = 1.0*dwdx10;
    JB[36] = -3.4444444444444402*dwdx10 + 1.0*dwdx8 - 1.0*dwdx9;
    JB[38] = 1.0*dwdx8;
    JB[39] = -1.0*dwdx8;
    JB[42] = -3.4444444444444402*dwdx12;
    JB[43] = -1.0*dwdx11;
    JB[45] = 1.0*dwdx11 + 1.0*dwdx12;
    JB[52] = 1.0*dwdx13;
    JB[54] = 1.0*dwdx13;
    JB[55] = -1.0*dwdx13;
    JB[60] = 1.0*dwdx14;
    JB[62] = 1.0*dwdx14;
    JB[63] = -1.0*dwdx14;
}