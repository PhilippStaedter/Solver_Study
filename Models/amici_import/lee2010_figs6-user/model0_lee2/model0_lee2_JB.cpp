#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_lee2(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 + 1.0*dwdx1 + 1.0*dwdx2 + 1.0*dwdx3;
    JB[1] = -1.0*dwdx2;
    JB[2] = -1.0*dwdx1;
    JB[3] = -1.0*dwdx0;
    JB[4] = -1.0*dwdx3;
    JB[5] = 1.0*dwdx2;
    JB[6] = 1.0*dwdx0 + 1.0*dwdx3;
    JB[7] = 1.0*dwdx1;
    JB[9] = 1.0*dwdx4 - 1.0*dwdx5;
    JB[10] = -1.0*dwdx4 + 1.0*dwdx5;
    JB[14] = 1.0*dwdx4;
    JB[17] = -1.0*dwdx5;
    JB[18] = 1.0*dwdx6 - 1.0*dwdx7;
    JB[20] = -1.0*dwdx6 + 1.0*dwdx7;
    JB[25] = 1.0*dwdx6;
    JB[26] = -1.0*dwdx7;
    JB[27] = 1.0*dwdx8 - 1.0*dwdx9;
    JB[30] = -1.0*dwdx8 + 1.0*dwdx9;
    JB[32] = -1.0*dwdx9;
    JB[33] = 1.0*dwdx8;
    JB[36] = 1.0*dwdx10 - 1.0*dwdx11;
    JB[40] = -1.0*dwdx10 + 1.0*dwdx11;
    JB[42] = 1.0*dwdx10;
    JB[43] = -1.0*dwdx11;
    JB[45] = 1.0*dwdx12;
    JB[46] = -1.0*dwdx12;
    JB[50] = 1.0*dwdx12;
    JB[54] = 1.0*dwdx13 + 1.0*dwdx14;
    JB[57] = -1.0*dwdx13;
    JB[58] = -1.0*dwdx14;
    JB[60] = 1.0*dwdx13 + 1.0*dwdx14;
    JB[63] = 1.0*dwdx15;
    JB[65] = -1.0*dwdx15;
    JB[70] = 1.0*dwdx15;
}